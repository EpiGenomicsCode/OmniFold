from dataclasses import dataclass, field
from typing import List, Literal, Generator, Any, Optional, Dict

# Sequence alphabets (from af3/constants.py)
prot_alphabet: set[str] = set("ARNDCQEGHILKMFPSTWYV") | set("arndcqeghilkmfpstwyv")
rna_alphabet: set[str] = set("ACGU") | set("acgu") # Added lowercase for robustness
dna_alphabet: set[str] = set("ACGT") | set("acgt") # Added lowercase for robustness
# Simplified CCD alphabet check for as_entity, full CCD dict will be in af3_models.py
ccd_alphabet_check: set[str] = set([chr(x) for x in range(ord("A"), ord("Z") + 1)]) | set(f"{a}" for a in range(10))
smiles_alphabet_check: set[str] = set(r"Hh+Nn27#Oo.[]58lLfF4PpSs(=C@c3-@\16/") # Expanded common SMILES chars

SequenceType = Literal["protein", "rna", "dna", "ligand_smiles", "ligand_ccd", "unknown"]

@dataclass
class SequenceInfo:
    original_name: str
    sequence: str
    molecule_type: SequenceType
    chain_id: str
    molecule_type_confidence: float = 1.0 # For guessed types

@dataclass
class JobInput:
    name_stem: str # From input filename or user provided
    sequences: List[SequenceInfo] = field(default_factory=list)
    raw_input_type: Optional[Literal["fasta", "af3_json", "boltz_yaml"]] = None
    # Stores MSA paths found *in the input file itself*, mapping chain_id to its unpairedMsaPath
    input_msa_paths: Optional[Dict[str, str]] = field(default_factory=dict) 
    # Future fields:
    # ligand_info: Any = None # Could be part of SequenceInfo or separate
    # template_data: Any = None # Could be Dict[str, List[Template]] mapping chain_id to templates


# ID generator (from af3/cli_convert.py)
def idgen() -> Generator[str, None, None]:
    """Generate sequence ids in the same way af3 documents it (A-Z, AA-ZZ)."""
    letters = [chr(x) for x in range(ord("A"), ord("Z") + 1)]
    yield from letters
    for l1 in letters:
        for l2 in letters:
            yield f"{l2}{l1}" # Corrected to match AA, AB, ... BA, BB, ...

def as_entity(seq: str, chain_id: str, original_name: str) -> SequenceInfo:
    """
    Guess the type of a sequence and return a SequenceInfo object.
    Adapted from af3_extras.cli_convert.as_entity.
    """
    s_upper = seq.upper() # Ensure case-insensitivity for alphabet checks
    sset = set(s_upper)

    # Check order is important: DNA/RNA are subsets of protein if not careful.
    # Check for more specific nucleic acid alphabets first.
    if sset.issubset(rna_alphabet):
        return SequenceInfo(original_name=original_name, sequence=s_upper, molecule_type="rna", chain_id=chain_id)
    if sset.issubset(dna_alphabet):
        return SequenceInfo(original_name=original_name, sequence=s_upper, molecule_type="dna", chain_id=chain_id)
    if sset.issubset(prot_alphabet):
        return SequenceInfo(original_name=original_name, sequence=s_upper, molecule_type="protein", chain_id=chain_id)
    
    # Ligand checks (simplified from af3_extras for now)
    # CCD codes are typically short (1-3 chars) and alphanumeric.
    if len(s_upper) <= 3 and sset.issubset(ccd_alphabet_check):
        # This is a basic check. Real CCD validation would use the full CCD dictionary.
        return SequenceInfo(original_name=original_name, sequence=s_upper, molecule_type="ligand_ccd", chain_id=chain_id)
    
    # SMILES check is harder and less definitive without a full parser.
    # Check for common SMILES characters. This is a heuristic.
    # AlphaFold 3 example for SMILES: "CC(C)C[C@H](C(=O)O)N" for Leucine
    # A simple heuristic: if it contains characters not in protein/RNA/DNA and common in SMILES.
    if sset.issubset(smiles_alphabet_check): # Using the expanded smiles_alphabet_check
         # More robust check might involve trying to parse with RDKit if available, or more specific regex.
        return SequenceInfo(original_name=original_name, sequence=seq, molecule_type="ligand_smiles", chain_id=chain_id)

    return SequenceInfo(original_name=original_name, sequence=seq, molecule_type="unknown", chain_id=chain_id, molecule_type_confidence=0.5) 