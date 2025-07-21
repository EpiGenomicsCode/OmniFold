from dataclasses import dataclass, field
from typing import List, Literal, Generator, Any, Optional, Dict, Tuple
import string

# Sequence alphabets (from af3/constants.py)
prot_alphabet: set[str] = set("ARNDCQEGHILKMFPSTWYV") | set("arndcqeghilkmfpstwyv")
rna_alphabet: set[str] = set("ACGU") | set("acgu") 
dna_alphabet: set[str] = set("ACGT") | set("acgt") 
# Simplified CCD alphabet check for as_entity, full CCD dict will be in af3_models.py
ccd_alphabet_check: set[str] = set([chr(x) for x in range(ord("A"), ord("Z") + 1)]) | set(f"{a}" for a in range(10))
smiles_alphabet_check: set[str] = set(r"Hh+Nn27#Oo.[]58lLfF4PpSs(=C@c3-@\\16/)09*:%BbRrIiMmZzGgQq") 

SequenceType = Literal["protein", "rna", "dna", "ligand_smiles", "ligand_ccd", "unknown"]

@dataclass
class SequenceInfo:
    original_name: str
    sequence: str
    molecule_type: SequenceType
    chain_id: str
    molecule_type_confidence: float = 1.0 

@dataclass
class JobInput:
    """Unified representation of a prediction job input."""
    name_stem: str
    sequences: List[SequenceInfo]
    output_dir: str 
    raw_input_type: Literal["fasta", "af3_json", "boltz_yaml"]
    input_msa_paths: Dict[str, str] = field(default_factory=dict) 
    constraints: Optional[List[Dict[str, Any]]] = None 
    has_msa: bool = False 
    af3_data_json: Optional[str] = None
    original_af3_config_path: Optional[str] = None
    original_boltz_config_path: Optional[str] = None
    # Can be a dict of {chain_id: path} or {"unpaired": {chain_id: path}, "paired": {chain_id: path}}
    protein_id_to_a3m_path: Dict[str, str] = field(default_factory=dict) 
    protein_id_to_pqt_path: Dict[str, str] = field(default_factory=dict)
    boltz_csv_msa_dir: Optional[str] = None
    template_store_path: Optional[str] = None
    model_seeds: Optional[List[int]] = None
    num_model_seeds_from_input: Optional[int] = 1
    bonded_atom_pairs: Optional[List[Tuple[int, int]]] = None 
    is_boltz_config: bool = False 
    is_af3_msa_config_only: bool = False 

    def __post_init__(self):
        pass 

# ID generator (from af3/cli_convert.py)
def idgen() -> Generator[str, None, None]:
    """Generate sequence ids in the same way af3 documents it (A-Z, AA-ZZ)."""
    letters = [chr(x) for x in range(ord("A"), ord("Z") + 1)]
    yield from letters
    for l1 in letters:
        for l2 in letters:
            yield f"{l2}{l1}" 

def as_entity(seq: str, chain_id: str, original_name: str) -> SequenceInfo:
    """
    Guess the type of a sequence and return a SequenceInfo object.
    Adapted from af3_extras.cli_convert.as_entity.
    """
    s_upper = seq.upper() 
    sset = set(s_upper)

    if sset.issubset(rna_alphabet):
        return SequenceInfo(original_name=original_name, sequence=s_upper, molecule_type="rna", chain_id=chain_id)
    if sset.issubset(dna_alphabet):
        return SequenceInfo(original_name=original_name, sequence=s_upper, molecule_type="dna", chain_id=chain_id)
    if sset.issubset(prot_alphabet):
        return SequenceInfo(original_name=original_name, sequence=s_upper, molecule_type="protein", chain_id=chain_id)
    
    if len(s_upper) <= 3 and sset.issubset(ccd_alphabet_check):
        return SequenceInfo(original_name=original_name, sequence=s_upper, molecule_type="ligand_ccd", chain_id=chain_id)
    
    if sset.issubset(smiles_alphabet_check): 
        return SequenceInfo(original_name=original_name, sequence=seq, molecule_type="ligand_smiles", chain_id=chain_id)

    return SequenceInfo(original_name=original_name, sequence=seq, molecule_type="unknown", chain_id=chain_id, molecule_type_confidence=0.5) 