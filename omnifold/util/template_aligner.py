import gemmi
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Tuple

def template_seq_and_index(cif_path: str, chain_id: str) -> Tuple[str, Dict[int, int]]:
    """
    Extracts a template sequence and a 0-based index mapping from a CIF file.
    
    Args:
        cif_path: Path to the mmCIF file.
        chain_id: The chain ID to extract.

    Returns:
        A tuple containing:
        - The full template sequence.
        - A dictionary mapping 0-based sequence indices to 1-based residue indices (as in CIF).
    """
    st = gemmi.read_structure(cif_path)
    model_idx = 0  # Assuming first model
    
    # Check if the chain exists in the model
    chain_names = [c.name for c in st[model_idx]]
    if chain_id not in chain_names:
        raise ValueError(f"Chain '{chain_id}' not found in model {model_idx} of {cif_path}. Available chains: {', '.join(chain_names)}")

    chain = st[model_idx][chain_id]
    
    polymer = chain.get_polymer()
    if polymer is None:
        raise ValueError(f"Chain {chain_id} is not a polymer in {cif_path}")
    
    scheme = polymer.get_poly_seq_scheme()

    sequence = ""
    # gemmi uses 1-based seqid for CIFs. We need to map our 0-based sequence index
    # to the original residue number from the file.
    seq_to_res_idx = {}
    for i, r in enumerate(scheme):
        sequence += r.mon
        if r.seqid: # Ensure seqid is not None
            seq_to_res_idx[i] = r.seqid.num - 1 # AF3 templates are 0-indexed

    return sequence, seq_to_res_idx

def kalign_pair(q_seq: str, t_seq: str) -> Tuple[str, str]:
    """
    Aligns two sequences using Kalign and returns the aligned strings.
    """
    with tempfile.NamedTemporaryFile('w+', delete=False, suffix=".fasta") as f:
        f.write(f">query\n{q_seq}\n>template\n{t_seq}\n")
        name = f.name
    
    try:
        # Note: kalign must be in the system's PATH
        aln = subprocess.check_output(
            ["kalign", "-i", name, "-o", "/dev/stdout", "-format", "clu"],
            text=True, stderr=subprocess.PIPE
        )
        
        q_aln_lines = [line for line in aln.splitlines() if line.startswith('query')]
        t_aln_lines = [line for line in aln.splitlines() if line.startswith('template')]

        if not q_aln_lines or not t_aln_lines:
            raise RuntimeError("Kalign output did not contain aligned sequences.")

        q_aln = "".join(q_aln_lines[0].split()[1:])
        t_aln = "".join(t_aln_lines[0].split()[1:])

    finally:
        Path(name).unlink()
        
    return q_aln, t_aln

def build_mapping(q_aln: str, t_aln: str, t_seq_to_idx: Dict[int, int]) -> Dict[int, int]:
    """
    Builds a 0-based query index to 0-based template index mapping from aligned sequences.
    """
    q2t = {}
    q_idx = -1
    t_idx = -1
    for q_char, t_char in zip(q_aln, t_aln):
        if q_char != '-':
            q_idx += 1
        if t_char != '-':
            t_idx += 1
        
        if q_char != '-' and t_char != '-':
            if t_idx in t_seq_to_idx:
                q2t[q_idx] = t_seq_to_idx[t_idx]
            else:
                # This can happen if the alignment extends beyond the mapped residues
                pass
    return q2t 