import gemmi
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Tuple

def template_seq_and_index(cif_path: str, chain_id: str) -> Tuple[str, Dict[int, int]]:
    """
    Extracts the full template sequence and a mapping from sequence index to residue index from a CIF file.
    """
    st = gemmi.read_structure(cif_path)
    # Ensure chain_id is valid for this structure
    if chain_id not in [ch.name for ch in st[0]]:
        raise ValueError(f"Chain ID '{chain_id}' not found in structure {cif_path}")
        
    block = st[0][chain_id]
    scheme = block.get_poly_seq_scheme()
    seq = ''.join(r.mon for r in scheme)
    
    # gemmi uses 1-based seqid, AF3 uses 0-based templateIndices
    seq_to_idx = {i: scheme[i].seqid.num - 1 for i in range(len(scheme))}
    return seq, seq_to_idx

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