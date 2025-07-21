import gemmi
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Tuple

def template_seq_and_index(cif_path: str, chain_id: str) -> Tuple[str, Dict[int, int]]:
    """
    Return (sequence, mapping) where mapping[q_idx] = template_idx (0‑based).
    template_idx counts every residue in _pdbx_poly_seq_scheme, so positions
    with missing atoms are still addressable by AlphaFold‑3.
    """
    st = gemmi.read_structure(cif_path)
    chain = st[0][chain_id]                       # first model only

    polymer = chain.get_polymer()                 # None for waters/ligands
    if polymer is None:
        raise ValueError(f"{cif_path}:{chain_id} is not a polymer")

    scheme = polymer.get_poly_seq_scheme()
    seq            : str           = ""
    seq_to_res_idx : Dict[int, int] = {}

    for i, step in enumerate(scheme):             # `step` is PolySeqScheme
        aa = gemmi.to_one_letter_code(step.mon_id) or 'X'
        seq += aa
        # seq_id is 1‑based; convert to 0‑based as AF‑3 expects
        seq_to_res_idx[i] = step.seq_id - 1

    return seq, seq_to_res_idx

def kalign_pair(q_seq: str, t_seq: str) -> Tuple[str, str]:
    """Return (aligned_query, aligned_template) using Kalign‑3."""
    with tempfile.NamedTemporaryFile("w+", suffix=".fasta", delete=False) as f:
        f.write(f">q\n{q_seq}\n>t\n{t_seq}\n")
        tmp = f.name
    try:
        aln = subprocess.check_output(
            ["kalign", "-i", tmp, "-o", "/dev/stdout", "-format", "clu"],
            text=True, stderr=subprocess.PIPE
        )
    finally:
        Path(tmp).unlink(missing_ok=True)

    # grab the sequences from CLUSTAL output
    q_aln = "".join(re.findall(r"^q\s+([A-Z\-]+)", aln, re.M))
    t_aln = "".join(re.findall(r"^t\s+([A-Z\-]+)", aln, re.M))
    if not q_aln or not t_aln:
        raise RuntimeError("Kalign did not return aligned sequences")
    return q_aln, t_aln


def build_mapping(q_aln: str, t_aln: str, t_seq_to_idx: Dict[int, int]) -> Dict[int, int]:
    """Map alignment columns back to template indices."""
    q2t   : Dict[int, int] = {}
    qi = ti = -1
    for qc, tc in zip(q_aln, t_aln):
        if qc != "-": qi += 1
        if tc != "-": ti += 1
        if qc != "-" and tc != "-":
            q2t[qi] = t_seq_to_idx[ti]    # noqa: E501
    return q2t