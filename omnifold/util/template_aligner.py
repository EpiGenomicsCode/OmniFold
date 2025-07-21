import gemmi
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Tuple

def template_seq_and_index(cif_path: str, chain_id: str) -> Tuple[str, Dict[int, int]]:
    """
    Extract full SEQRES‑style sequence for `chain_id` from `_pdbx_poly_seq_scheme`
    and build a dict mapping: query_index (0‑based) → template_index (0‑based).
    The template_index counts every residue in the scheme, even if unresolved,
    which is what AlphaFold‑3 expects for `templateIndices`.
    """
    doc   = gemmi.cif.read_file(cif_path)
    block = doc.sole_block()

    loop = block.find_loop('_pdbx_poly_seq_scheme.mon_id')
    if loop is None:
        raise ValueError(f"{cif_path} lacks _pdbx_poly_seq_scheme")

    # Column indices (mandatory in PDBx/mmCIF spec)
    tag = {name: i for i, name in enumerate(loop.tags)}
    mon_col   = tag['mon_id']
    seq_col   = tag['seq_id']
    # Auth chain ID is in auth_asym_id; fall back to asym_id if absent
    auth_col  = tag.get('auth_asym_id', tag.get('asym_id'))

    seq: str = ""
    m  : Dict[int, int] = {}

    for row in range(loop.length()):
        if loop[auth_col][row] != chain_id:
            continue
        res3   = loop[mon_col][row]
        res1   = gemmi.to_one_letter_code(res3) or 'X'
        seq_id = int(loop[seq_col][row])        # 1‑based
        seq += res1
        m[len(seq)-1] = seq_id - 1              # convert to 0‑based

    if not seq:
        raise ValueError(f"Chain {chain_id} not found in {cif_path}")

    return seq, m
    
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