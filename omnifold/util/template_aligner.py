import gemmi
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Tuple

TAGSET = [  # full tags – Gemmi’s Table takes *only* absolute names
    '_pdbx_poly_seq_scheme.mon_id',
    '_pdbx_poly_seq_scheme.auth_asym_id',
    '_pdbx_poly_seq_scheme.asym_id',
    '_pdbx_poly_seq_scheme.seq_id',
]

FALLBACK = [  # older files → _entity_poly_seq
    '_entity_poly_seq.mon_id',
    '_entity_poly_seq.entity_id',
    '_entity_poly_seq.num',
]

def template_seq_and_index(cif_path: str, chain_id: str) -> Tuple[str, Dict[int, int]]:
    doc   = gemmi.cif.read_file(cif_path)
    block = doc.sole_block()

    # choose the category present in this file
    tags = TAGSET if block.find_value(TAGSET[0]) else FALLBACK
    # gemmi.cif.Table expects the column tags as separate positional arguments,
    # not as a single list. Unpack the list with *.
    table = gemmi.cif.Table(block, *tags)

    seq, mapping = "", {}
    for row in table:
        asym = row[1] or row[2]  # auth_asym_id first; fallback to asym_id/entity_id
        if asym != chain_id:
            continue
        aa = gemmi.to_one_letter_code(row[0]) or 'X'
        seq += aa
        mapping[len(seq) - 1] = int(row[-1]) - 1  # 0-based for AF3

    if not seq:
        raise ValueError(f"{chain_id} missing in {cif_path}")

    return seq, mapping

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