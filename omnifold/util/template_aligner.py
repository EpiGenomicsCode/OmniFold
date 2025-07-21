import gemmi
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Tuple

import gemmi
from typing import Dict, Tuple

import gemmi
from typing import Dict, Tuple

def template_seq_and_index(cif_path: str, chain_id: str) -> Tuple[str, Dict[int, int]]:
    """
    Return SEQRES-style sequence for `chain_id` and a mapping
    query_idx (0-based) → template_idx (0-based, includes unresolved residues).
    """
    doc   = gemmi.cif.read_file(cif_path)
    block = doc.sole_block()

    # Prefer the modern table; fall back to the legacy one
    prefix = '_pdbx_poly_seq_scheme.' if block.find_value('_pdbx_poly_seq_scheme.mon_id') \
             else '_entity_poly_seq.'
    tbl = gemmi.cif.Table(block, prefix,
                          ['mon_id', 'auth_asym_id', 'asym_id', 'seq_id'])

    seq      : str            = ""
    mapping  : Dict[int, int] = {}

    for row in tbl:
        asym = row['auth_asym_id'] or row['asym_id']
        if asym != chain_id:
            continue
        res1 = gemmi.to_one_letter_code(row['mon_id']) or 'X'
        seq += res1
        mapping[len(seq) - 1] = int(row['seq_id']) - 1   # AF3 expects 0-based

    if not seq:
        raise ValueError(f"{chain_id} not present in {cif_path}")

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