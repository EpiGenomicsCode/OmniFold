import gemmi
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Tuple

import gemmi
from typing import Dict, Tuple

def template_seq_and_index(cif_path: str, chain_id: str) -> Tuple[str, Dict[int, int]]:
    """
    Read _pdbx_poly_seq_scheme (or _entity_poly_seq as fallback) and return:
        sequence: full SEQRES string for that chain
        mapping : {query_idx (0-based) -> template_idx (0-based)}
    """
    doc   = gemmi.cif.read_file(cif_path)
    block = doc.sole_block()

    # 1) Try the authoritative table first
    loop = block.find_loop('_pdbx_poly_seq_scheme.mon_id')
    if loop is None:
        # PDB entries < 2019 use _entity_poly_seq instead
        loop = block.find_loop('_entity_poly_seq.mon_id')
    if loop is None:
        raise ValueError(f"{cif_path} has neither poly_seq_scheme nor entity_poly_seq")

    # If we got a Column, grab its parent loop
    if isinstance(loop, gemmi.cif.Column):
        loop = loop.loop

    # Convenience helpers
    get = loop.get_column

    mon_ids  = get('mon_id')           # 3-letter codes
    asym_ids = get('auth_asym_id') if loop.find_tag('auth_asym_id') \
               else get('asym_id')     # fallback
    seq_ids  = get('seq_id')           # always present

    seq = ""
    mapping: Dict[int, int] = {}

    for i in range(loop.length()):
        if asym_ids[i] != chain_id:
            continue
        res1 = gemmi.to_one_letter_code(mon_ids[i]) or 'X'
        seq += res1
        mapping[len(seq) - 1] = int(seq_ids[i]) - 1  # 0-based for AF-3

    if not seq:
        raise ValueError(f"Chain {chain_id} not found in {cif_path}")

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