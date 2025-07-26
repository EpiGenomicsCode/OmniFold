import gemmi
import re
import subprocess
import logging
from pathlib import Path
from typing import Dict, Tuple

logger = logging.getLogger(__name__)

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

def build_mapping(aligned_query: str, aligned_template: str, q_start_offset: int = 0, t_start_offset: int = 0) -> Dict[int, int]:
    """
    Given two aligned sequences, build a mapping from query to template indices.
    Accounts for an offset if sequence segments were aligned instead of full sequences.
    """
    mapping = {}
    q_pos = 0
    t_pos = 0
    for q_char, t_char in zip(aligned_query, aligned_template):
        q_gap = q_char == '-'
        t_gap = t_char == '-'
        
        if not q_gap and not t_gap:
            # Add offsets to map back to original full-sequence coordinates
            mapping[q_pos + q_start_offset] = t_pos + t_start_offset
        
        if not q_gap:
            q_pos += 1
        if not t_gap:
            t_pos += 1
            
    return mapping


def template_seq_and_index(cif_path: str, chain_id: str) -> Tuple[str, Dict[int, int]]:
    st = gemmi.read_structure(cif_path)
    if not st:
        raise ValueError(f"Could not read structure from {cif_path}")
    
    # A structure can have multiple models, we'll work on the first one.
    model = st[0]
    chain = model.find_chain(chain_id)
    if not chain:
        raise ValueError(f"Chain {chain_id} not found in {cif_path}")
        
    polymer = chain.get_polymer()
    if not polymer:
        raise ValueError(f"Could not get polymer for chain {chain_id} in {cif_path}")

    seq = polymer.make_one_letter_sequence()
    
    mapping = {}
    for i, res in enumerate(polymer):
        # res.seqid is gemmi.SeqId (num, ins_code), .num is what we need
        # The API gives us 1-based index, we need 0-based for AF3
        mapping[i] = res.seqid.num - 1

    if not seq:
        raise ValueError(f"Chain {chain_id} missing or empty in {cif_path}")

    return seq, mapping

def kalign_pair(q_seq: str, t_seq: str) -> Tuple[str, str]:
    """Return (aligned_query, aligned_template) using Kalign‑3 via stdin."""
    fasta = f">q\n{q_seq}\n>t\n{t_seq}\n"

    try:
        # Read sequences from STDIN (no -i) and write CLUSTAL to STDOUT.
        # --format is aliased to -f
        aln = subprocess.check_output(
            ["kalign", "--format", "clu"],
            input=fasta,
            text=True
        )
    except subprocess.CalledProcessError as e:
        logger.error(f"Kalign failed while reading from stdin.")
        logger.error(f"Kalign stderr:\n{e.stderr}")
        raise

    # grab the sequences from CLUSTAL output
    q_aln = "".join(re.findall(r"^q\s+([A-Z\-]+)", aln, re.M))
    t_aln = "".join(re.findall(r"^t\s+([A-Z\-]+)", aln, re.M))
    if not q_aln or not t_aln:
        raise RuntimeError("Kalign returned no aligned sequences")
    return q_aln, t_aln