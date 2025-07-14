#!/usr/bin/env python
"""Utility to run ColabFold (MMseqs2) once and stage the outputs for
AlphaFold-3, Chai-1, and Boltz-1.

Features
--------
1. Generates unpaired and paired MSAs via the public ColabFold API.
2. Saves chain-specific aligned Parquet MSAs (.aligned.pqt) that are
   read natively by Chai-1.
3. Optionally saves the raw A3M files so downstream tools (e.g. AF3 or
   ad-hoc scripts) can consume them directly.
4. Optionally fetches template hits (pdb70.m8) which can be provided to
   AlphaFold-3 and Chai-1.

Example
-------
python generate_colabfold_msas.py input.fasta \
       --out_dir precomputed_msas \
       --include_templates \
       --write_a3m
"""
from __future__ import annotations

import argparse
import logging
import sys
import json
import os
from pathlib import Path
from typing import List, Tuple

# Use the local copy of the ColabFold helper to avoid depending on the external chai_lab package.
from .colabfold import generate_colabfold_msas  # type: ignore

# -----------------------------------------------------------------------------
# FASTA parsing helpers – minimal implementation to avoid heavy dependencies.
# -----------------------------------------------------------------------------

def _read_fasta(path: Path) -> Tuple[List[str], List[str]]:
    """Return (headers, sequences) lists from a FASTA file.

    This tiny parser keeps only the first token of each header line and
    upper-cases sequences. It does **not** validate that residues belong to a
    particular alphabet – that responsibility is delegated to downstream code
    (ColabFold/MMseqs2 will fail if the sequence contains invalid symbols)."""
    headers: List[str] = []
    seqs: List[str] = []
    try:
        with path.open() as handle:
            current_header: str | None = None
            current_seq: List[str] = []
            for line in handle:
                line = line.rstrip()
                if not line:
                    continue  # skip empty lines
                if line.startswith(">"):
                    # Commit the previous record if present.
                    if current_header is not None:
                        headers.append(current_header)
                        seqs.append("".join(current_seq).upper())
                    current_header = line[1:].split()[0]
                    current_seq = []
                else:
                    current_seq.append(line)
            # Flush last record.
            if current_header is not None:
                headers.append(current_header)
                seqs.append("".join(current_seq).upper())
    except FileNotFoundError:
        logging.error("FASTA file %s not found.", path)
        sys.exit(1)
    except Exception as exc:
        logging.error("Failed to read FASTA %s – %s", path, exc)
        sys.exit(1)

    if not headers:
        logging.error("No sequences found in FASTA %s", path)
        sys.exit(1)
    return headers, seqs


# -----------------------------------------------------------------------------
# Main driver
# -----------------------------------------------------------------------------

def _build_arg_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Generate ColabFold MSAs once and stage for AF3/Chai/ Boltz.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    p.add_argument("fasta", type=str, help="Input FASTA file containing one or more chains.")
    p.add_argument("--out_dir", required=True, type=str, help="Destination directory for the generated MSAs.")
    p.add_argument("--msa_server_url", type=str, default="https://api.colabfold.com", help="ColabFold MMseqs2 API endpoint.")
    p.add_argument("--include_templates", action="store_true", help="Fetch template hits (.m8) in addition to MSAs.")
    p.add_argument("--write_a3m", action="store_true", help="Save raw A3M files to <out_dir>/a3ms for debugging and AF3 consumption.")
    p.add_argument("--log_level", choices=["DEBUG", "INFO", "WARNING", "ERROR"], default="INFO", help="Logging verbosity.")
    return p


def main(argv: List[str] | None = None) -> None:  # noqa: D401 – imperative mood is fine here
    args = _build_arg_parser().parse_args(argv)

    # Ensure output directory is a Path object for consistent path handling
    args.out_dir = Path(args.out_dir)
    args.out_dir.mkdir(exist_ok=True, parents=True)

    logging.basicConfig(level=getattr(logging, args.log_level), format="%(asctime)s - %(levelname)s - %(message)s")

    fasta_path = Path(args.fasta).resolve()
    out_dir = args.out_dir

    headers, sequences = _read_fasta(fasta_path)
    logging.info("Parsed %d chain(s) from %s", len(sequences), fasta_path)

    logging.info("Requesting MSAs from ColabFold – this may take a while…")
    try:
        msa_mapping = generate_colabfold_msas(
            protein_seqs=sequences,
            msa_dir=out_dir,
            msa_server_url=args.msa_server_url,
            search_templates=args.include_templates,
            write_a3m_to_msa_dir=args.write_a3m,
        )
    except Exception as exc:
        logging.error("ColabFold MSA generation failed: %s", exc, exc_info=True)
        sys.exit(1)

    # ---------------------------------------------------------------------
    # Summarise results
    # ---------------------------------------------------------------------
    logging.info("MSA generation complete. Results:")

    # Create a JSON mapping file to make outputs easier to parse downstream
    # This is the primary way clients should consume the outputs.
    manifest_path = args.out_dir / "msa_map.json"
    manifest = {}
    for header, seq in zip(headers, sequences):
        pqt_path = msa_mapping.get(seq)
        if pqt_path:
            logging.info(f"  {header} -> {pqt_path}")
            manifest[header] = str(pqt_path)
        else:
            logging.warning(f"  {header} -> No MSA generated")

    manifest_path.write_text(json.dumps(manifest, indent=2))
    logging.info(f"Wrote MSA manifest to {manifest_path}")

    if args.include_templates:
        m8_path = out_dir / "all_chain_templates.m8"
        if m8_path.exists():
            logging.info("Template hits saved to %s", m8_path)
        else:
            logging.warning("--include_templates was set but template file not found – something went wrong.")

    if args.write_a3m:
        logging.info("Raw A3M files are in %s", out_dir / "a3ms")


if __name__ == "__main__":
    main() 