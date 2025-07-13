#!/usr/bin/env python3
"""
Create Boltz-style CSV MSAs from AlphaFold-3 paired+unpaired A3Ms
Usage:  af3_to_boltz_csv.py  --chains A,B   --msa_root /path/msas   --out /path/csv
"""
import re, argparse, textwrap, pathlib, random, gzip, sys, os
import logging

TAX_RE = re.compile(r"(?:OX|TaxID)=(\d+)")
LOWER  = str.maketrans('', '', 'abcdefghijklmnopqrstuvwxyz')

def read_a3m(path):
    hdr, seq = None, []
    # Try to open with gzip first, then as plain text
    try:
        open_func = gzip.open if path.name.endswith(".gz") else open
        mode = "rt" if path.name.endswith(".gz") else "r" # gzip needs "rt" for text mode
        with open_func(path, mode) as f:
            for line in f:
                line=line.rstrip()
                if line.startswith(">"):
                    if hdr is not None:
                        yield hdr, "".join(seq)
                    hdr, seq = line, []
                else:
                    seq.append(line)
            if hdr:
                yield hdr, "".join(seq)
    except FileNotFoundError:
        print(f"Warning: File not found {path}, skipping.")
        return
    except Exception as e:
        print(f"Warning: Could not read {path} (error: {e}), skipping.")
        return

def extract_taxid(header):
    m = TAX_RE.search(header)
    return m.group(1) if m else None

def main(opts):
    chains = [c.strip() for c in opts.chains.split(",")]
    paired_by_tax = {}
    
    msa_paired_root = opts.msa_root / "paired" # This line is no longer correct if uniprot.a3m is per chain
    # msa_combined_root = opts.msa_root / "combined" # This was already effectively replaced by json_extracted_unpaired_msa_dir

    if not opts.msa_root.is_dir():
        print(f"Error: MSA root directory {opts.msa_root} does not exist.")
        return
    # The checks for msa_paired_root and msa_combined_root might need to be re-evaluated or removed
    # if not msa_paired_root.is_dir():
    #     print(f"Warning: Paired MSA directory {msa_paired_root} does not exist. Paired MSAs might be missing.")
    # if not msa_combined_root.is_dir():
    #     print(f"Warning: Combined/Unpaired MSA directory {msa_combined_root} does not exist. Unpaired MSAs might be missing.")

    for c in chains:
        # Corrected path to the paired A3M file, assuming it's named 'uniprot.a3m' inside each chain-specific directory
        chain_specific_msa_dir = opts.msa_root / f"chain_{c}"
        paired_a3m_path = chain_specific_msa_dir / "uniprot.a3m"
        
        if not paired_a3m_path.exists(): # Try .gz as a fallback
            paired_a3m_path = chain_specific_msa_dir / "uniprot.a3m.gz"
        
        if paired_a3m_path.exists():
            query_seq_header_processed = False
            for hdr, seq in read_a3m(paired_a3m_path):
                if not query_seq_header_processed: # This is the query sequence header
                    query_seq_header_processed = True
                    continue # Skip the query sequence itself from taxid processing
                tax = extract_taxid(hdr)
                if not tax: continue
                paired_by_tax.setdefault(tax, {})[c] = seq
        else:
            # Updated warning message to reflect the new path structure
            print(f"Warning: Paired A3M file for chain {c} not found at {chain_specific_msa_dir / 'uniprot.a3m'} (or .gz).")

    paired = [(tax, seqs) for tax,seqs in paired_by_tax.items()
                               if len(seqs)==len(chains)]
    if opts.shuffle_paired:
        random.shuffle(paired)
    paired = paired[:opts.max_paired]

    rows = {c: [] for c in chains}
    current_sequences_for_chain = {c: set() for c in chains}

    for k,(tax,seqs) in enumerate(paired):
        for c in chains:
            rows[c].append((k, seqs[c]))
            current_sequences_for_chain[c].add(seqs[c])

    for c in chains:
        seen = {s for _,s in rows[c]} # Sequences already added from paired MSAs
        # Read unpaired A3M from the new dedicated directory for JSON-extracted unpaired MSAs
        unpaired_a3m_path_from_json = opts.json_extracted_unpaired_msa_dir / f"msa_{c}.a3m"
        if not unpaired_a3m_path_from_json.exists():
            # Try .gz as a fallback, though typically JSON extraction won't create .gz
            unpaired_a3m_path_from_json = opts.json_extracted_unpaired_msa_dir / f"msa_{c}.a3m.gz"

        if unpaired_a3m_path_from_json.exists():
            query_seq_unpaired_header_processed = False
            for hdr, seq in read_a3m(unpaired_a3m_path_from_json):
                if not query_seq_unpaired_header_processed: 
                    query_seq_unpaired_header_processed = True
                    # The first sequence in this file is the query as per AF3 JSON structure.
                    # We skip it because the query sequence isn't part of the MSA features for Boltz in this CSV format.
                    # Boltz gets the query sequence from the main input YAML.
                    continue 
                if seq in seen: continue # Skip if already added from paired MSAs
                if len(rows[c]) >= opts.max_total: break
                rows[c].append((-1, seq))
                # current_sequences_for_chain[c].add(seq) # No longer strictly needed to track here if we only read from json_extracted for unpaired
        else:
            # This case should ideally not happen if MSAManager successfully extracted them.
            # If it does, it means the primary source for unpaired (the JSON field) was missing or extraction failed for this chain.
            logger_print(f"Warning: Unpaired A3M for chain {c} not found at {unpaired_a3m_path_from_json} (expected from JSON extraction). This chain might have fewer/no unpaired sequences in CSV.")

    opts.out.mkdir(parents=True, exist_ok=True)
    for c in chains:
        out_path = opts.out / f"{c}.csv"
        with out_path.open("w") as f:
            f.write("key,sequence\n")
            for k,seq in rows[c]:
                f.write(f"{k},{seq}\n")
        print(f"Written Boltz CSV for chain {c} to {out_path}")

# Helper function to print to stderr for script context, as logger might not be configured.
def logger_print(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

if __name__ == "__main__":
    p = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
        description=textwrap.dedent("""\
        Convert AF3 MSAs to Boltz CSV.
        Expects AF3 MSA output structure:
        <msa_root>/chain_<CHAIN_ID>/uniprot.a3m[.gz] (for paired MSAs)
        And unpaired MSAs provided via --json_extracted_unpaired_msa_dir which should contain <msa_CHAIN_ID>.a3m files.
        Output CSVs will be named <CHAIN_ID>.csv in the --out directory.
        """))
    p.add_argument("--chains", required=True, help="Comma-separated list of chain IDs (e.g., A,B). Must match <CHAIN_ID> in filenames.")
    p.add_argument("--msa_root", type=pathlib.Path, required=True, help="Root directory of AF3's 'msas' output (e.g., .../experiment_name/msas), used for finding 'chain_<CHAIN_ID>/uniprot.a3m'.")
    p.add_argument("--json_extracted_unpaired_msa_dir", type=pathlib.Path, required=True, help="Directory containing unpaired A3Ms extracted from AF3 JSON (e.g., msa_A.a3m, msa_B.a3m).")
    p.add_argument("--out", type=pathlib.Path, required=True, help="Output directory for Boltz CSV files.")
    p.add_argument("--max_paired", type=int, default=256, help="Maximum number of paired sequences to retain.")
    p.add_argument("--max_total",  type=int, default=4096, help="Maximum total number of sequences per chain (paired + unpaired).")
    p.add_argument("--shuffle_paired", action='store_true', help="Shuffle paired MSAs before truncation. Useful for variability if many taxa are common.")

    main(p.parse_args()) 

logger = logging.getLogger(__name__)

def convert_a3m_to_boltz_csv(protein_to_a3m_path: dict, output_csv_dir: str):
    """
    Converts multiple A3M files to a simplified Boltz-compatible CSV format.
    This version does not perform pairing; it assumes the input A3M is final.

    Args:
        protein_to_a3m_path: Dictionary mapping protein IDs to their A3M file paths.
        output_csv_dir: Path to the directory where output CSV files will be saved.
    """
    os.makedirs(output_csv_dir, exist_ok=True)
    for protein_id, a3m_path in protein_to_a3m_path.items():
        output_csv_path = os.path.join(output_csv_dir, f"{protein_id.split('|')[0]}.csv") # Clean header
        
        try:
            with open(output_csv_path, 'w') as csv_file:
                csv_file.write("msa_sequence,deletion_counts\n")
                
                query_seq = None
                with open(a3m_path, 'r') as a3m_file:
                    for i, line in enumerate(a3m_file):
                        if line.startswith('>'):
                            continue # Skip header
                        
                        sequence = line.strip()
                        if i == 1: # First sequence is the query
                            query_seq = sequence
                            continue

                        if not query_seq:
                            raise ValueError("Could not determine query sequence from A3M file.")

                        # Calculate deletion counts based on the query
                        deletion_counts = [0] * len(query_seq)
                        msa_seq_no_gaps = []
                        query_idx = 0
                        for msa_char in sequence:
                            if msa_char == '-':
                                if query_idx < len(deletion_counts):
                                     deletion_counts[query_idx] += 1
                            else:
                                msa_seq_no_gaps.append(msa_char)
                                query_idx += 1
                        
                        final_msa_seq = "".join(msa_seq_no_gaps)
                        final_deletions = ",".join(map(str, deletion_counts))
                        csv_file.write(f'"{final_msa_seq}","{final_deletions}"\n')
            logger.info(f"Successfully converted {a3m_path} to {output_csv_path}")

        except Exception as e:
            logger.error(f"Failed to convert {a3m_path} for protein {protein_id}: {e}", exc_info=True)
            # Create an empty file to signify failure but allow pipeline to continue
            with open(output_csv_path, 'w') as f:
                f.write("msa_sequence,deletion_counts\n") 