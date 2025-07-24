import argparse
import json
import logging
import os
import sys
import pickle
import tempfile
import requests
import gemmi
import pandas as pd
from pathlib import Path

# --- Start of aggressive path correction ---
import sys
import os
# Find the correct site-packages for the current conda env
conda_env_path = os.environ.get("CONDA_PREFIX")
if conda_env_path:
    site_packages = next(Path(conda_env_path).glob("lib/python*/site-packages"), None)
    if site_packages and str(site_packages) not in sys.path:
        sys.path.insert(0, str(site_packages))
# --- End of aggressive path correction ---

# Add project root to sys.path to allow imports from omnifold
project_root = Path(__file__).resolve().parent.parent.parent
sys.path.append(str(project_root))

from omnifold.util.template_aligner import template_seq_and_index, kalign_pair, build_mapping
from omnifold.util.template_export import TemplateExport

logger = logging.getLogger(__name__)


def get_query_sequences(query_file: Path) -> dict[str, str]:
    """Reads an AF3-style JSON and returns a chain_id -> sequence map."""
    with open(query_file, 'r') as f:
        query_data = json.load(f)
    
    query_sequences = {}
    for entry in query_data.get("sequences", []):
        protein_info = entry.get("protein", {})
        if "id" in protein_info and "sequence" in protein_info:
            chain_ids = protein_info["id"]
            sequence = protein_info["sequence"]
            if isinstance(chain_ids, list):
                for chain_id in chain_ids:
                    query_sequences[chain_id] = sequence
            else:  # Handle single string ID
                query_sequences[chain_ids] = sequence
                
    return query_sequences


def process_templates(colabfold_dir: Path, af3_json_path: Path, template_store_dir: Path):
    """
    Processes ColabFold templates to create a template store.
    This function combines logic from msa_manager._process_colabfold_templates
    """
    logger.info("Starting template processing...")
    
    # Setup paths
    template_store_dir.mkdir(exist_ok=True)
    pdb_dir = template_store_dir / "pdb"
    pdb_dir.mkdir(exist_ok=True)
    hits_m8_path = template_store_dir / "hits.m8"
    mapping_pkl_path = template_store_dir / "mapping.pkl"

    m8_file_path = colabfold_dir / "all_chain_templates.m8"
    msa_map_path = colabfold_dir / "msa_map.json"

    if not m8_file_path.exists() or not msa_map_path.exists():
        logger.error(f"Required files not found in {colabfold_dir}. Cannot process templates.")
        return

    # Load inputs
    query_sequences = get_query_sequences(af3_json_path)
    
    with open(msa_map_path, 'r') as f:
        raw_msa_map = json.load(f)

    chain_id_map = {}
    for key, value in raw_msa_map.items():
        original_chain_id = key.split('|')[0]
        file_hash = Path(value).stem.split('.')[0]
        chain_id_map[file_hash] = original_chain_id
    
    logger.info(f"Built chain ID map: {chain_id_map}")

    template_hits = pd.read_csv(m8_file_path, sep='\t', header=None)
    
    # Sort by E-value (col 10) and take top 4 for each query sequence
    template_hits = template_hits.sort_values(by=10, ascending=True)
    top_templates = template_hits.groupby(0).head(4)
    logger.info(f"Filtered to top {len(top_templates)} templates.")

    all_exports = []
    
    if hits_m8_path.exists():
        hits_m8_path.unlink()

    for index, row in top_templates.iterrows():
        try:
            query_id = str(row[0])
            subject_id = str(row[1])
            pdb_id, template_chain_id = subject_id.split('_')
            
            # Definitive column indices based on MMseqs2 documentation
            # MMseqs2 default format: query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits
            # Using 1-based indexing from docs, so subtract 1 for 0-based DataFrame
            q_start = int(row[6]) - 1
            q_end = int(row[7])
            t_start = int(row[8]) - 1
            t_end = int(row[9])
            e_value = float(row[10])
            identity = float(row[2])

            original_query_chain_id = chain_id_map.get(query_id)
            if not original_query_chain_id:
                logger.warning(f"Could not find original chain ID for query hash {query_id}. Skipping.")
                continue

            full_query_sequence = query_sequences.get(original_query_chain_id)
            if not full_query_sequence:
                 logger.warning(f"Could not find query sequence for chain {original_query_chain_id}. Skipping.")
                 continue

            logger.info(f"Processing hit {pdb_id}_{template_chain_id} for query {original_query_chain_id}...")

            # --- 1. Download and extract single-chain template CIF ---
            full_cif_path = pdb_dir / f"{pdb_id.lower()}.cif"
            if not full_cif_path.exists():
                rcsb_url = f"https://files.rcsb.org/download/{pdb_id}.cif"
                logger.info(f"Downloading template from {rcsb_url}")
                try:
                    response = requests.get(rcsb_url, timeout=60)
                    response.raise_for_status()
                    
                    raw_bytes = response.content
                    trimmed_bytes = raw_bytes.lstrip()
                    if not trimmed_bytes.lower().startswith(b"data_"):
                        logger.warning(f"Downloaded content for {pdb_id} does not look like a CIF file. Skipping.")
                        continue
                    
                    clean_bytes = bytes((b if b < 0x80 else 0x3F) for b in raw_bytes)
                    with tempfile.NamedTemporaryFile(dir=full_cif_path.parent, delete=False, mode='wb') as tmp_file:
                        tmp_file.write(clean_bytes)
                        tmp_file.flush()
                        os.fsync(tmp_file.fileno())
                    os.rename(tmp_file.name, full_cif_path)
                    logger.info(f"Successfully downloaded and sanitized {full_cif_path.name}")
                    
                except requests.exceptions.RequestException as e:
                    logger.error(f"Failed to download {rcsb_url}: {e}")
                    continue

            single_chain_cif_path = pdb_dir / "single_chain" / f"{subject_id.lower()}.cif"
            single_chain_cif_path.parent.mkdir(exist_ok=True)
            
            if not single_chain_cif_path.exists():
                logger.info(f"Extracting chain {template_chain_id} from {full_cif_path}")
                try:
                    st = gemmi.read_structure(str(full_cif_path))
                    if not st: raise ValueError("Structure is empty.")
                    
                    model_to_modify = st[0]
                    chains_to_remove = [ch.name for ch in model_to_modify if ch.name != template_chain_id]

                    if len(chains_to_remove) == len(model_to_modify):
                        raise ValueError(f"Chain {template_chain_id} not found in {full_cif_path}")

                    for chain_name in chains_to_remove:
                        model_to_modify.remove_chain(chain_name)

                    doc = st.make_mmcif_document()
                    doc.write_file(str(single_chain_cif_path))
                except Exception as e:
                    logger.error(f"FAILED to extract chain for {subject_id}: {e}", exc_info=True)
                    continue

            # --- 2. Get full template sequence ---
            full_template_sequence, _ = template_seq_and_index(str(single_chain_cif_path), template_chain_id)

            # --- 3. Align SEGMENTS and build mapping ---
            # Use 0-based indices directly for slicing
            query_segment = full_query_sequence[q_start:q_end]
            template_segment = full_template_sequence[t_start:t_end]

            if not query_segment or not template_segment:
                logger.warning("Empty query or template segment, skipping.")
                continue

            q_aligned, t_aligned = kalign_pair(query_segment, template_segment)
            
            # The mapping must be built with offsets to be relative to the full sequence
            mapping = build_mapping(q_aligned, t_aligned, q_start_offset=q_start, t_start_offset=t_start)

            if not mapping:
                logger.warning(f"No alignment mapping generated for {subject_id}. Skipping.")
                continue

            # Add a quality filter
            coverage = len(mapping) / len(full_query_sequence)
            if coverage < 0.3 or identity < 0.2:
                logger.info(f"Skipping low-quality template {subject_id} (Coverage: {coverage:.2f}, Identity: {identity:.2f})")
                continue

            export = TemplateExport(
                pdb_id=pdb_id,
                chain_id=template_chain_id,
                cif_path=single_chain_cif_path,
                query_idx_to_template_idx=mapping,
                e_value=e_value,
                hit_from_chain=original_query_chain_id
            )
            all_exports.append(export)
            
            # --- 4. Append to hits.m8 for Chai ---
            new_m8_line = '\t'.join(map(str, row)) + '\n'
            with open(hits_m8_path, 'a') as f_out:
                f_out.write(new_m8_line)

        except Exception as e:
            logger.error(f"Failed processing row {row.to_list()}: {e}", exc_info=True)
            continue
            
    with open(mapping_pkl_path, "wb") as pkl_f:
        pickle.dump(all_exports, pkl_f)
    logger.info(f"Successfully created template store at {template_store_dir} with {len(all_exports)} templates.")


def inject_templates_into_config(af3_json_path: Path, template_store_dir: Path, output_json_path: Path):
    """Injects template data from the store into the AF3 config."""
    logger.info("Injecting templates into AF3 config...")
    
    mapping_pkl_path = template_store_dir / "mapping.pkl"
    if not mapping_pkl_path.exists():
        logger.error(f"mapping.pkl not found in {template_store_dir}. Cannot inject templates.")
        return

    with open(mapping_pkl_path, "rb") as f:
        template_exports = pickle.load(f)

    with open(af3_json_path, 'r') as f:
        af3_config = json.load(f)

    templates_by_chain = {}
    for export in template_exports:
        if export.hit_from_chain not in templates_by_chain:
            templates_by_chain[export.hit_from_chain] = []
        
        template_entry = {
            "mmcifPath": str(export.cif_path),
            "queryIndices": list(export.query_idx_to_template_idx.keys()),
            "templateIndices": list(export.query_idx_to_template_idx.values())
        }
        templates_by_chain[export.hit_from_chain].append(template_entry)

    for protein_entry in af3_config.get("sequences", []):
        protein_info = protein_entry.get("protein", {})
        chain_ids = protein_info.get("id", [])
        if not isinstance(chain_ids, list):
            chain_ids = [chain_ids]
        
        for chain_id in chain_ids:
            if chain_id in templates_by_chain:
                protein_info["templates"] = templates_by_chain[chain_id]
                logger.info(f"Injected {len(templates_by_chain[chain_id])} templates for chain {chain_id}.")
                break 

    with open(output_json_path, 'w') as f:
        json.dump(af3_config, f, indent=2)
    logger.info(f"Successfully wrote final AF3 config with templates to {output_json_path}")


def main():
    parser = argparse.ArgumentParser(description="Generate a template store and final AF3 JSON from ColabFold outputs.")
    parser.add_argument("--colabfold_dir", type=Path, required=True, help="Path to ColabFold output dir (containing .m8 and msa_map.json).")
    parser.add_argument("--af3_json_path", type=Path, required=True, help="Path to the initial AF3 input JSON.")
    parser.add_argument("--template_store_dir", type=Path, required=True, help="Path to output the template store (hits.m8, mapping.pkl, pdb/).")
    parser.add_argument("--output_json_path", type=Path, required=True, help="Path to write the final AF3 JSON with injected templates.")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

    process_templates(args.colabfold_dir, args.af3_json_path, args.template_store_dir)
    inject_templates_into_config(args.af3_json_path, args.template_store_dir, args.output_json_path)


if __name__ == "__main__":
    main() 