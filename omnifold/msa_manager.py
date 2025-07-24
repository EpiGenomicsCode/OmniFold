import json
import os
import subprocess
import logging
import tempfile # For temporary FASTA
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, List
import sys # Added for sys.executable
import pickle
import hashlib
import gzip
import requests

from .config_generator import ConfigGenerator
from .util.definitions import JobInput, SequenceInfo
from .util.file_converters import af3_json_to_chai_fasta, job_input_to_chai_fasta
from .util.msa_utils import extract_all_protein_a3ms_from_af3_json # Import the utility
from .util.template_export import TemplateExport
from .util.template_aligner import template_seq_and_index, kalign_pair, build_mapping


logger = logging.getLogger(__name__)


class MSAManager:
    """
    Manages the generation of Multiple Sequence Alignments (MSAs) using
    either the AlphaFold 3 pipeline or a ColabFold-like (MMseqs2) pipeline.
    """

    def __init__(self, config: Dict[str, Any], output_dir: str):
        """
        Initializes the MSA_Manager.

        Args:
            config: Global configuration dictionary containing paths and settings.
            output_dir: The *base* output directory for MSA-related files for this job.
                      Intermediate files will go into a subdirectory.
        """
        self.config = config
        self.output_dir = Path(output_dir) 
        self.msa_tmp_dir = self.output_dir / "msa_intermediate_files"
        self.msa_tmp_dir.mkdir(parents=True, exist_ok=True)
        self.config_generator = ConfigGenerator()
        self.job_input: Optional[JobInput] = None 
        logger.info(f"MSA_Manager initialized. Intermediate dir: {self.msa_tmp_dir}")

    def _get_internal_name_from_af3_json(self, json_path: Path) -> Optional[str]:
        """Parses an AF3 JSON and returns the value of its 'name' field."""
        try:
            with open(json_path, 'r') as f:
                data = json.load(f)
            return data.get("name")
        except Exception as e:
            logger.error(f"Failed to parse internal name from {json_path}: {e}")
            return None

    def _run_command(self, cmd: List[str], cwd: Optional[str] = None) -> Tuple[int, str, str]:
        """Runs a shell command and returns exit code, stdout, stderr."""
        logger.info(f"Running command: {' '.join(cmd)}")
        try:
            process = subprocess.run(
                cmd, 
                capture_output=True, 
                text=True, 
                check=False, 
                cwd=cwd
            )
            if process.stdout:
                logger.info(f"Command stdout:\n{process.stdout.strip()}")

            # Always log stderr, as many tools (like our MSA script) write
            # informational logs here. Use INFO level if the command succeeds,
            # WARNING if it fails.
            if process.stderr:
                log_level = logging.WARNING if process.returncode != 0 else logging.INFO
                logger.log(log_level, f"Command stderr:\n{process.stderr.strip()}")
            
            return process.returncode, process.stdout, process.stderr
        except FileNotFoundError:
            logger.error(f"Command not found: {cmd[0]}. Ensure it's installed and in PATH.")
            return -1, "", f"Command not found: {cmd[0]}"
        except Exception as e:
             logger.error(f"Error running command {' '.join(cmd)}: {e}", exc_info=True)
             return -1, "", str(e)

    def _process_template_metadata(self, af3_msa_output_actual_dir: Path) -> Optional[str]:
        """
        Processes intermediate template metadata JSONs to create the final
        template store for Chai-1 and Boltz-2.
        """
        logger.info("Starting template processing for Chai-1 and Boltz-2.")
        template_store_dir = self.output_dir / "templates"
        template_store_dir.mkdir(parents=True, exist_ok=True)
        pdb_dir = template_store_dir / "pdb"
        pdb_dir.mkdir(parents=True, exist_ok=True)

        m8_path = template_store_dir / "hits.m8"
        mapping_pkl_path = template_store_dir / "mapping.pkl"
        
        all_exports: List[TemplateExport] = []
        processed_cifs = set()
        seen_m8_hits = set()
        
        # This function is now defined inside where it's used to have access to scope
        def sha1(seq: str) -> str:
            return hashlib.sha1(seq.encode()).hexdigest()

        try:
            with open(m8_path, "a") as m8_fh:
                # Find all intermediate template metadata files
                metadata_files = list(af3_msa_output_actual_dir.glob("msas/chain_*/*_template_metadata.json"))
                if not metadata_files:
                    logger.info("No intermediate template metadata JSONs found. Skipping template processing.")
                    return None
                
                logger.info(f"Found {len(metadata_files)} template metadata files to process.")

                for metadata_file in metadata_files:
                    with open(metadata_file, 'r') as f:
                        template_data_list = json.load(f)
                    
                    for template_data in template_data_list:
                        pdb_id = template_data['pdb_id']
                        
                        # Define paths for both .cif (Boltz) and .cif.gz (Chai)
                        cif_name_boltz = f"{pdb_id}.cif"
                        cif_path_boltz = pdb_dir / cif_name_boltz

                        # Chai expects uppercase PDB IDs in filenames
                        cif_name_chai = f"{pdb_id.upper()}.cif.gz"
                        cif_path_chai = pdb_dir / cif_name_chai

                        if cif_name_boltz not in processed_cifs:
                            mmcif_content = template_data['mmcif_string']
                            
                            # Save .cif for Boltz
                            cif_path_boltz.write_text(mmcif_content)
                            
                            # Save .cif.gz for Chai
                            with gzip.open(cif_path_chai, 'wt', encoding='utf-8') as f_gz:
                                f_gz.write(mmcif_content)

                            processed_cifs.add(cif_name_boltz)

                        # Create TemplateExport for Boltz-2, pointing to the .cif file
                        export = TemplateExport(
                            pdb_id=pdb_id,
                            chain_id=template_data['chain_id'],
                            cif_path=cif_path_boltz.resolve(),
                            query_idx_to_template_idx=template_data['query_to_template_map'],
                            e_value=template_data.get('e_value', 999.0),
                            hit_from_chain=template_data['hit_from_chain']
                        )
                        all_exports.append(export)

                        # De-duplicate hits for the m8 file
                        m8_key = (export.pdb_id, export.chain_id, export.hit_from_chain)
                        if m8_key in seen_m8_hits:
                            continue
                        seen_m8_hits.add(m8_key)

                        # Write line for Chai-1 hits.m8 file
                        query_id = export.hit_from_chain
                        subject_id = f"{pdb_id}_{template_data['chain_id']}"
                        
                        # Calculate correct 1-based coordinates
                        q_indices = sorted([int(k) for k in template_data['query_to_template_map'].keys()])
                        t_indices = sorted([v for k, v in template_data['query_to_template_map'].items()])
                        
                        q_start = q_indices[0] + 1 if q_indices else 0
                        q_end = q_indices[-1] + 1 if q_indices else 0
                        s_start = t_indices[0] + 1 if t_indices else 0
                        s_end = t_indices[-1] + 1 if t_indices else 0

                        # Simplified stats
                        ident, length, mism, gapopen = 0, 0, 0, 0
                        evalue = template_data.get('e_value', 999.0)
                        bitscore = 0

                        logger.debug(
                            f"Writing template to m8: PDB={pdb_id}_{template_data['chain_id']}, "
                            f"QueryID={query_id}, E-value={evalue}"
                        )
                        m8_fh.write(
                            f"{query_id}\t{subject_id}\t"
                            f"{ident:.1f}\t{length}\t{mism}\t{gapopen}\t"
                            f"{q_start}\t{q_end}\t{s_start}\t{s_end}\t"
                            f"{evalue:.3e}\t{bitscore}\n"
                        )
                
            # After processing all files, save the final pickle for Boltz
            # TODO: Honor user's 'force' flag when building Boltz templates from this pickle.
            with open(mapping_pkl_path, "wb") as pkl_f:
                pickle.dump(all_exports, pkl_f)
            
            logger.info(f"Successfully created template store at: {template_store_dir}")
            logger.info(f"  - {len(processed_cifs)} mmCIF files written to {pdb_dir}")
            logger.info(f"  - {len(all_exports)} hits written to {m8_path} for Chai-1")
            logger.info(f"  - {len(all_exports)} template export objects pickled to {mapping_pkl_path} for Boltz-2")

            return str(template_store_dir.resolve())

        except Exception as e:
            logger.error(f"Failed to process template metadata and create template store: {e}", exc_info=True)
            return None

    def _process_colabfold_templates(self, colabfold_output_dir: Path) -> Optional[str]:
        """
        Processes ColabFold's m8 template file to create the canonical template store.
        """
        logger.info("Starting ColabFold template processing.")
        m8_file = colabfold_output_dir / "all_chain_templates.m8"
        if not m8_file.is_file():
            logger.warning(f"ColabFold template file not found at {m8_file}. This might be expected if no templates were found. Skipping template processing.")
            return None

        template_store_dir = self.output_dir / "templates"
        template_store_dir.mkdir(parents=True, exist_ok=True)
        pdb_dir = template_store_dir / "pdb"
        pdb_dir.mkdir(parents=True, exist_ok=True)

        hits_m8_path = template_store_dir / "hits.m8"
        mapping_pkl_path = template_store_dir / "mapping.pkl"
        
        all_exports: List[TemplateExport] = []
        
        try:
            # Get query sequences to perform alignments
            query_sequences = {seq.chain_id: seq.sequence for seq in self.job_input.sequences if seq.molecule_type == 'protein'}

            # The m8 file from colabfold has SHA1s. We need to map them back to chain IDs for Chai.
            msa_map_path = colabfold_output_dir / "msa_map.json"
            if not msa_map_path.is_file():
                logger.error(f"Could not find msa_map.json in {colabfold_output_dir}. Cannot remap template query IDs.")
                return None
            with open(msa_map_path, 'r') as f:
                msa_map = json.load(f)
            
            # The PQT filename is the hash. We create a map from hash -> chain_id
            hash_to_chain_id = {Path(pqt_path).stem.split('.')[0]: header.split('|')[0] for header, pqt_path in msa_map.items()}
            logger.debug(f"Loaded {len(hash_to_chain_id)} chain hashes for template remapping.")

            with open(m8_file, 'r') as f_in, open(hits_m8_path, 'w') as f_out:
                unique_pdb_ids = set()
                num_hits_processed = 0
                for line in f_in:
                    cols = line.strip().split('\t')
                    if len(cols) < 12:
                        continue
                    
                    query_hash, subject_id, pident_str, length_str, mismatch, gapopen, qstart, qend, sstart, send, evalue_str, bitscore = cols[:12]
                    
                    chain_id_for_hit = hash_to_chain_id.get(query_hash)
                    if not chain_id_for_hit:
                        logger.warning(f"Could not find chain ID for template hit with query hash {query_hash}. Skipping.")
                        continue
                        
                    pdb_id, template_chain_id = subject_id.split('_')
                    unique_pdb_ids.add(pdb_id.lower())

                    # --- New Alignment Logic ---
                    # Define file paths
                    full_cif_path = pdb_dir / f"{pdb_id.lower()}.cif"
                    single_chain_cif_path = pdb_dir / f"{pdb_id.lower()}_{template_chain_id}.cif"

                    # Download the full mmCIF if we don't already have it
                    if not full_cif_path.is_file():
                        url = f"https://files.rcsb.org/download/{pdb_id.lower()}.cif"
                        try:
                            response = requests.get(url, timeout=30)
                            response.raise_for_status()
                            import gzip, io
                            raw_bytes = response.content
                            # If the server sent gzipped data but with .cif name, detect and decompress
                            if raw_bytes[:2] == b"\x1f\x8b":
                                try:
                                    raw_bytes = gzip.decompress(raw_bytes)
                                except OSError:
                                    pass  # not valid gzip after all

                            # --- Sanity-check that we actually got an mmCIF ---
                            # A valid mmCIF ALWAYS starts with a data_ block name.
                            trimmed = raw_bytes.lstrip()  # drop leading whitespace/newlines
                            if not trimmed.lower().startswith(b"data_"):
                                snippet = trimmed[:40].replace(b"\n", b" ")
                                logger.warning(
                                    f"Download for {pdb_id.lower()}.cif does not look like CIF; "
                                    f"starts with {snippet!r}. Skipping this template.")
                                continue  # skip this hit entirely

                            clean_bytes = bytes((b if b < 0x80 else 0x3F) for b in raw_bytes)

                            # --- Robust atomic write to bypass potential NFS caching issues ---
                            import tempfile
                            import os
                            temp_fd, temp_path_str = -1, ""
                            try:
                                temp_fd, temp_path_str = tempfile.mkstemp(
                                    dir=full_cif_path.parent, 
                                    prefix=f"{pdb_id.lower()}_", 
                                    suffix=".tmp"
                                )
                                with os.fdopen(temp_fd, 'wb') as f_temp:
                                    temp_fd = -1  # fd is now managed by the file object
                                    f_temp.write(clean_bytes)
                                    f_temp.flush()
                                    os.fsync(f_temp.fileno())  # Force write to disk

                                # Atomically move the file to its final destination
                                os.rename(temp_path_str, full_cif_path)
                                temp_path_str = "" # Prevent cleanup
                            finally:
                                # Cleanup in case of error
                                if temp_fd != -1:
                                    os.close(temp_fd)
                                if temp_path_str and os.path.exists(temp_path_str):
                                    os.unlink(temp_path_str)
                            # --- End robust write ---

                            logger.info(
                                f"Saved mmCIF {full_cif_path.name} (chain extraction pending) – {len(clean_bytes)} bytes"
                            )
                        except requests.exceptions.RequestException as e:
                            logger.warning(f"Failed to download {url}: {e}. Skipping template {subject_id}.")
                            continue

                    # Extract the single chain into its own mmCIF (AF3 requires one chain per file)
                    if not single_chain_cif_path.is_file():
                        try:
                            import gemmi  # Local import to avoid heavy dependency at module import time
                            # --- BEGIN DIAGNOSTIC PRINTS ---
                            print(f"DEBUG: Attempting to use gemmi version {gemmi.__version__}", file=sys.stderr)
                            print(f"DEBUG: Gemmi module loaded from: {gemmi.__file__}", file=sys.stderr)
                            # --- END DIAGNOSTIC PRINTS ---
                            st = gemmi.read_structure(str(full_cif_path))
                            print(f"DEBUG: Type of object returned by read_structure: {type(st)}", file=sys.stderr)
                            model = st[0]

                            chains_to_remove = [ch.name for ch in model if ch.name != template_chain_id]

                            if len(chains_to_remove) == len(model):
                                raise ValueError(f"Chain {template_chain_id} not found in {full_cif_path}")

                            for chain_name in chains_to_remove:
                                model.remove_chain(chain_name)

                            # Write the single-chain mmCIF using the correct method
                            st.write_cif(str(single_chain_cif_path))
                        except Exception as e:
                            logger.warning(f"Failed to extract chain {template_chain_id} from {full_cif_path}: {e}. Skipping template {subject_id}.")
                            continue

                    try:
                        query_sequence = query_sequences[chain_id_for_hit]
                        template_sequence, t_seq_to_idx = template_seq_and_index(str(single_chain_cif_path), template_chain_id)

                        # Align the full sequences to derive an exact mapping
                        q_aln, t_aln = kalign_pair(query_sequence, template_sequence)
                        q2t_map = build_mapping(q_aln, t_aln, t_seq_to_idx)

                        # Quality filters
                        identity = float(pident_str)
                        coverage = len(q2t_map) / len(query_sequence)
                        if coverage < 0.3 or identity < 20.0 or float(evalue_str) > 1e-3:
                            logger.debug(
                                f"Skipping low-quality template {subject_id} for chain {chain_id_for_hit} "
                                f"(Coverage: {coverage:.2f}, Identity: {identity:.2f}, E-value: {evalue_str})"
                            )
                            continue

                    except Exception as e:
                        logger.warning(f"Failed to process alignment for template {subject_id} for chain {chain_id_for_hit}: {e}")
                        continue
                    # --- End Alignment Logic ---

                    # Write the entry for Chai’s hits.m8 (query_id is the real chain ID, not the hash)
                    new_m8_line = (
                        f"{chain_id_for_hit}\t{subject_id}\t{pident_str}\t{length_str}\t{mismatch}\t{gapopen}\t"
                        f"{qstart}\t{qend}\t{sstart}\t{send}\t{evalue_str}\t{bitscore}\n"
                    )
                    f_out.write(new_m8_line)

                    # Store for Boltz & AF3
                    export = TemplateExport(
                        pdb_id=pdb_id.lower(),
                        chain_id=template_chain_id,
                        cif_path=single_chain_cif_path.resolve(),
                        query_idx_to_template_idx=q2t_map,
                        e_value=float(evalue_str),
                        hit_from_chain=chain_id_for_hit
                    )
                    all_exports.append(export)
                    num_hits_processed += 1
            
            logger.info(f"Processed and kept {num_hits_processed} high-quality template hits from ColabFold.")

            # Downloading is now handled within the loop
            # logger.info(f"Fetching {len(unique_pdb_ids)} unique mmCIF files from RCSB.")
            # for pdb_id in unique_pdb_ids:
            #     cif_path = pdb_dir / f"{pdb_id}.cif"
            #     if not cif_path.exists():
            #         url = f"https://files.rcsb.org/download/{pdb_id}.cif"
            #         try:
            #             response = requests.get(url, timeout=30)
            #             response.raise_for_status()
            #             cif_path.write_text(response.text)
            #         except requests.exceptions.RequestException as e:
            #             logger.warning(f"Failed to download {url}: {e}")

            with open(mapping_pkl_path, "wb") as pkl_f:
                pickle.dump(all_exports, pkl_f)

            logger.info(f"Successfully created template store from ColabFold at: {template_store_dir}")
            return str(template_store_dir.resolve())

        except Exception as e:
            logger.error(f"Failed to process ColabFold templates: {e}", exc_info=True)
            return None

    def generate_msa(self) -> Optional[Dict[str, Any]]:
        """
        Generates MSA based on the configuration.

        Relies on self.job_input being set by the caller (Orchestrator).

        Returns:
            A dictionary containing MSA results (e.g., {"af3_data_json": path}
            or {"protein_id_to_a3m_path": {id: path}}), or None if MSA generation failed.
        """
        if self.job_input is None:
            logger.error("Cannot generate MSA: job_input has not been set.")
            return None
            
        if self.job_input.has_msa and \
           (self.job_input.af3_data_json or \
            self.job_input.protein_id_to_a3m_path or \
            self.job_input.input_msa_paths):
            logger.info("MSA data found or indicated in job_input. Skipping MSA generation step.")
            return {}

        msa_method = self.config.get("msa_method_preference", "alphafold3").lower()
        logger.info(f"MSA generation requested using method: {msa_method}")

        if msa_method == "alphafold3":
            return self._run_alphafold3_msa_pipeline()
        elif msa_method == "colabfold":
            return self._run_colabfold_msa_pipeline()
        else:
            logger.error(f"Unsupported MSA method configured: {msa_method}")
            return None

    def _generate_temp_af3_json_for_msa(self) -> Optional[Path]:
        """Generates a temporary AF3 JSON suitable for the data pipeline."""
        if not self.job_input:
            return None
        
        temp_job_input_obj = JobInput(
            name_stem=self.job_input.name_stem + "_msa_gen", 
            sequences=self.job_input.sequences,
            raw_input_type=self.job_input.raw_input_type, 
            output_dir=str(self.msa_tmp_dir), 
            input_msa_paths={}, 
            constraints=None, 
            has_msa=False, 
            af3_data_json=None, 
            protein_id_to_a3m_path={}, 
            model_seeds=None, 
            bonded_atom_pairs=None, 
            is_boltz_config=False,
            is_af3_msa_config_only=True 
        )
        
        try:
            temp_json_filename = f"{temp_job_input_obj.name_stem}_af3_msa_pipeline_input.json" 
            temp_json_path = self.config_generator._generate_af3_json_from_job_input(
                temp_job_input_obj, 
                self.msa_tmp_dir, 
                filename=temp_json_filename 
            )
            if temp_json_path:
                 logger.info(f"Generated temporary AF3 JSON for MSA pipeline: {temp_json_path}")
            return temp_json_path
        except Exception as e:
            logger.error(f"Failed to generate temporary AF3 JSON for MSA pipeline: {e}", exc_info=True)
            return None

    def _run_alphafold3_msa_pipeline(self) -> Optional[Dict[str, Any]]:
        """
        Runs the AlphaFold 3 data pipeline to generate MSAs.
        Uses user's original AF3 JSON if provided and suitable, otherwise generates a temporary one.
        """
        logger.info("Starting AlphaFold 3 MSA pipeline.")
        if not self.job_input: return None

        af3_sif_path = self.config.get("alphafold3_sif_path")
        if not af3_sif_path or not Path(af3_sif_path).is_file():
            logger.error("AlphaFold 3 SIF path not configured or not found.")
            return None

        temp_input_json_path: Optional[Path] = None
        internal_json_name_stem: Optional[str] = None

        if self.job_input.raw_input_type == "af3_json" and \
           self.job_input.original_af3_config_path and \
           not self.job_input.has_msa: 
            
            original_path = Path(self.job_input.original_af3_config_path)
            if original_path.is_file():
                temp_input_json_path = original_path
                internal_json_name_stem = self._get_internal_name_from_af3_json(temp_input_json_path)
                if not internal_json_name_stem:
                    logger.error(f"Could not read internal 'name' from user-provided AF3 JSON: {temp_input_json_path}. Cannot determine AF3 output dir name.")
                    return None
                logger.info(f"Using user-provided AF3 JSON for MSA pipeline: {temp_input_json_path}. Internal name: {internal_json_name_stem}")
            else:
                logger.warning(f"Original AF3 config path {original_path} not found, will generate a temporary one.")
        
        if not temp_input_json_path: 
            temp_input_json_path = self._generate_temp_af3_json_for_msa()
            if not temp_input_json_path: return None
            internal_json_name_stem = self.job_input.name_stem + "_msa_gen"
            logger.info(f"Generated temporary AF3 JSON for MSA pipeline. Internal name: {internal_json_name_stem}")

        if not temp_input_json_path or not internal_json_name_stem:
            logger.error("Failed to determine input JSON or internal name for AF3 MSA pipeline.")
            return None
            
        input_json_filename = temp_input_json_path.name

        af3_data_pipeline_host_output_base = self.msa_tmp_dir

        container_input_json = f"/app/input/{input_json_filename}"
        container_output_dir_in_af3 = "/app/output" 
        container_db_dir = "/data/public_databases" 
        container_model_dir = "/data/models"     

        cmd = ["singularity", "run", "--nv"] 

        # 1. pipeline.py
        custom_pipeline_py_host_path_str = "omnifold/singularity_af3/alphafold3_venv/lib/python3.11/site-packages/alphafold3/data/pipeline.py"
        custom_pipeline_py_host_path = Path(custom_pipeline_py_host_path_str)
        container_pipeline_py_path = "/alphafold3_venv/lib/python3.11/site-packages/alphafold3/data/pipeline.py"
        if custom_pipeline_py_host_path.is_file():
            cmd.extend(["-B", f"{custom_pipeline_py_host_path.resolve()}:{container_pipeline_py_path}:ro"])
            logger.info(f"Binding custom pipeline.py: {custom_pipeline_py_host_path.resolve()} -> {container_pipeline_py_path}")
        else:
            logger.warning(f"Custom pipeline.py not found at {custom_pipeline_py_host_path.resolve()}. Using SIF default.")

        # 2. run_alphafold.py
        custom_run_alphafold_py_host_path_str = "omnifold/singularity_af3/alphafold3/run_alphafold.py"
        custom_run_alphafold_py_host_path = Path(custom_run_alphafold_py_host_path_str)
        container_run_alphafold_py_path = "/alphafold3/run_alphafold.py"
        if custom_run_alphafold_py_host_path.is_file():
            cmd.extend(["-B", f"{custom_run_alphafold_py_host_path.resolve()}:{container_run_alphafold_py_path}:ro"])
            logger.info(f"Binding custom run_alphafold.py: {custom_run_alphafold_py_host_path.resolve()} -> {container_run_alphafold_py_path}")
        else:
            logger.warning(f"Custom run_alphafold.py not found at {custom_run_alphafold_py_host_path.resolve()}. Using SIF default.")

        # 3. templates.py
        custom_templates_py_host_path_str = "omnifold/singularity_af3/alphafold3_venv/lib/python3.11/site-packages/alphafold3/data/templates.py"
        custom_templates_py_host_path = Path(custom_templates_py_host_path_str)
        container_templates_py_path = "/alphafold3_venv/lib/python3.11/site-packages/alphafold3/data/templates.py"
        if custom_templates_py_host_path.is_file():
            cmd.extend(["-B", f"{custom_templates_py_host_path.resolve()}:{container_templates_py_path}:ro"])
            logger.info(f"Binding custom templates.py: {custom_templates_py_host_path.resolve()} -> {container_templates_py_path}")
        else:
            logger.warning(f"Custom templates.py not found at {custom_templates_py_host_path.resolve()}. Using SIF default.")

        # 4. tools/hmmsearch.py
        custom_hmmsearch_py_host_path_str = "omnifold/singularity_af3/alphafold3_venv/lib/python3.11/site-packages/alphafold3/data/tools/hmmsearch.py"
        custom_hmmsearch_py_host_path = Path(custom_hmmsearch_py_host_path_str)
        container_hmmsearch_py_path = "/alphafold3_venv/lib/python3.11/site-packages/alphafold3/data/tools/hmmsearch.py"
        if custom_hmmsearch_py_host_path.is_file():
            cmd.extend(["-B", f"{custom_hmmsearch_py_host_path.resolve()}:{container_hmmsearch_py_path}:ro"])
            logger.info(f"Binding custom hmmsearch.py: {custom_hmmsearch_py_host_path.resolve()} -> {container_hmmsearch_py_path}")
        else:
            logger.warning(f"Custom hmmsearch.py not found at {custom_hmmsearch_py_host_path.resolve()}. Using SIF default.")

        cmd.extend(["-B", f"{temp_input_json_path.parent.resolve()}:/app/input:ro"])
        cmd.extend(["-B", f"{af3_data_pipeline_host_output_base.resolve()}:{container_output_dir_in_af3}"])
        db_dir_host = self.config.get("alphafold3_database_dir")
        if db_dir_host and Path(db_dir_host).is_dir():
            cmd.extend(["-B", f"{Path(db_dir_host).resolve()}:{container_db_dir}:ro"])
        else:
            logger.error(f"AlphaFold DB directory not configured or found: {db_dir_host}. Cannot run AF3 MSA pipeline.")
            return None
        model_dir_host = self.config.get("alphafold3_model_weights_dir")
        if model_dir_host and Path(model_dir_host).is_dir():
             cmd.extend(["-B", f"{Path(model_dir_host).resolve()}:{container_model_dir}:ro"])
        else:
             logger.error(f"AlphaFold model weights directory not configured or found: {model_dir_host}. Cannot run AF3 MSA pipeline.")
             return None

        cmd.append(af3_sif_path)

        run_script_args = [
            f"--json_path={container_input_json}",
            f"--output_dir={container_output_dir_in_af3}", 
            f"--model_dir={container_model_dir}",
            "--run_data_pipeline=True",
            "--run_inference=False", 
        ]
        if db_dir_host and Path(db_dir_host).is_dir():
            run_script_args.append(f"--db_dir={container_db_dir}")
        
        cmd.extend(run_script_args)

        exit_code, stdout, stderr = self._run_command(cmd)

        if exit_code != 0:
            logger.error(f"AlphaFold 3 MSA pipeline failed with exit code {exit_code}.")
            return None
        

        # Determine the name_prefix from the input JSON file, which dictates the subdirectory AF3 creates
        input_json_path_obj = Path(temp_input_json_path)

        try:
            with open(input_json_path_obj, 'r') as f_json:
                input_data = json.load(f_json)
            name_from_json_field = input_data.get("name")

            if name_from_json_field and isinstance(name_from_json_field, str) and name_from_json_field.strip():
                # Use the "name" field from the JSON file, sanitize, and lowercase
                sanitized_name = "".join(c if c.isalnum() or c in ['_','-'] else '_' for c in name_from_json_field)
                name_prefix = sanitized_name.lower()
                logger.info(f"Derived name_prefix '{name_prefix}' from JSON 'name' field ('{name_from_json_field}'), sanitized and lowercased.")
            else:
                # Fallback to the input JSON filename's stem, and lowercase
                name_prefix = input_json_path_obj.stem.lower()
                logger.info(f"Derived name_prefix '{name_prefix}' from input JSON filename stem ('{input_json_path_obj.stem}'), lowercased (no valid 'name' field found in JSON).")
        except Exception as e_json:
            # Fallback if JSON parsing fails or 'name' field access has issues
            logger.warning(f"Could not robustly read 'name' field from {input_json_path_obj} (Error: {e_json}). "
                           f"Falling back to filename stem ('{input_json_path_obj.stem}'), lowercased.")
            name_prefix = input_json_path_obj.stem.lower()

        logger.info(f"Final name_prefix for AF3 output directory determination: '{name_prefix}'")
        
        # This is the base directory on the host where AF3 outputs will be written to by Singularity
        # (e.g., /path/to/msa_intermediate_files)
        af3_msa_output_base_dir_host = self.msa_tmp_dir

        # Path to the AF3-generated _data.json file, considering AF3 creates a subdirectory named 'name_prefix'
        output_data_json_path = af3_msa_output_base_dir_host / name_prefix / f"{name_prefix}_data.json"

        if output_data_json_path.is_file():
            logger.info(f"AlphaFold 3 MSA data JSON found at: {output_data_json_path.resolve()}")
            results = {"af3_data_json": str(output_data_json_path.resolve())}

            af3_msa_output_actual_dir = output_data_json_path.parent # e.g. .../msa_intermediate_files/job_name_msa_gen/

            # --- Process Templates for Chai/Boltz ---
            template_store_path = self._process_template_metadata(af3_msa_output_actual_dir)
            if template_store_path:
                results["template_store_path"] = template_store_path
            else:
                logger.warning("Template processing failed or produced no output. Chai/Boltz may not use templates.")
            # --- End Template Processing ---


            # --- Extract Unpaired MSAs from JSON for af3_to_boltz_csv.py and general use ---
            json_extracted_unpaired_msas_output_dir = af3_msa_output_actual_dir / "json_extracted_unpaired_msas"
            json_extracted_unpaired_msas_output_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"Extracting unpaired MSAs from {output_data_json_path} to {json_extracted_unpaired_msas_output_dir}")
            
            protein_ids_in_job = []
            if self.job_input and self.job_input.sequences:
                protein_ids_in_job = [seq.chain_id for seq in self.job_input.sequences if seq.molecule_type == 'protein']

            # Call the extraction utility
            # Note: extract_all_protein_a3ms_from_af3_json returns a dict {protein_id: path}
            # This utility handles cases where a protein entity in JSON might have multiple IDs (homomers)
            extracted_json_unpaired_a3m_paths: Optional[Dict[str, str]] = extract_all_protein_a3ms_from_af3_json(
                json_path=str(output_data_json_path.resolve()),
                output_dir=str(json_extracted_unpaired_msas_output_dir.resolve())
            )

            if extracted_json_unpaired_a3m_paths is not None:
                logger.info(f"Successfully extracted {len(extracted_json_unpaired_a3m_paths)} unpaired A3M files from JSON.")
                results["protein_id_to_json_unpaired_a3m_path"] = extracted_json_unpaired_a3m_paths
            else:
                logger.error(f"Failed to extract unpaired A3Ms from JSON {output_data_json_path}. This might affect Boltz CSV generation.")
            # ---

            chai1_sif_path_str = self.config.get("chai1_sif_path")
            if chai1_sif_path_str and Path(chai1_sif_path_str).is_file():
                logger.info("Chai-1 SIF found. Attempting A3M to PQT conversion.")
                source_msas_dir = af3_msa_output_actual_dir / "msas"
                target_pqt_dir = af3_msa_output_actual_dir / "msas_forChai"
                
                if not source_msas_dir.is_dir():
                    logger.warning(f"AF3 MSA 'msas' directory not found at {source_msas_dir}. Skipping PQT conversion.")
                else:
                    target_pqt_dir.mkdir(parents=True, exist_ok=True)
                    logger.info(f"Created PQT output directory: {target_pqt_dir}")
                    
                    processed_any_pqt = False
                    for chain_msa_dir in source_msas_dir.iterdir():
                        if chain_msa_dir.is_dir():
                            chain_name = chain_msa_dir.name # e.g., "chain_A"
                            # Define container paths for binding
                            container_input_a3ms = "/input_a3ms"
                            container_output_pqts_dir = "/output_pqts"
                            # Output filename is determined by chai tool using a hash of the sequence.
                            # So we only provide the output directory to the tool.

                            pqt_cmd = [
                                "singularity", "exec",
                                "-B", f"{chain_msa_dir.resolve()}:{container_input_a3ms}:ro",
                                "-B", f"{target_pqt_dir.resolve()}:{container_output_pqts_dir}",
                                chai1_sif_path_str,
                                "chai-lab", "a3m-to-pqt",
                                container_input_a3ms,
                                "--output-directory", container_output_pqts_dir
                            ]
                            
                            logger.info(f"Running PQT conversion for {chain_name}: {' '.join(pqt_cmd)}")
                            # Store current .pqt files to check for new ones later
                            pqt_files_before = set(f.name for f in target_pqt_dir.glob("*.pqt"))
                            exit_code, stdout, stderr = self._run_command(pqt_cmd)

                            if exit_code == 0:
                                pqt_files_after = set(f.name for f in target_pqt_dir.glob("*.pqt"))
                                new_pqt_files = pqt_files_after - pqt_files_before
                                if new_pqt_files:
                                    logger.info(f"Successfully converted A3M to PQT for {chain_name}. New PQT file(s): {', '.join(new_pqt_files)} in {target_pqt_dir}")
                                    processed_any_pqt = True
                                else:
                                    # This case might occur if the tool ran successfully (exit 0) but didn't produce output for some reason (e.g. empty input)
                                    logger.warning(f"PQT conversion for {chain_name} exited successfully but no new .pqt file was found in {target_pqt_dir}.")
                                    logger.debug(f"PQT conversion STDOUT for {chain_name}:\n{stdout}")
                                    logger.debug(f"PQT conversion STDERR for {chain_name}:\n{stderr}")
                            else:
                                logger.error(f"PQT conversion failed for {chain_name}. Exit code: {exit_code}")
                                logger.error(f"PQT conversion STDOUT for {chain_name}:\n{stdout}")
                    
                    # After attempting all chains, check if any .pqt files exist in the target directory.
                    pqt_files_in_target = list(target_pqt_dir.glob("*.pqt"))
                    if pqt_files_in_target:
                        logger.info(f"PQT conversion successful. Found {len(pqt_files_in_target)} .pqt file(s) in {target_pqt_dir}: {[f.name for f in pqt_files_in_target]}")
                        results["chai_pqt_msa_dir"] = str(target_pqt_dir.resolve())
                        
                        # --- Convert AF3 JSON to Chai FASTA ---
                        chai_fasta_path = af3_msa_output_actual_dir / "chai_input.fasta"
                        logger.info(f"Attempting to convert AF3 JSON {output_data_json_path} to Chai FASTA {chai_fasta_path}")
                        conversion_successful = af3_json_to_chai_fasta(
                            json_path=output_data_json_path, 
                            fasta_path=chai_fasta_path
                        )
                        if conversion_successful:
                            results["chai_fasta_path"] = str(chai_fasta_path.resolve())
                            logger.info(f"Successfully generated Chai FASTA: {chai_fasta_path}")
                        else:
                            logger.error(f"Failed to generate Chai FASTA from AF3 JSON {output_data_json_path}.")
                        # --- End of AF3 JSON to Chai FASTA Conversion ---
                    else:
                        logger.warning(f"No .pqt files found in {target_pqt_dir} after attempting conversion for all chains. Skipping Chai FASTA generation.")
            else:
                logger.info("Chai-1 SIF not provided or not found. Skipping PQT conversion and Chai FASTA generation.")
            
            # --- Run af3_to_boltz_csv.py script ---
            if self.job_input and self.job_input.sequences:
                # af3_msa_output_actual_dir is already defined above
                source_af3_paired_msas_root_dir = af3_msa_output_actual_dir / "msas" # This is the <msa_root> for paired UniProt A3Ms
                boltz_csv_output_dir = af3_msa_output_actual_dir / "boltz_csv_msas"
                boltz_csv_output_dir.mkdir(parents=True, exist_ok=True)

                # chain_ids_for_script was already defined as protein_ids_in_job
                
                if protein_ids_in_job and source_af3_paired_msas_root_dir.is_dir() and json_extracted_unpaired_msas_output_dir.is_dir():
                    script_path = Path(__file__).parent / "util" / "af3_to_boltz_csv.py"
                    if not script_path.is_file():
                        logger.error(f"af3_to_boltz_csv.py script not found at {script_path}. Cannot generate Boltz CSVs.")
                    else:
                        cmd_boltz_csv = [
                            sys.executable, str(script_path),
                            "--chains", ",".join(protein_ids_in_job),
                            "--msa_root", str(source_af3_paired_msas_root_dir.resolve()), # For uniprot_*.a3m files
                            "--json_extracted_unpaired_msa_dir", str(json_extracted_unpaired_msas_output_dir.resolve()), # For msa_*.a3m files
                            "--out", str(boltz_csv_output_dir.resolve()),
                            "--max_paired", str(self.config.get("boltz_csv_max_paired_msas", 256)),
                            "--max_total", str(self.config.get("boltz_csv_max_total_msas", 4096))
                        ]
                        if self.config.get("boltz_csv_shuffle_paired_msas", False):
                            cmd_boltz_csv.append("--shuffle_paired")
                        
                        logger.info(f"Running af3_to_boltz_csv.py: {' '.join(cmd_boltz_csv)}")
                        exit_code_csv, stdout_csv, stderr_csv = self._run_command(cmd_boltz_csv)

                        if exit_code_csv == 0:
                            # Check if any CSV files were actually created
                            if any(boltz_csv_output_dir.glob("*.csv")):
                                logger.info(f"Successfully generated Boltz CSV MSAs in {boltz_csv_output_dir}")
                                results["boltz_csv_msa_dir"] = str(boltz_csv_output_dir.resolve())
                            else:
                                logger.warning(f"af3_to_boltz_csv.py ran successfully but no CSV files found in {boltz_csv_output_dir}. Stdout: {stdout_csv}. Stderr: {stderr_csv}")
                        else:
                            logger.error(f"af3_to_boltz_csv.py failed with exit code {exit_code_csv}. Stdout: {stdout_csv}. Stderr: {stderr_csv}")
                elif not protein_ids_in_job:
                    logger.info("No protein chains found in job input. Skipping Boltz CSV generation.")
                elif not source_af3_paired_msas_root_dir.is_dir():
                    logger.warning(f"AF3 'msas' directory for Boltz CSV (paired) conversion not found at {source_af3_paired_msas_root_dir}. Check AF3 output. Skipping CSV gen.")
                elif not json_extracted_unpaired_msas_output_dir.is_dir():
                     logger.warning(f"Directory for JSON-extracted unpaired MSAs ({json_extracted_unpaired_msas_output_dir}) not found or extraction failed. Skipping Boltz CSV gen.")
            # --- End of af3_to_boltz_csv.py script run ---

            return results
        else:
            logger.error(f"AlphaFold 3 MSA output data JSON not found at expected location: {output_data_json_path}")
            logger.info(f"Check AF3 pipeline stdout/stderr in debug logs if needed.")
            return None

    def _generate_temp_fasta(self) -> Optional[Path]:
         """Creates a temporary FASTA file from job_input sequences."""
         if not self.job_input or not self.job_input.sequences:
             logger.error("No sequences available to create FASTA.")
             return None
         
         fasta_path = self.msa_tmp_dir / f"{self.job_input.name_stem}_colabfold_input.fasta"
         try:
             with open(fasta_path, 'w') as f:
                 for seq_info in self.job_input.sequences:
                     f.write(f">{seq_info.chain_id}|{seq_info.original_name}\n") 
                     f.write(f"{seq_info.sequence}\n")
             logger.info(f"Generated temporary FASTA for ColabFold MSA: {fasta_path}")
             return fasta_path
         except IOError as e:
             logger.error(f"Failed to write temporary FASTA file {fasta_path}: {e}")
             return None


    def _run_colabfold_msa_pipeline(self) -> Optional[Dict[str, Any]]:
        """
        Runs the ColabFold MSA generation script.
        This generates both A3M and Chai-compatible PQT files.
        """
        logger.info("Starting ColabFold MSA pipeline.")
        if not self.job_input:
            return None

        temp_fasta_path = self._generate_temp_fasta()
        if not temp_fasta_path:
            return None

        colabfold_output_dir = self.msa_tmp_dir / f"{self.job_input.name_stem}_colabfold"
        colabfold_output_dir.mkdir(exist_ok=True)

        # Use sys.executable to ensure we use the same python interpreter
        # that is running the main application.
        cmd = [
            sys.executable,
            "-m", "omnifold.util.generate_colabfold_msas",
            str(temp_fasta_path),
            "--out_dir", str(colabfold_output_dir),
            "--write_a3m", # Always write A3M for Boltz/AF3 consumption
            "--include_templates" # Always search for templates
        ]
        
        # Note: We are not exposing the --include_templates flag from the underlying
        # script to the main pipeline config yet. This could be added later.

        exit_code, stdout, stderr = self._run_command(cmd)

        if exit_code != 0:
            logger.error(f"ColabFold MSA script failed with exit code {exit_code}.")
            logger.error(f"Stderr: {stderr}")
            return None

        # --- Parse the output ---
        # The script is expected to produce a manifest file mapping chain IDs
        # to their .pqt files. It also produces .a3m files that we need to find.
        manifest_path = colabfold_output_dir / "msa_map.json"
        if not manifest_path.is_file():
            logger.error(f"ColabFold MSA script did not produce the expected manifest file: {manifest_path}")
            return None

        with open(manifest_path, 'r') as f:
            pqt_manifest = json.load(f)

        protein_id_to_pqt_path = {
            header.split('|')[0]: path for header, path in pqt_manifest.items()
        }
        
        # The A3M files are in a subdirectory and are named by hash.
        # We need to map them back to the protein IDs. We can do this by
        # leveraging the PQT manifest, as the base names match.
        a3m_dir = colabfold_output_dir / "a3ms"
        protein_id_to_unpaired_a3m_path = {}
        protein_id_to_paired_a3m_path = {}

        for header, pqt_path_str in pqt_manifest.items():
            chain_id = header.split('|')[0]
            # The PQT filename is <HASH>.aligned.pqt. We need to extract the hash.
            pqt_filename = Path(pqt_path_str).name
            file_hash = pqt_filename.split(".")[0]
            
            # The A3M files are named <HASH>.pair.a3m and <HASH>.single.a3m
            pair_a3m = a3m_dir / f"{file_hash}.pair.a3m"
            single_a3m = a3m_dir / f"{file_hash}.single.a3m"
            
            if pair_a3m.is_file():
                protein_id_to_paired_a3m_path[chain_id] = str(pair_a3m)
            if single_a3m.is_file():
                protein_id_to_unpaired_a3m_path[chain_id] = str(single_a3m)

        if not protein_id_to_unpaired_a3m_path or not protein_id_to_pqt_path:
            logger.error("Failed to map any A3M or PQT files from ColabFold output.")
            return None
            
        logger.info(f"Successfully generated MSAs using ColabFold for {len(protein_id_to_unpaired_a3m_path)} chains.")
        
        # Create the correctly formatted FASTA for Chai-1 using the new utility
        chai_fasta_path = colabfold_output_dir / f"{self.job_input.name_stem}_chai_input.fasta"
        if not job_input_to_chai_fasta(self.job_input, chai_fasta_path):
            logger.error("Failed to generate the correctly formatted FASTA for Chai-1.")
            # We can still proceed without it if Chai-1 is not being run.
            chai_fasta_path = None
        
        # --- Process Templates ---
        template_store_path = self._process_colabfold_templates(colabfold_output_dir)
        
        final_results = {
            "source": "colabfold",
            "protein_id_to_a3m_path": {
                "unpaired": protein_id_to_unpaired_a3m_path,
                "paired": protein_id_to_paired_a3m_path
            },
            "protein_id_to_pqt_path": protein_id_to_pqt_path,
            "chai_fasta_path": str(chai_fasta_path) if chai_fasta_path else None,
            "template_store_path": template_store_path
        }
        logger.info(f"--- MSAManager ColabFold Results ---\n{json.dumps(final_results, indent=2)}\n------------------------------------")
        return final_results

