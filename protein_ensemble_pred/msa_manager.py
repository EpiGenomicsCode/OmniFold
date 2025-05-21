import json
import os
import subprocess
import logging
import tempfile # For temporary FASTA
from pathlib import Path
from typing import Any, Dict, Optional, Tuple, List

from .config_generator import ConfigGenerator
from .util.definitions import JobInput, SequenceInfo

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
            logger.debug(f"Command stdout: {process.stdout}")
            if process.returncode != 0:
                 logger.warning(f"Command stderr: {process.stderr}")
            return process.returncode, process.stdout, process.stderr
        except FileNotFoundError:
            logger.error(f"Command not found: {cmd[0]}. Ensure it's installed and in PATH.")
            return -1, "", f"Command not found: {cmd[0]}"
        except Exception as e:
             logger.error(f"Error running command {' '.join(cmd)}: {e}", exc_info=True)
             return -1, "", str(e)

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
        custom_pipeline_py_host_path_str = "protein_ensemble_pred/singularity_af3/alphafold3_venv/lib/python3.11/site-packages/alphafold3/data/pipeline.py"
        custom_pipeline_py_host_path = Path(custom_pipeline_py_host_path_str)
        container_pipeline_py_path = "/alphafold3_venv/lib/python3.11/site-packages/alphafold3/data/pipeline.py"
        if custom_pipeline_py_host_path.is_file():
            cmd.extend(["-B", f"{custom_pipeline_py_host_path.resolve()}:{container_pipeline_py_path}:ro"])
            logger.info(f"Binding custom pipeline.py: {custom_pipeline_py_host_path.resolve()} -> {container_pipeline_py_path}")
        else:
            logger.warning(f"Custom pipeline.py not found at {custom_pipeline_py_host_path.resolve()}. Using SIF default.")

        # 2. run_alphafold.py
        custom_run_alphafold_py_host_path_str = "protein_ensemble_pred/singularity_af3/alphafold3/run_alphafold.py"
        custom_run_alphafold_py_host_path = Path(custom_run_alphafold_py_host_path_str)
        container_run_alphafold_py_path = "/alphafold3/run_alphafold.py"
        if custom_run_alphafold_py_host_path.is_file():
            cmd.extend(["-B", f"{custom_run_alphafold_py_host_path.resolve()}:{container_run_alphafold_py_path}:ro"])
            logger.info(f"Binding custom run_alphafold.py: {custom_run_alphafold_py_host_path.resolve()} -> {container_run_alphafold_py_path}")
        else:
            logger.warning(f"Custom run_alphafold.py not found at {custom_run_alphafold_py_host_path.resolve()}. Using SIF default.")

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
        

        expected_output_subdirectory = af3_data_pipeline_host_output_base / internal_json_name_stem
        expected_output_filename = f"{internal_json_name_stem}_data.json"
        output_data_json_path = expected_output_subdirectory / expected_output_filename
        
        if output_data_json_path.is_file():
            logger.info(f"AlphaFold 3 MSA data JSON found at: {output_data_json_path.resolve()}")
            results = {"af3_data_json": str(output_data_json_path.resolve())}

            chai1_sif_path_str = self.config.get("chai1_sif_path")
            if chai1_sif_path_str and Path(chai1_sif_path_str).is_file():
                logger.info("Chai-1 SIF found. Attempting A3M to PQT conversion.")
                af3_msa_output_base_dir = output_data_json_path.parent
                source_msas_dir = af3_msa_output_base_dir / "msas"
                target_pqt_dir = af3_msa_output_base_dir / "msas_forChai"
                
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
                                "chai", "a3m-to-pqt",
                                "-i", container_input_a3ms, # Input is the DIRECTORY
                                "-o", container_output_pqts_dir  # Output is also the DIRECTORY
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
                    
                    if processed_any_pqt:
                        results["chai_pqt_msa_dir"] = str(target_pqt_dir.resolve())
            else:
                logger.info("Chai-1 SIF not provided or not found. Skipping A3M to PQT conversion.")
            
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
        Runs a ColabFold/MMseqs2 MSA pipeline using the Boltz container.
        Assumes Boltz container has a CLI mode for local MMseqs2 search.
        """
        logger.info("Starting ColabFold (MMseqs2) MSA pipeline via Boltz container.")
        if not self.job_input: return None

        boltz_sif_path = self.config.get("boltz_sif_path")
        if not boltz_sif_path or not Path(boltz_sif_path).is_file():
            logger.error("Boltz SIF path not configured or not found.")
            return None
            
        temp_fasta_path = self._generate_temp_fasta()
        if not temp_fasta_path:
             return None
        input_fasta_filename = temp_fasta_path.name
        
        colabfold_msa_output_dir = self.msa_tmp_dir / f"{self.job_input.name_stem}_colabfold_a3m_out"
        colabfold_msa_output_dir.mkdir(parents=True, exist_ok=True)

        container_input_fasta = f"/app/input/{input_fasta_filename}"
        container_output_dir = "/app/output"

        cmd = ["singularity", "exec"]
        
        cmd.extend(["-B", f"{temp_fasta_path.parent.resolve()}:/app/input:ro"])
        cmd.extend(["-B", f"{colabfold_msa_output_dir.resolve()}:{container_output_dir}"])
        
        cmd.append(boltz_sif_path)
        
        cmd.extend([
            "boltz",
            "predict",
            container_input_fasta, 
            "--output-dir", container_output_dir,
            "--use_msa_server", 
            "--msa_only",      
        ])
        
        msa_server_url = self.config.get("colabfold_msa_server_url")
        if msa_server_url:
            cmd.extend(["--msa_server_url", msa_server_url])
            logger.info(f"Using provided MSA server URL: {msa_server_url}")
        else:
             logger.info("Using Boltz internal default MSA server URL.")

        exit_code, stdout, stderr = self._run_command(cmd)

        if exit_code != 0:
            logger.error(f"ColabFold/Boltz MSA pipeline failed with exit code {exit_code}.")
            return None
            
        protein_id_to_a3m_path: Dict[str, str] = {}
        found_a3m = False
        try:
            for filename in os.listdir(colabfold_msa_output_dir):
                 if filename.endswith(".a3m"):
                     parts = filename.split('.a3m')[0]
                     chain_id = parts 
                     if any(s.chain_id == chain_id for s in self.job_input.sequences):
                         abs_path = str((colabfold_msa_output_dir / filename).resolve())
                         protein_id_to_a3m_path[chain_id] = abs_path
                         found_a3m = True
                         logger.info(f"Found generated A3M for chain {chain_id}: {abs_path}")
                     else:
                          logger.warning(f"Found A3M file {filename} but corresponding chain ID {chain_id} not in original input. Skipping.")
            
            if not found_a3m:
                 logger.error(f"ColabFold/Boltz MSA pipeline ran but no A3M files found in {colabfold_msa_output_dir}")
                 return None
                 
            return {"protein_id_to_a3m_path": protein_id_to_a3m_path}
        
        except Exception as e:
             logger.error(f"Error processing ColabFold/Boltz MSA output directory {colabfold_msa_output_dir}: {e}", exc_info=True)
             return None

