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
        self.output_dir = Path(output_dir) # Ensure it's a Path object
        self.msa_tmp_dir = self.output_dir / "msa_intermediate_files"
        self.msa_tmp_dir.mkdir(parents=True, exist_ok=True)
        # Instantiate ConfigGenerator internally if needed for temp AF3 JSON
        self.config_generator = ConfigGenerator()
        self.job_input: Optional[JobInput] = None # Set by Orchestrator before calling generate_msa
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
                check=False, # Don't raise exception on non-zero exit
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
            
        # Check if MSA was already provided in a way that bypasses generation
        # (e.g., input was data.json or Boltz YAML with MSA paths)
        if self.job_input.has_msa and \
           (self.job_input.af3_data_json or \
            self.job_input.protein_id_to_a3m_path or \
            self.job_input.input_msa_paths):
            logger.info("MSA data found or indicated in job_input. Skipping MSA generation step.")
            # Return empty dict as signal that generation was skipped but okay
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
        
        # Create a new JobInput that *doesn't* contain any MSA paths or existing MSA data for the temp config
        temp_job_input_obj = JobInput(
            name_stem=self.job_input.name_stem + "_msa_gen", # Ensure a unique name for temp files
            sequences=self.job_input.sequences,
            raw_input_type=self.job_input.raw_input_type, # Keep original type for context, though we force AF3 JSON out
            output_dir=str(self.msa_tmp_dir), # Provide the intermediate dir as output_dir context
            input_msa_paths={}, # Explicitly empty
            constraints=None, # Not relevant for AF3 MSA pipeline
            has_msa=False, # Explicitly False for MSA generation run
            af3_data_json=None, # Explicitly None
            protein_id_to_a3m_path={}, # Explicitly empty
            model_seeds=None, # Not relevant for MSA pipeline only
            bonded_atom_pairs=None, # Not relevant for MSA pipeline only
            is_boltz_config=False,
            is_af3_msa_config_only=True # Mark this as an MSA-only config
        )
        
        # Generate the JSON in the intermediate directory
        try:
            # Use a distinct name for the temp file
            temp_json_filename = f"{temp_job_input_obj.name_stem}_af3_msa_pipeline_input.json" # More descriptive name
            temp_json_path = self.config_generator._generate_af3_json_from_job_input(
                temp_job_input_obj, 
                self.msa_tmp_dir, # Save in intermediate dir
                filename=temp_json_filename # Pass the explicit filename
            )
            if temp_json_path:
                 logger.info(f"Generated temporary AF3 JSON for MSA pipeline: {temp_json_path}")
                 # Manually check and remove unpairedMsa/pairedMsa fields just in case?
                 # This is belt-and-suspenders if _generate_af3_json respects empty input_msa_paths
                 # with open(temp_json_path, 'r+') as f:
                 #     data = json.load(f)
                 #     modified = False
                 #     for entity in data.get("sequences", []):
                 #         for key in ["protein", "rna", "dna"]:
                 #              if key in entity and isinstance(entity[key], dict):
                 #                  if "unpairedMsaPath" in entity[key]: del entity[key]["unpairedMsaPath"]; modified=True
                 #                  if "pairedMsaPath" in entity[key]: del entity[key]["pairedMsaPath"]; modified=True
                 #                  if "pairedMsa" in entity[key]: del entity[key]["pairedMsa"]; modified=True 
                 #     if modified:
                 #          f.seek(0)
                 #          json.dump(data, f, indent=2)
                 #          f.truncate()
                 #          logger.info(f"Ensured MSA path fields removed from temp JSON: {temp_json_path}")
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
           not self.job_input.has_msa: # original_af3_config_path exists and it did NOT have MSAs
            
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
        
        if not temp_input_json_path: # Fallback or if input was not AF3 JSON
            temp_input_json_path = self._generate_temp_af3_json_for_msa()
            if not temp_input_json_path: return None
            # The name used inside this temp JSON is job_input.name_stem + "_msa_gen"
            internal_json_name_stem = self.job_input.name_stem + "_msa_gen"
            logger.info(f"Generated temporary AF3 JSON for MSA pipeline. Internal name: {internal_json_name_stem}")

        if not temp_input_json_path or not internal_json_name_stem:
            logger.error("Failed to determine input JSON or internal name for AF3 MSA pipeline.")
            return None
            
        input_json_filename = temp_input_json_path.name

        # 2. Define output directory for this specific AF3 data pipeline run
        # The host output base is self.msa_tmp_dir
        # AF3 will create a subdir named 'internal_json_name_stem' inside the container's output dir.
        # So, the *container's* output dir mapping should point to self.msa_tmp_dir.
        # And the final data will be in self.msa_tmp_dir / internal_json_name_stem / ...
        af3_data_pipeline_host_output_base = self.msa_tmp_dir

        # 3. Construct Singularity command
        container_input_json = f"/app/input/{input_json_filename}"
        container_output_dir_in_af3 = "/app/output" # AF3 will write to /app/output/internal_json_name_stem
        container_db_dir = "/data/public_databases" # Match Runner's AF3 convention
        container_model_dir = "/data/models"      # Match Runner's AF3 convention

        cmd = ["singularity", "run", "--nv"] # Using run as per Runner setup for AF3

        # Bindings (absolute paths on host)
        cmd.extend(["-B", f"{temp_input_json_path.parent.resolve()}:/app/input:ro"])
        cmd.extend(["-B", f"{af3_data_pipeline_host_output_base.resolve()}:{container_output_dir_in_af3}"])
        # Bind DBs - Use the path from the main config
        db_dir_host = self.config.get("alphafold3_database_dir")
        if db_dir_host and Path(db_dir_host).is_dir():
            cmd.extend(["-B", f"{Path(db_dir_host).resolve()}:{container_db_dir}:ro"])
        else:
            logger.error(f"AlphaFold DB directory not configured or found: {db_dir_host}. Cannot run AF3 MSA pipeline.")
            return None
        # Bind Models - Use the path from the main config
        model_dir_host = self.config.get("alphafold3_model_weights_dir")
        if model_dir_host and Path(model_dir_host).is_dir():
             cmd.extend(["-B", f"{Path(model_dir_host).resolve()}:{container_model_dir}:ro"])
        else:
             logger.error(f"AlphaFold model weights directory not configured or found: {model_dir_host}. Cannot run AF3 MSA pipeline.")
             return None

        cmd.append(af3_sif_path)

        # Command arguments for run_alphafold.py inside container
        run_script_args = [
            f"--json_path={container_input_json}",
            f"--output_dir={container_output_dir_in_af3}", # Tell AF3 to write to /app/output
            f"--model_dir={container_model_dir}",
            "--run_data_pipeline=True",
            "--run_inference=False", 
        ]
        if db_dir_host and Path(db_dir_host).is_dir(): # db_dir_host is from self.config.get("alphafold3_database_dir")
            run_script_args.append(f"--db_dir={container_db_dir}")
        
        cmd.extend(run_script_args)

        # 4. Execute command
        exit_code, stdout, stderr = self._run_command(cmd)

        # 5. Check results
        if exit_code != 0:
            logger.error(f"AlphaFold 3 MSA pipeline failed with exit code {exit_code}.")
            return None
        
        # The stem used by AF3 for its output subdirectory and _data.json file is internal_json_name_stem.
        # This was derived either from the user's JSON or generated for the temp JSON.

        # af3_data_pipeline_host_output_base is the HOST path to what was container_output_dir_in_af3 in container
        # AF3 creates a subdirectory inside this based on internal_json_name_stem
        expected_output_subdirectory = af3_data_pipeline_host_output_base / internal_json_name_stem
        expected_output_filename = f"{internal_json_name_stem}_data.json"
        output_data_json_path = expected_output_subdirectory / expected_output_filename
        
        if output_data_json_path.is_file():
            logger.info(f"AlphaFold 3 MSA data JSON found at: {output_data_json_path.resolve()}")
            # Return the path to this comprehensive JSON
            return {"af3_data_json": str(output_data_json_path.resolve())}
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
                     # Include chain ID in header for potential parsing later
                     f.write(f">{seq_info.chain_id}|{seq_info.original_name}\n") 
                     # Write sequence in lines of 60 chars? MMseqs might not care.
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
            
        # 1. Generate temporary input FASTA
        temp_fasta_path = self._generate_temp_fasta()
        if not temp_fasta_path:
             return None
        input_fasta_filename = temp_fasta_path.name
        
        # 2. Define output directory for this specific ColabFold run
        colabfold_msa_output_dir = self.msa_tmp_dir / f"{self.job_input.name_stem}_colabfold_a3m_out"
        colabfold_msa_output_dir.mkdir(parents=True, exist_ok=True)

        # 3. Construct Singularity command (using exec for Boltz)
        container_input_fasta = f"/app/input/{input_fasta_filename}"
        container_output_dir = "/app/output"

        cmd = ["singularity", "exec"]
        
        # Bindings
        cmd.extend(["-B", f"{temp_fasta_path.parent.resolve()}:/app/input:ro"])
        cmd.extend(["-B", f"{colabfold_msa_output_dir.resolve()}:{container_output_dir}"])
        
        cmd.append(boltz_sif_path)
        
        # Command inside Boltz container
        # Calls 'boltz predict' but uses flags to only run the MSA step via server
        cmd.extend([
            "boltz",
            "predict",
            container_input_fasta, # Input FASTA file
            "--output-dir", container_output_dir,
            "--use_msa_server", # Tell Boltz to use MSA server
            "--msa_only",       # Tell Boltz to *only* run MSA step
        ])
        
        # Add optional MSA server URL if provided
        msa_server_url = self.config.get("colabfold_msa_server_url")
        if msa_server_url:
            cmd.extend(["--msa_server_url", msa_server_url])
            logger.info(f"Using provided MSA server URL: {msa_server_url}")
        else:
             logger.info("Using Boltz internal default MSA server URL.")

        # 4. Execute command
        exit_code, stdout, stderr = self._run_command(cmd)

        # 5. Check results
        if exit_code != 0:
            logger.error(f"ColabFold/Boltz MSA pipeline failed with exit code {exit_code}.")
            return None
            
        # 6. Parse output directory for A3M files
        protein_id_to_a3m_path: Dict[str, str] = {}
        found_a3m = False
        try:
            for filename in os.listdir(colabfold_msa_output_dir):
                 # Assume format like <input_stem>_<chain_id>.a3m or just <chain_id>.a3m
                 # Let's assume Boltz names them based on FASTA header ID, e.g., A.a3m, B.a3m
                 if filename.endswith(".a3m"):
                     parts = filename.split('.a3m')[0]
                     # Simplistic assumption: filename is the chain ID
                     chain_id = parts 
                     # Check if this chain_id was actually in our input
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

