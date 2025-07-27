import os
import logging
import threading
import re
from typing import Any, Dict, List, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from datetime import datetime
import json

from .input_handler import InputHandler
from .msa_manager import MSAManager
from .config_generator import ConfigGenerator
from .runner import Runner
from .util.definitions import JobInput
from .util.gpu_utils import assign_gpus_to_models, set_gpu_visibility
from .util.msa_utils import extract_all_protein_a3ms_from_af3_json
from .util.af3_to_boltz_csv import convert_a3m_to_boltz_csv
from .html_report.generate_report import run_report_generation

logger = logging.getLogger(__name__)

class Orchestrator:
    """
    Central controller that orchestrates the entire protein ensemble prediction pipeline.
    Manages the flow from input processing through MSA generation, config creation,
    and parallel model execution.
    """

    def __init__(self, config: Dict[str, Any], output_dir: str):
        """
        Initialize the Orchestrator.
        
        Configures paths, initializes components, and sets up a unique timestamped output directory.

        Args:
            config: Global configuration dictionary containing paths and settings
            output_dir: Base output directory for all results
        """
        self.config = config
        base_output_dir = Path(output_dir)

        # Check if the base output directory contains specific model output folders to avoid conflicts
        model_dirs_exist = any(
            (base_output_dir / model_dir).exists() for model_dir in ["alphafold3", "boltz", "chai1"]
        )

        if model_dirs_exist:
            # If so, create a new timestamped subdirectory for this run
            run_timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
            input_file_stem = Path(config["input_file"]).stem
            # Sanitize the filename stem to remove problematic characters like ':'
            sanitized_stem = re.sub(r'[^a-zA-Z0-9_-]', '_', input_file_stem)
            self.output_dir = base_output_dir / f"{sanitized_stem}_{run_timestamp}"
        else:
            # Otherwise, use the user-provided directory as is
            self.output_dir = base_output_dir

        # Ensure the final output directory exists *before* setting up logging
        os.makedirs(self.output_dir, exist_ok=True)

        # Now that the final output directory is determined, setup logging
        self._setup_logging()

        # Log the directory decision
        if model_dirs_exist:
            logger.info(f"Target output directory contains previous results. Creating a new unique subdirectory for this run: {self.output_dir}")
        else:
            logger.info(f"Using specified output directory: {self.output_dir}")
        
        # Initialize components and create all necessary subdirectories
        self.input_handler = InputHandler()
        self.msa_output_dir = self.output_dir / "msa_generation"
        self.af3_output_dir = self.output_dir / "alphafold3"
        self.boltz_output_dir = self.output_dir / "boltz"
        self.chai1_output_dir = self.output_dir / "chai1"
        
        for path in [self.msa_output_dir, self.af3_output_dir, self.boltz_output_dir, self.chai1_output_dir]:
            os.makedirs(path, exist_ok=True)

        self.msa_manager = MSAManager(self.config, str(self.msa_output_dir))
        self.config_generator = ConfigGenerator()
        self.runner = Runner(self.config)
        self.job_state_file = self.output_dir / "omnifold_job.json"

    def _setup_logging(self):
        """
        Configure logging to write to both file and console.
        Sets up file handler with formatting and adds it to root logger if not already present.
        """
        log_file = os.path.join(self.output_dir, "ensemble_prediction.log")
        file_handler = logging.FileHandler(log_file)
        file_handler.setFormatter(logging.Formatter(
            '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
        ))
        root_logger = logging.getLogger()
        if not any(isinstance(h, logging.FileHandler) and h.baseFilename == log_file for h in root_logger.handlers):
            root_logger.addHandler(file_handler)
        root_logger.setLevel(logging.INFO)

    def _run_model(self, model_name: str, config_path: str, output_dir: str, gpu_id: int) -> Tuple[int, str, str]:
        """
        Run a single model prediction.
        Sets GPU visibility and executes model prediction.
        
        Args:
            model_name: Name of the model to run
            config_path: Path to the model's input config file
            output_dir: Directory to write model outputs
            gpu_id: GPU ID to use for this model
            
        Returns:
            Tuple of (exit_code, stdout, stderr)
        """
        logger.info(f"Starting {model_name} prediction on GPU {gpu_id}")
        set_gpu_visibility(gpu_id)
        return self.runner.run_model(
            model_name=model_name,
            input_config_file_host_path=config_path,
            model_specific_output_dir_host_path=output_dir,
            gpu_id=gpu_id
        )

    def _run_a3m_to_boltz_csv_conversion(self, protein_id_to_a3m_path: Dict[str, Any]) -> Optional[str]:
        """
        Converts A3M files to the Boltz CSV format.
        
        Args:
            protein_id_to_a3m_path: A dictionary mapping protein IDs to their A3M file paths.
                                    Can be flat {id: path} or nested {"unpaired": {id: path}}.
        
        Returns:
            The path to the output directory containing the CSV files, or None on failure.
        """
        boltz_csv_output_dir = self.msa_output_dir / "boltz_csv_msas"
        boltz_csv_output_dir.mkdir(exist_ok=True)
        
        # Accommodate both flat and nested dicts
        if not protein_id_to_a3m_path:
            logger.warning("No A3M files provided for Boltz CSV conversion.")
            return None

        # Pass the *full* mapping (including both 'paired' and 'unpaired' sections when present)
        logger.info(f"Starting A3M to Boltz CSV conversion. Output dir: {boltz_csv_output_dir}")
        try:
            convert_a3m_to_boltz_csv(
                protein_to_a3m_path=protein_id_to_a3m_path,
                output_csv_dir=str(boltz_csv_output_dir)
            )
            logger.info("A3M to Boltz CSV conversion completed successfully.")
            return str(boltz_csv_output_dir)
        except Exception as e:
            logger.error(f"Error during A3M to Boltz CSV conversion: {e}", exc_info=True)
            return None

    def _write_job_state(self, state: Dict[str, Any]):
        """Saves the job state (like config paths) to a JSON file."""
        try:
            with open(self.job_state_file, 'w') as f:
                json.dump(state, f, indent=4)
            logger.info(f"Job state saved to {self.job_state_file}")
        except IOError as e:
            logger.error(f"Failed to write job state to {self.job_state_file}: {e}", exc_info=True)

    def _read_job_state(self) -> Optional[Dict[str, Any]]:
        """Loads the job state from a JSON file."""
        if not self.job_state_file.is_file():
            logger.error(f"Job state file not found: {self.job_state_file}. Cannot proceed with --gpu_only run.")
            return None
        try:
            with open(self.job_state_file, 'r') as f:
                state = json.load(f)
            logger.info(f"Job state successfully loaded from {self.job_state_file}")
            return state
        except (IOError, json.JSONDecodeError) as e:
            logger.error(f"Failed to read or parse job state file {self.job_state_file}: {e}", exc_info=True)
            return None

    def run_pipeline(self, input_path: str, msa_only: bool = False, gpu_only: bool = False) -> bool:
        """
        Run the complete prediction pipeline.
        Handles input parsing, MSA generation, config creation, and model execution.
        Manages GPU assignments and parallel/sequential execution.
        
        Args:
            input_path: Path to the input file (FASTA, etc.) or input directory for gpu_only mode.
            msa_only: If True, stops after MSA and config generation.
            gpu_only: If True, skips to the GPU execution phase.
            
        Returns:
            True if pipeline completed successfully, False otherwise
        """
        if gpu_only:
            logger.info(f"--- Starting GPU-only phase from directory: {input_path} ---")
            job_state = self._read_job_state()
            if not job_state:
                return False
            
            # The 'configs' dict is the main artifact we need from the state file
            configs = job_state.get("configs", {})
            if not configs:
                logger.error("Job state file does not contain necessary 'configs' information.")
                return False
        else:
            logger.info(f"--- Starting new pipeline run for input: {input_path} ---")
            configs = self._generate_msas_and_configs(input_path)
            if not configs:
                logger.error("Failed to generate MSAs and configs.")
                return False

        if msa_only:
            job_state_to_save = {"configs": configs}
            # Persist extra artefacts outside configs if useful for quick access
            if configs.get("chai_pqt_msa_dir"):
                job_state_to_save["chai_pqt_msa_dir"] = configs["chai_pqt_msa_dir"]
            self._write_job_state(job_state_to_save)
            logger.info("MSA-only phase complete. Pipeline halted as requested.")
            return True

        # --- This point is the start of the GPU phase ---
        
        models_to_run_info = []
        if "af3_config_path" in configs and self.config.get("alphafold3_sif_path"):
            models_to_run_info.append(("alphafold3", configs["af3_config_path"], self.af3_output_dir))
        if "boltz_config_path" in configs and self.config.get("boltz1_sif_path"):
            models_to_run_info.append(("boltz1", configs["boltz_config_path"], self.boltz_output_dir))

        if "chai_config_path" in configs and self.config.get("chai1_sif_path"):
            chai_fasta_path = configs["chai_config_path"]
            if chai_fasta_path and Path(chai_fasta_path).is_file():
                logger.info(f"Chai-1 will use FASTA: {chai_fasta_path}")
                
                # In gpu_only mode, we must retrieve the PQT path from the job state.
                if gpu_only:
                    pqt_dir = configs.get("chai_pqt_msa_dir")
                    if pqt_dir and Path(pqt_dir).is_dir():
                        self.config["current_chai1_msa_pqt_dir"] = pqt_dir
                        logger.info(f"Loaded Chai-1 PQT MSA directory from job state: {pqt_dir}")
                    else:
                        logger.warning("Running Chai-1 in --gpu_only mode, but no PQT MSA directory found in job state. Chai-1 may fail if it requires local MSAs.")
                
                models_to_run_info.append(("chai1", chai_fasta_path, self.chai1_output_dir))
            else:
                logger.warning("Chai-1 SIF is provided, but no suitable Chai-1 FASTA input was found/generated. Skipping Chai-1 execution.")
        
        if not models_to_run_info:
            logger.info("No models to run after configuration generation and checks.")
            return True

        model_names_for_gpu_assignment = [info[0] for info in models_to_run_info]
        gpu_assignments = assign_gpus_to_models(model_names_for_gpu_assignment, force_sequential=self.config.get("run_sequentially", False))

        if not gpu_assignments:
            logger.error("Failed to assign GPUs to models. Aborting execution.")
            return False

        unique_gpu_ids = set(filter(None, gpu_assignments.values()))
        if self.config.get("run_sequentially", False) or len(unique_gpu_ids) <= 1 and len(models_to_run_info) > 1:
            max_workers = 1
            logger.info(f"Executing models sequentially with max_workers=1.")
        else:
            max_workers = len(unique_gpu_ids) if unique_gpu_ids else 1
            logger.info(f"Executing models potentially in parallel with max_workers={max_workers} (based on unique GPUs).")

        results = {}
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            
            for model_name, config_path_or_fasta, model_output_dir in models_to_run_info:
                model_gpu = gpu_assignments.get(model_name)
                
                if model_gpu is not None:
                    futures.append(
                        executor.submit(
                            self._run_model,
                            model_name,
                            config_path_or_fasta,
                            str(model_output_dir),
                            model_gpu
                        )
                    )
                else:
                    logger.warning(f"No GPU assigned for {model_name}, skipping execution.")

            if not futures:
                logger.error("No models were submitted for execution (e.g., due to no GPU assignment or no SIFs).")
                return False
            
            num_submitted = len(futures)
            processed_models_map = {id(future): models_to_run_info[i][0] for i, future in enumerate(futures)} 
            
            for future in as_completed(futures):
                model_name_completed = processed_models_map.get(id(future), f"unknown_model__{id(future)}")
                try:
                    result = future.result()
                    results[model_name_completed] = result
                    logger.info(f"Received result for {model_name_completed}")
                except Exception as e:
                    logger.error(f"Model execution for {model_name_completed} failed: {e}", exc_info=True)
                    if model_name_completed not in results:
                        results[model_name_completed] = (-1, "", f"Future failed with exception: {e}")
                finally:
                    if model_name_completed == "chai1" and "current_chai1_msa_pqt_dir" in self.config:
                        del self.config["current_chai1_msa_pqt_dir"]
                        logger.info("Cleaned up temporary Chai-1 PQT MSA directory path from config.")
        
        if len(results) != num_submitted:
            logger.error("Mismatch in submitted vs completed models. Some models may have failed or not reported results properly.")
            return False

        success = True
        for model_name, (exit_code, stdout, stderr) in results.items():
            if exit_code != 0:
                logger.error(f"{model_name} failed with exit code {exit_code}")
                logger.error(f"STDOUT:\n{stdout}")
                logger.error(f"STDERR:\n{stderr}")
                success = False
            else:
                logger.info(f"{model_name} completed successfully")
                logger.debug(f"{model_name} STDOUT:\n{stdout}")
        
        self._write_summary(results)
        
        # --- Generate Final HTML Report ---
        # Generate report if at least one model succeeded
        any_success = any(res[0] == 0 for res in results.values())
        if any_success:
            try:
                logger.info("Generating final HTML report for successful models...")
                run_report_generation(Path(self.output_dir))
            except Exception as e:
                logger.error(f"Failed to generate final HTML report: {e}", exc_info=True)
        else:
            logger.warning("Skipping final report generation as all model predictions failed.")

        # The overall pipeline success depends on ALL models succeeding.
        # This remains unchanged to correctly signal job status.
        return success

    def _generate_msas_and_configs(self, input_file: str) -> Optional[Dict[str, str]]:
        """
        Handles the first part of the pipeline: input parsing, MSA generation, and config creation.
        This function is called for a full run or an MSA-only run.
        """
        try:
            logger.info(f"Parsing input file: {input_file}")
            job_input = self.input_handler.parse_input(input_file)
            if not job_input:
                logger.error("Failed to parse input file.")
                return None
            self.msa_manager.job_input = job_input
            
            msa_result = None
            if not job_input.has_msa:
                logger.info("MSA not found in input, generating alignments...")
                msa_result = self.msa_manager.generate_msa()
                
                if not msa_result:
                    logger.error("MSA generation failed.")
                    return None

                # Existing handling when MSAManager provides detailed A3M mapping
                if "protein_id_to_a3m_path" in msa_result:
                    job_input.protein_id_to_a3m_path = msa_result["protein_id_to_a3m_path"]
                    logger.info("Updated job_input with MSA paths.")

                    # --- New: Convert A3Ms to Boltz CSV format ---
                    boltz_csv_dir = self._run_a3m_to_boltz_csv_conversion(job_input.protein_id_to_a3m_path)
                    if boltz_csv_dir:
                        job_input.boltz_csv_msa_dir = boltz_csv_dir
                        logger.info(f"Updated job_input with boltz_csv_msa_dir: {job_input.boltz_csv_msa_dir}")

                # If the MSAManager already produced Boltz-style CSV MSAs (e.g. AF3 pipeline), just propagate the path.
                elif "boltz_csv_msa_dir" in msa_result:
                    job_input.boltz_csv_msa_dir = msa_result["boltz_csv_msa_dir"]
                    logger.info(f"Using Boltz CSV MSAs generated by MSAManager at: {job_input.boltz_csv_msa_dir}")

                if "template_store_path" in msa_result and msa_result["template_store_path"]:
                    job_input.template_store_path = msa_result["template_store_path"]
                    logger.info(f"Templates generated by MSAManager at: {job_input.template_store_path}")
                    self.config["template_store_path"] = job_input.template_store_path

                if "protein_id_to_pqt_path" in msa_result:
                    job_input.protein_id_to_pqt_path = msa_result["protein_id_to_pqt_path"]
                    logger.info(f"Updated job_input with PQT paths for {len(job_input.protein_id_to_pqt_path)} chains.")

                if "af3_data_json" in msa_result:
                    job_input.af3_data_json = msa_result["af3_data_json"]

                if "chai_fasta_path" in msa_result:
                    # This path is now crucial for the job state.
                    self.config["current_chai1_fasta_path"] = msa_result["chai_fasta_path"]
                    logger.info(f"Chai-1 will use FASTA generated by MSA provider: {msa_result['chai_fasta_path']}")

                if "chai_pqt_msa_dir" in msa_result:
                    self.config["current_chai1_msa_pqt_dir"] = msa_result["chai_pqt_msa_dir"]
                    logger.info(f"Chai-1 will use PQT MSAs from: {self.config['current_chai1_msa_pqt_dir']}")

                logger.info(f"--- JobInput State After MSA ---\n{job_input}\n----------------------------------")

            if job_input.model_seeds is None and self.config.get("default_seed") is not None:
                default_seed_val = self.config.get("default_seed")
                job_input.model_seeds = [int(default_seed_val)]
                logger.info(f"Propagating CLI default_seed ({default_seed_val}) to job_input.model_seeds for ConfigGenerator.")
                if job_input.num_model_seeds_from_input is None:
                    job_input.num_model_seeds_from_input = 1 

            logger.info("Generating model configurations...")
            configs = self.config_generator.generate_configs(job_input, Path(self.output_dir), self.config)
            
            if not configs:
                logger.error("Failed to generate model configurations.")
                return None
            
            # Add the special Chai-1 FASTA path to the configs dict so it gets saved in the job state
            if self.config.get("current_chai1_fasta_path"):
                configs["chai_config_path"] = self.config.get("current_chai1_fasta_path")
            
            # Also save the PQT directory if it exists
            if self.config.get("current_chai1_msa_pqt_dir"):
                configs["chai_pqt_msa_dir"] = self.config.get("current_chai1_msa_pqt_dir")

            return configs

        except Exception as e:
            logger.error(f"MSA and config generation failed with unexpected error: {e}", exc_info=True)
            return None
            
    def _write_summary(self, results: Dict[str, Tuple[int, str, str]]):
        """
        Write a summary of the prediction results to a file.
        Includes execution status, error details if any, and output directory paths.
        
        Args:
            results: Dictionary mapping model names to their execution results
        """
        summary_path = os.path.join(self.output_dir, "prediction_summary.txt")
        try:
            with open(summary_path, "w") as f:
                f.write("OmniFold Prediction Summary\n")
                f.write("================================\n\n")
                
                model_order = ["alphafold3", "boltz1", "chai1"]
                for model_name in model_order:
                    if model_name in results:
                        exit_code, stdout, stderr = results[model_name]
                        f.write(f"{model_name.upper()}:\n")
                        f.write(f"  Status: {'Success' if exit_code == 0 else 'Failed'}\n")
                        f.write(f"  Exit Code: {exit_code}\n")
                        if exit_code != 0:
                            f.write(f"  Error Snippet (see logs for full details):\n  ---\n{stderr[:1000]}...\n  ---\n")
                        f.write("\n")
                
                f.write("\nOutput Directories:\n")
                f.write(f"AlphaFold3: {os.path.relpath(self.af3_output_dir, self.output_dir)}\n")
                f.write(f"Boltz-1: {os.path.relpath(self.boltz_output_dir, self.output_dir)}\n")
                f.write(f"Chai-1: {os.path.relpath(self.chai1_output_dir, self.output_dir)}\n")
            logger.info(f"Prediction summary written to: {summary_path}")
        except IOError as e:
            logger.error(f"Failed to write summary file {summary_path}: {e}")