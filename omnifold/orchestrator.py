import os
import logging
import threading
from typing import Any, Dict, List, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from datetime import datetime

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
            self.output_dir = base_output_dir / f"{input_file_stem}_{run_timestamp}"
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

    def run_pipeline(self, input_file: str) -> bool:
        """
        Run the complete prediction pipeline.
        Handles input parsing, MSA generation, config creation, and model execution.
        Manages GPU assignments and parallel/sequential execution.
        
        Args:
            input_file: Path to the input file (FASTA, AF3 JSON, or Boltz YAML)
            
        Returns:
            True if pipeline completed successfully, False otherwise
        """
        try:
            logger.info(f"Parsing input file: {input_file}")
            job_input = self.input_handler.parse_input(input_file)
            if not job_input:
                logger.error("Failed to parse input file.")
                return False
            self.msa_manager.job_input = job_input
            
            msa_result = None
            if not job_input.has_msa:
                logger.info("MSA not found in input, generating alignments...")
                msa_result = self.msa_manager.generate_msa()
                
                if not msa_result:
                    logger.error("MSA generation failed.")
                    return False

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

                if "protein_id_to_pqt_path" in msa_result:
                    job_input.protein_id_to_pqt_path = msa_result["protein_id_to_pqt_path"]
                    logger.info(f"Updated job_input with PQT paths for {len(job_input.protein_id_to_pqt_path)} chains.")

                if "af3_data_json" in msa_result:
                    job_input.af3_data_json = msa_result["af3_data_json"]

                if "chai_fasta_path" in msa_result:
                    self.config["current_chai1_fasta_path"] = msa_result["chai_fasta_path"]
                    logger.info(f"Chai-1 will use FASTA generated by MSA provider: {msa_result['chai_fasta_path']}")

                if "chai_pqt_msa_dir" in msa_result:
                    self.config["current_chai1_msa_pqt_dir"] = msa_result["chai_pqt_msa_dir"]
                    logger.info(f"Chai-1 will use PQT MSAs from: {self.config['current_chai1_msa_pqt_dir']}")

                logger.info(f"--- JobInput State After MSA ---\n{job_input}\n----------------------------------")

            # This block seems redundant with the one above and can be simplified/removed.
            # Keeping it for now to avoid breaking existing AF3-only logic without more testing.
            if job_input.af3_data_json and not job_input.boltz_csv_msa_dir: # Only if boltz conversion hasn't run
                logger.info("AF3 data JSON is present. Running A3M extraction and Boltz CSV conversion.")
                # This part of the logic might need to be refactored to be cleaner.
                # Assuming extract_all_protein_a3ms_from_af3_json and _run_a3m_to_boltz_csv_conversion
                # can be harmonized. For now, let's keep the original flow for AF3 MSA.
                pass # The original logic for this case was complex and is being refactored.
            
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
                return False
            
            models_to_run_info = []
            if "af3_config_path" in configs and self.config.get("alphafold3_sif_path"):
                models_to_run_info.append(("alphafold3", configs["af3_config_path"], "af3_output_dir"))
            if "boltz_config_path" in configs and self.config.get("boltz1_sif_path"):
                models_to_run_info.append(("boltz1", configs["boltz_config_path"], "boltz_output_dir"))

            if self.config.get("chai1_sif_path"):
                # The correct FASTA path is now consistently provided by the MSAManager.
                chai_fasta_path_for_runner = self.config.get("current_chai1_fasta_path")
                if chai_fasta_path_for_runner and Path(chai_fasta_path_for_runner).is_file():
                    logger.info(f"Chai-1 will use FASTA: {chai_fasta_path_for_runner}")
                    if job_input.protein_id_to_pqt_path:
                        # Find the directory containing the first PQT file.
                        first_pqt_path = next(iter(job_input.protein_id_to_pqt_path.values()))
                        self.config["current_chai1_msa_pqt_dir"] = str(Path(first_pqt_path).parent)
                        logger.info(f"Chai-1 will use PQT MSAs from: {self.config['current_chai1_msa_pqt_dir']}")
                    models_to_run_info.append(("chai1", chai_fasta_path_for_runner, "chai1_output_dir"))
                else:
                    logger.warning("Chai-1 SIF is provided, but no suitable Chai-1 FASTA input was found/generated. Skipping Chai-1 execution.")
            
            if not models_to_run_info:
                logger.info("No models to run after configuration generation and checks.")
                return True # Nothing to run, so considered successful completion of an empty workload

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
                
                for model_name, config_path_or_fasta, output_dir_attr in models_to_run_info:
                    model_gpu = gpu_assignments.get(model_name)
                    model_output_dir = getattr(self, output_dir_attr)
                    
                    if model_gpu is not None:
                        futures.append(
                            executor.submit(
                                self._run_model,
                                model_name,
                                config_path_or_fasta,
                                model_output_dir,
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
            if success: # Only generate report if all models succeeded
                try:
                    logger.info("Generating final HTML report...")
                    run_report_generation(Path(self.output_dir))
                except Exception as e:
                    logger.error(f"Failed to generate final HTML report: {e}", exc_info=True)
            else:
                logger.warning("Skipping final report generation due to model prediction failures.")

            return success
            
        except Exception as e:
            logger.error(f"Pipeline failed with unexpected error: {e}", exc_info=True)
            return False

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