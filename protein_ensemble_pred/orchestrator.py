import os
import logging
import threading
from typing import Any, Dict, List, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

from .input_handler import InputHandler
from .msa_manager import MSAManager
from .config_generator import ConfigGenerator
from .runner import Runner
from .util.gpu_utils import assign_gpus_to_models, set_gpu_visibility
from .util.msa_utils import extract_all_protein_a3ms_from_af3_json
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
        
        Configures paths, initializes components, and sets up output directories.

        Args:
            config: Global configuration dictionary containing paths and settings
            output_dir: Base output directory for all results
        """
        self.config = config
        self.output_dir = output_dir
        self.input_handler = InputHandler()
        self.msa_output_dir = os.path.join(output_dir, "msa_generation")
        os.makedirs(self.msa_output_dir, exist_ok=True)
        self.msa_manager = MSAManager(self.config, str(self.msa_output_dir))
        self.config_generator = ConfigGenerator()
        self.runner = Runner(self.config)

        os.makedirs(output_dir, exist_ok=True)
        self.af3_output_dir = os.path.join(output_dir, "alphafold3")
        self.boltz_output_dir = os.path.join(output_dir, "boltz")
        os.makedirs(self.af3_output_dir, exist_ok=True)
        os.makedirs(self.boltz_output_dir, exist_ok=True)

        self.chai1_output_dir = os.path.join(output_dir, "chai1")
        os.makedirs(self.chai1_output_dir, exist_ok=True)

        self._setup_logging()

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
                
                if msa_result.get("af3_data_json"):
                    job_input.af3_data_json = msa_result["af3_data_json"]
                
                if isinstance(msa_result.get("protein_id_to_json_unpaired_a3m_path"), dict):
                    job_input.protein_id_to_a3m_path = msa_result["protein_id_to_json_unpaired_a3m_path"]
                    logger.info(f"Updated job_input.protein_id_to_a3m_path with {len(job_input.protein_id_to_a3m_path)} JSON-extracted A3Ms.")
                elif isinstance(msa_result.get("protein_id_to_a3m_path"), dict):
                    job_input.protein_id_to_a3m_path = msa_result["protein_id_to_a3m_path"]
                    logger.info(f"Updated job_input.protein_id_to_a3m_path with {len(job_input.protein_id_to_a3m_path)} A3Ms from general MSA results.")
                elif msa_result.get("a3m_path") and not job_input.protein_id_to_a3m_path:
                    first_protein_id = job_input.sequences[0].chain_id if job_input.sequences else "unknown_protein_1"
                    job_input.protein_id_to_a3m_path = {first_protein_id: msa_result["a3m_path"]}
                    logger.info(f"Updated job_input.protein_id_to_a3m_path with single A3M: {msa_result['a3m_path']}")

                if msa_result.get("boltz_csv_msa_dir"):
                    job_input.boltz_csv_msa_dir = msa_result["boltz_csv_msa_dir"]
                    logger.info(f"Updated job_input with boltz_csv_msa_dir: {job_input.boltz_csv_msa_dir}")

                if msa_result.get("chai_pqt_msa_dir"):
                    self.config["current_chai1_msa_pqt_dir"] = msa_result["chai_pqt_msa_dir"]
                    logger.info(f"Chai-1 will use PQT MSA from AF3 MSA stage: {msa_result['chai_pqt_msa_dir']}")
                if msa_result.get("chai_fasta_path"):
                    logger.info(f"Chai-1 FASTA generated by MSA_Manager: {msa_result['chai_fasta_path']}")
            
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

            chai_fasta_path_for_runner = None
            if self.config.get("chai1_sif_path"):
                if msa_result and msa_result.get("chai_fasta_path"):
                    chai_fasta_path_for_runner = msa_result["chai_fasta_path"]
                    logger.info(f"Chai-1 will use FASTA generated from AF3 MSA stage: {chai_fasta_path_for_runner}")
                    if msa_result.get("chai_pqt_msa_dir"):
                        self.config["current_chai1_msa_pqt_dir"] = msa_result["chai_pqt_msa_dir"]
                        logger.info(f"Chai-1 will use PQT MSA from AF3 MSA stage: {msa_result['chai_pqt_msa_dir']}")
                    else:
                        logger.info("Chai-1: No PQT MSA directory found from AF3 MSA stage. Will rely on other MSA settings for Chai-1 (server or chai1_msa_directory).")
                    models_to_run_info.append(("chai1", chai_fasta_path_for_runner, "chai1_output_dir"))
                else:
                    logger.warning("Chai-1 SIF is provided, but no suitable Chai-1 FASTA input was found/generated from the AF3 MSA stage. Skipping Chai-1 execution.")
            
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
                f.write("Protein Ensemble Prediction Summary\n")
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