import os
import logging
import threading
from typing import Any, Dict, List, Optional, Tuple
from concurrent.futures import ThreadPoolExecutor, as_completed

from .input_handler import InputHandler
from .msa_manager import MSAManager
from .config_generator import ConfigGenerator
from .runner import Runner
from .util.gpu_utils import assign_gpus_to_models, set_gpu_visibility
from .util.msa_utils import extract_all_protein_a3ms_from_af3_json

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

        # Create output directory structure
        os.makedirs(output_dir, exist_ok=True)
        self.af3_output_dir = os.path.join(output_dir, "alphafold3")
        self.boltz_output_dir = os.path.join(output_dir, "boltz")
        os.makedirs(self.af3_output_dir, exist_ok=True)
        os.makedirs(self.boltz_output_dir, exist_ok=True)

        # Setup logging
        self._setup_logging()

    def _setup_logging(self):
        """Configure logging to write to both file and console."""
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
                if isinstance(msa_result.get("protein_id_to_a3m_path"), dict):
                    job_input.protein_id_to_a3m_path = msa_result["protein_id_to_a3m_path"]
                elif msa_result.get("a3m_path") and not job_input.protein_id_to_a3m_path:
                    first_protein_id = job_input.sequences[0].chain_id if job_input.sequences else "unknown_protein_1"
                    job_input.protein_id_to_a3m_path = {first_protein_id: msa_result["a3m_path"]}

            # Extract A3Ms for Boltz if an AF3 JSON is available (either from input or MSA step)
            # and we don't yet have the per-protein A3M path dictionary.
            if job_input.af3_data_json and not job_input.protein_id_to_a3m_path:
                logger.info("AF3 JSON found, extracting A3Ms for Boltz for each protein...")
                protein_id_to_a3m_path_dict = extract_all_protein_a3ms_from_af3_json(
                    job_input.af3_data_json, 
                    self.msa_output_dir 
                )
                if protein_id_to_a3m_path_dict is not None: # None indicates critical error
                    job_input.protein_id_to_a3m_path = protein_id_to_a3m_path_dict
                    logger.info(f"A3Ms for Boltz extracted: {len(protein_id_to_a3m_path_dict)} files.")
                    if not protein_id_to_a3m_path_dict:
                        logger.warning("No A3Ms were extracted from the AF3 JSON, though no critical error occurred (e.g., no proteins or no MSAs in JSON).")
                else:
                    logger.error("Critical error extracting A3Ms from AF3 JSON for Boltz. Cannot proceed.")
                    return False

            # For now, we assume ConfigGenerator will handle missing MSAs for specific proteins if Boltz allows it.
            if not job_input.protein_id_to_a3m_path and \
               not job_input.is_boltz_config and \
               not job_input.is_af3_msa_config_only:
                logger.warning("No per-protein A3M paths found or generated. Boltz prediction might fail or run MSA-free if not configured otherwise.")
                # Not returning False here, as Boltz might be run MSA-free or AF3 only is intended.

            logger.info("Generating model configurations...")
            configs = self.config_generator.generate_configs(job_input)
            if not configs:
                logger.error("Failed to generate model configurations.")
                return False
            
            gpu_assignments = assign_gpus_to_models(2)
            
            results = {}
            with ThreadPoolExecutor(max_workers=len(gpu_assignments)) as executor:
                futures = []
                
                if "af3_config_path" in configs:
                    af3_gpu = gpu_assignments.get("alphafold3")
                    if af3_gpu is not None:
                        futures.append(
                            executor.submit(
                                self._run_model,
                                "alphafold3",
                                configs["af3_config_path"],
                                self.af3_output_dir,
                                af3_gpu
                            )
                        )
                    else:
                        logger.warning("No GPU assigned for AlphaFold3, skipping execution.")
                else:
                    logger.info("No AlphaFold3 configuration generated, skipping execution.")

                if "boltz_config_path" in configs:
                    boltz_gpu = gpu_assignments.get("boltz1")
                    if boltz_gpu is not None:
                        futures.append(
                            executor.submit(
                                self._run_model,
                                "boltz1",
                                configs["boltz_config_path"],
                                self.boltz_output_dir,
                                boltz_gpu
                            )
                        )
                    else:
                        logger.warning("No GPU assigned for Boltz-1, skipping execution.")
                else:
                    logger.info("No Boltz-1 configuration generated, skipping execution.")

                if not futures:
                    logger.error("No models were submitted for execution.")
                    return False
                
                num_submitted = len(futures)
                processed_models = set()
                for future in as_completed(futures):
                    try:
                        model_name_order = [m for m in ["alphafold3", "boltz1"] if m in gpu_assignments and m not in processed_models]
                        model_name = model_name_order[0] if model_name_order else f"unknown_model_{len(processed_models)}"                        
                        result = future.result()
                        results[model_name] = result
                        processed_models.add(model_name)
                        logger.info(f"Received result for {model_name}")
                    except Exception as e:
                        logger.error(f"A model execution failed: {e}", exc_info=True)
                        return False
            
            if len(results) != num_submitted:
                logger.error("Mismatch in submitted vs completed models. Some models may have failed silently.")
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
            
            return success
            
        except Exception as e:
            logger.error(f"Pipeline failed with unexpected error: {e}", exc_info=True)
            return False

    def _write_summary(self, results: Dict[str, Tuple[int, str, str]]):
        """
        Write a summary of the prediction results.
        
        Args:
            results: Dictionary mapping model names to their execution results
        """
        summary_path = os.path.join(self.output_dir, "prediction_summary.txt")
        try:
            with open(summary_path, "w") as f:
                f.write("Protein Ensemble Prediction Summary\n")
                f.write("================================\n\n")
                
                model_order = ["alphafold3", "boltz1"]
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
            logger.info(f"Prediction summary written to: {summary_path}")
        except IOError as e:
            logger.error(f"Failed to write summary file {summary_path}: {e}")