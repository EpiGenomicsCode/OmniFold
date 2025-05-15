import os
import subprocess
import logging
from typing import Any, Dict, Optional, Tuple

logger = logging.getLogger(__name__)

class Runner:
    """
    Manages the execution of model prediction processes inside Singularity containers.
    """

    def __init__(self, config: Dict[str, Any]):
        """
        Initializes the Runner.

        Args:
            config: Global application configuration dictionary.
        """
        self.config = config
        logger.info("Runner initialized.")

    def _construct_base_singularity_cmd(self, sif_path: str, binds: Dict[str, str], gpu_id: Optional[int], use_run: bool = False) -> list[str]:
        """
        Constructs the base singularity command including common options and binds.

        Args:
            sif_path: Path to the Singularity image file.
            binds: Dictionary mapping host paths to container paths.
            gpu_id: Optional GPU ID to use.
            use_run: Whether to use 'singularity run' instead of 'singularity exec'.

        Returns:
            List of command components for the base singularity command.
        """
        cmd = ["singularity", "run" if use_run else "exec"]
        if gpu_id is not None:
            cmd.append("--nv")
        
        for host_path, container_path in binds.items():
            if os.path.exists(host_path):
                cmd.extend(["--bind", f"{host_path}:{container_path}"])
            else:
                logger.warning(f"Host path for binding does not exist, skipping bind: {host_path}")

        if not use_run:
            cmd.append(sif_path)
        else:
            # For 'run', the SIF path comes after all binds
            cmd.append(sif_path)
            
        return cmd

    def run_model(
        self,
        model_name: str,
        input_config_file_host_path: str,
        model_specific_output_dir_host_path: str,
        gpu_id: Optional[int] = None
    ) -> Tuple[int, str, str]:
        """
        Runs a specific model (AlphaFold3 or Boltz-1) inside its Singularity container.

        Args:
            model_name: Name of the model ("alphafold3" or "boltz1").
            input_config_file_host_path: Absolute path to the model-specific input config file (JSON/YAML) on the host.
            model_specific_output_dir_host_path: Absolute path to the dedicated output directory for this model on the host.
            gpu_id: The specific GPU ID to assign for this model run. If None, GPU is not specifically assigned.

        Returns:
            A tuple (exit_code, stdout, stderr).
        """
        logger.info(f"Preparing to run model: {model_name} with GPU ID: {gpu_id}")
        os.makedirs(model_specific_output_dir_host_path, exist_ok=True)

        sif_path = ""
        model_command = []
        binds = {}

        if model_name == "alphafold3":
            sif_path = self.config.get("alphafold3_sif_path")
            if not sif_path or not os.path.exists(sif_path):
                return -1, "", "AlphaFold3 SIF path not configured or not found."

            # Set up binds for AF3
            binds = {
                input_config_file_host_path: "/data/af_input/fold_input.json",
                model_specific_output_dir_host_path: "/data/af_output",
                self.config["alphafold3_model_weights_dir"]: "/data/models",
            }
            
            # Add database directory if provided
            if self.config.get("alphafold3_database_dir"):
                binds[self.config["alphafold3_database_dir"]] = "/data/public_databases"

            # Construct the run_alphafold.py command
            model_command = [
                "--json_path=/data/af_input/fold_input.json",
                "--model_dir=/data/models",
                "--output_dir=/data/af_output",
                "--run_data_pipeline=False",
                "--run_inference=True",
            ]

            # Add database directory if provided
            if self.config.get("alphafold3_database_dir"):
                model_command.append("--db_dir=/data/public_databases")

            # Add AlphaFold3 specific parameters from config if they are set
            if self.config.get("af3_num_recycles") is not None:
                model_command.extend(["--num_recycles", str(self.config["af3_num_recycles"])])
            if self.config.get("af3_num_diffusion_samples") is not None:
                model_command.extend(["--num_diffusion_samples", str(self.config["af3_num_diffusion_samples"])])
            if self.config.get("af3_num_seeds") is not None:
                model_command.extend(["--num_seeds", str(self.config["af3_num_seeds"])])
            if self.config.get("af3_save_embeddings"):
                model_command.append("--save_embeddings")
            if self.config.get("af3_max_template_date"):
                model_command.extend(["--max_template_date", self.config["af3_max_template_date"]])
            if self.config.get("af3_conformer_max_iterations") is not None:
                model_command.extend(["--conformer_max_iterations", str(self.config["af3_conformer_max_iterations"])])
            if self.config.get("af3_buckets"):
                for bucket_val in self.config["af3_buckets"]:
                    model_command.extend(["--buckets", str(bucket_val)])

        elif model_name == "boltz1":
            sif_path = self.config.get("boltz1_sif_path")
            if not sif_path or not os.path.exists(sif_path):
                return -1, "", "Boltz-1 SIF path not configured or not found."

            # Set up binds for Boltz
            binds = {
                input_config_file_host_path: input_config_file_host_path,  # Keep original path
                model_specific_output_dir_host_path: model_specific_output_dir_host_path,  # Keep original path
            }

            model_command = [
                "boltz", "predict",
                input_config_file_host_path,
                "--out_dir", model_specific_output_dir_host_path,
                "--recycling_steps", "3",
                "--sampling_steps", "200",
                "--diffusion_samples", "1",
                "--step_scale", "1.638",
                "--output_format", "mmcif"
            ]

            # Add Boltz-1 specific parameters from config
            if self.config.get("boltz_recycling_steps") is not None:
                model_command.extend(["--recycling_steps", str(self.config["boltz_recycling_steps"])])
            if self.config.get("boltz_sampling_steps") is not None:
                model_command.extend(["--sampling_steps", str(self.config["boltz_sampling_steps"])])
            if self.config.get("boltz_diffusion_samples") is not None:
                model_command.extend(["--diffusion_samples", str(self.config["boltz_diffusion_samples"])])
            if self.config.get("boltz_step_scale") is not None:
                model_command.extend(["--step_scale", str(self.config["boltz_step_scale"])])
            if self.config.get("boltz_no_potentials"):
                model_command.append("--no_potentials")
            if self.config.get("boltz_write_full_pae"):
                model_command.append("--write_full_pae")
            if self.config.get("boltz_write_full_pde"):
                model_command.append("--write_full_pde")
            if self.config.get("boltz_output_format"):
                model_command.extend(["--output_format", self.config["boltz_output_format"]])

            # Only use MSA server if explicitly configured
            if self.config.get("use_msa_server", False):
                model_command.append("--use_msa_server")
                if self.config.get("colabfold_msa_server_url"):
                    model_command.append("--msa_server_url")
                    model_command.append(self.config["colabfold_msa_server_url"])

        else:
            return -1, "", f"Unsupported model_name: {model_name}"

        # Construct the full singularity command
        # Use 'run' for AF3, 'exec' for Boltz
        use_run = (model_name == "alphafold3")
        singularity_base_cmd = self._construct_base_singularity_cmd(sif_path, binds, gpu_id, use_run)
        full_cmd = singularity_base_cmd + model_command

        logger.info(f"Executing command for {model_name}: {' '.join(full_cmd)}")

        env = os.environ.copy()
        if gpu_id is not None:
            env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
            logger.info(f"Setting CUDA_VISIBLE_DEVICES={gpu_id} for {model_name}")

        try:
            process = subprocess.run(
                full_cmd, 
                capture_output=True, 
                text=True, 
                check=False,
                env=env
            )
            logger.info(f"{model_name} execution finished with exit code: {process.returncode}")
            if process.stdout:
                logger.debug(f"{model_name} stdout:\n{process.stdout}")
            if process.stderr:
                if process.returncode == 0:
                    logger.debug(f"{model_name} stderr:\n{process.stderr}")
                else:
                    logger.error(f"{model_name} stderr:\n{process.stderr}")
            
            return process.returncode, process.stdout, process.stderr
        except FileNotFoundError:
            error_msg = "Singularity command not found. Ensure Singularity is installed and in PATH."
            logger.error(error_msg)
            return -1, "", error_msg
        except Exception as e:
            error_msg = f"An unexpected error occurred while trying to run {model_name}: {e}"
            logger.critical(error_msg, exc_info=True)
            return -1, "", str(e) 