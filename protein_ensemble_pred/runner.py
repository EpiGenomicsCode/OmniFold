import os
import subprocess
import logging
from typing import Any, Dict, Optional, Tuple
from pathlib import Path

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

            binds = {
                input_config_file_host_path: "/data/af_input/fold_input.json",
                model_specific_output_dir_host_path: "/data/af_output",
                self.config["alphafold3_model_weights_dir"]: "/data/models",
            }
            
            if self.config.get("alphafold3_database_dir"):
                binds[self.config["alphafold3_database_dir"]] = "/data/public_databases"

            model_command = [
                "--json_path=/data/af_input/fold_input.json",
                "--model_dir=/data/models",
                "--output_dir=/data/af_output",
                "--run_data_pipeline=False",
                "--run_inference=True",
            ]

            if self.config.get("alphafold3_database_dir"):
                model_command.append("--db_dir=/data/public_databases")

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

            job_output_root_host_path = Path(input_config_file_host_path).parent.parent

            container_job_output_root = "/data/job_output"
            
            container_config_path = str(Path(container_job_output_root) / Path(input_config_file_host_path).relative_to(job_output_root_host_path))
            container_model_out_dir = str(Path(container_job_output_root) / Path(model_specific_output_dir_host_path).relative_to(job_output_root_host_path))

            binds = {
                str(job_output_root_host_path): container_job_output_root,
            }

            model_command = [
                "boltz", "predict",
                container_config_path,
                "--out_dir", container_model_out_dir,
            ]
            
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

            if self.config.get("use_msa_server", False):
                model_command.append("--use_msa_server")
                if self.config.get("colabfold_msa_server_url"):
                    model_command.append("--msa_server_url")
                    model_command.append(self.config["colabfold_msa_server_url"])

        elif model_name == "chai1":
            sif_path = self.config.get("chai1_sif_path")
            if not sif_path or not os.path.exists(sif_path):
                return -1, "", "Chai-1 SIF path not configured or not found."

            chai_fasta_host_path = input_config_file_host_path # Orchestrator passes FASTA path here
            if not os.path.exists(chai_fasta_host_path):
                return -1, "", f"Chai-1 input FASTA file not found: {chai_fasta_host_path}"

            container_fasta_path = "/data/input.fasta"
            container_output_dir_path = "/data/output"
            container_msa_dir_path = "/data/msas" 

            binds = {
                chai_fasta_host_path: container_fasta_path,
                model_specific_output_dir_host_path: container_output_dir_path,
            }

            model_command = ["chai-lab", "fold"]
            model_command.extend(["--fasta_file", container_fasta_path])
            model_command.extend(["--output_dir", container_output_dir_path])

            # MSA Directory Handling
            # Priority: 1. PQTs from AF3 MSA stage, 2. User-specified general MSA dir for Chai
            current_pqt_msa_dir = self.config.get("current_chai1_msa_pqt_dir")
            user_msa_dir = self.config.get("chai1_msa_directory")

            if not self.config.get("chai1_use_msa_server"):
                msa_dir_to_bind_host = None
                if current_pqt_msa_dir and os.path.isdir(current_pqt_msa_dir):
                    msa_dir_to_bind_host = current_pqt_msa_dir
                    logger.info(f"Chai-1: Using PQT MSA directory (from AF3 MSA): {msa_dir_to_bind_host}")
                elif user_msa_dir and os.path.isdir(user_msa_dir):
                    msa_dir_to_bind_host = user_msa_dir
                    logger.info(f"Chai-1: Using user-specified MSA directory: {msa_dir_to_bind_host}")
                
                if msa_dir_to_bind_host:
                    binds[msa_dir_to_bind_host] = container_msa_dir_path
                    model_command.extend(["--msa_directory", container_msa_dir_path])
                elif not current_pqt_msa_dir and not user_msa_dir:
                     logger.info("Chai-1: No local MSA directory provided and not using MSA server. Chai will run MSA-free or use internal defaults.")
                elif current_pqt_msa_dir and not os.path.isdir(current_pqt_msa_dir):
                    logger.warning(f"Chai-1: PQT MSA directory {current_pqt_msa_dir} not found. Will proceed without local MSAs.")
                elif user_msa_dir and not os.path.isdir(user_msa_dir):
                    logger.warning(f"Chai-1: User-specified MSA directory {user_msa_dir} not found. Will proceed without local MSAs.")
            else:
                logger.info("Chai-1: Configured to use MSA server.")
                if self.config.get("chai1_msa_server_url"):
                    model_command.extend(["--msa-server-url", self.config["chai1_msa_server_url"]])

            # Add other chai-lab fold arguments from config
            # Boolean flags (presence of flag means True, absence means False for chai-lab CLI typically)
            if self.config.get("chai1_use_msa_server"):
                 model_command.append("--use-msa-server") # Added above if logic allows, this ensures it if logic was bypassed by server priority
            if self.config.get("chai1_use_templates_server"):
                model_command.append("--use-templates-server")
            
            # Valued arguments
            # Note: For paths like template_hits_path and constraint_path, we are NOT binding them here.
            # If the user provides these to the CLI, they are passed as string arguments to chai-lab fold.
            # Chai-1/Singularity would need to handle access if these are host paths.
            # This Runner only explicitly binds the core input/output/MSA dirs.
            passthrough_args = {
                "chai1_template_hits_path": "--template_hits_path",
                "chai1_constraint_path": "--constraint_path",
                "chai1_msa_server_url": "--msa-server-url", # Already handled if use_msa_server is true, but safe to list
                "chai1_recycle_msa_subsample": "--recycle_msa_subsample",
                "chai1_num_trunk_recycles": "--num_trunk_recycles",
                "chai1_num_diffn_timesteps": "--num_diffn_timesteps",
                "chai1_num_diffn_samples": "--num_diffn_samples",
                "chai1_num_trunk_samples": "--num_trunk_samples",
                # chai1_seed and chai1_device are handled specially below due to defaults/gpu_id interaction
            }
            for config_key, cli_flag in passthrough_args.items():
                if self.config.get(config_key) is not None:
                    # Special handling for msa-server-url to avoid duplication if already added
                    if cli_flag == "--msa-server-url" and "--use-msa-server" in model_command:
                        if cli_flag not in model_command: # only add if not already there by the server block
                             model_command.extend([cli_flag, str(self.config[config_key])])
                    else:
                        model_command.extend([cli_flag, str(self.config[config_key])])

            # Seed: Use chai1_seed if provided, else orchestrator's default_seed
            seed_to_use = self.config.get("chai1_seed") if self.config.get("chai1_seed") is not None else self.config.get("default_seed")
            if seed_to_use is not None:
                model_command.extend(["--seed", str(seed_to_use)])

            # Device argument for Chai-1
            chai_device_arg = self.config.get("chai1_device")
            if chai_device_arg is None and gpu_id is not None: # If user didn't specify, but a GPU is assigned
                 chai_device_arg = "cuda:0" # Chai will see the assigned GPU as cuda:0 due to CUDA_VISIBLE_DEVICES
            if chai_device_arg:
                 model_command.extend(["--device", chai_device_arg])

            # Boolean flags that map to value strings for chai-lab run_inference function
            # Example: --use_esm_embeddings=False (if default is True in run_inference)
            # Assuming chai-lab fold CLI directly takes these boolean-like valued flags
            if not self.config.get("chai1_use_esm_embeddings", True):
                 model_command.append("--use_esm_embeddings=False")
            if not self.config.get("chai1_low_memory", True):
                 model_command.append("--low_memory=False")

        else:
            return -1, "", f"Unsupported model_name: {model_name}"

        use_run = (model_name == "alphafold3") # chai-lab fold is an argument to exec
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