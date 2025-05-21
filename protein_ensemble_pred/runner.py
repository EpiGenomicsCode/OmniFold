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

        elif model_name == "chai1":
            sif_path = self.config.get("chai1_sif_path")
            if not sif_path or not os.path.exists(sif_path):
                return -1, "", "Chai-1 SIF path not configured or not found."

            chai_fasta_host_path = input_config_file_host_path
            if not os.path.exists(chai_fasta_host_path):
                return -1, "", f"Chai-1 input FASTA file not found: {chai_fasta_host_path}"

            container_fasta_path = "/data/input.fasta"
            container_output_dir = "/data/output"
            container_msa_dir = "/data/msas"

            binds = {
                chai_fasta_host_path: f"{container_fasta_path}:ro",
                model_specific_output_dir_host_path: container_output_dir
            }

            model_command = [
                "chai-lab",
                "fold",
                container_fasta_path,
                container_output_dir,
            ]

            chai_pqt_msa_dir_host = self.config.get("current_chai1_msa_pqt_dir")
            user_msa_dir_host = self.config.get("chai1_msa_directory")

            if chai_pqt_msa_dir_host and os.path.isdir(chai_pqt_msa_dir_host):
                logger.info(f"Chai-1: Using PQT MSA directory (from AF3 MSA): {chai_pqt_msa_dir_host}")
                binds[chai_pqt_msa_dir_host] = f"{container_msa_dir}:ro"
                model_command.extend(["--msa-directory", container_msa_dir])
            elif user_msa_dir_host and os.path.isdir(user_msa_dir_host):
                logger.info(f"Chai-1: Using user-provided MSA directory: {user_msa_dir_host}")
                binds[user_msa_dir_host] = f"{container_msa_dir}:ro"
                model_command.extend(["--msa-directory", container_msa_dir])
            elif self.config.get("chai1_use_msa_server", True):
                logger.info("Chai-1: Using MSA server.")
                model_command.append("--use-msa-server")
                if self.config.get("colabfold_msa_server_url"):
                    model_command.extend(["--msa-server-url", self.config["colabfold_msa_server_url"]])
            else:
                logger.info("Chai-1: No local MSAs provided and MSA server is disabled.")

            optional_chai_params = {
                "recycle-msa-subsample": "chai1_recycle_msa_subsample",
                "num-trunk-recycles": "chai1_num_trunk_recycles",
                "num-diffn-timesteps": "chai1_num_diffn_timesteps",
                "num-diffn-samples": "chai1_num_diffn_samples",
                "num-trunk-samples": "chai1_num_trunk_samples",
                "seed": "chai1_seed",
                "device": "chai1_device",
                "use-templates-server": "chai1_use_templates_server",
            }

            for cli_opt, config_key in optional_chai_params.items():
                if config_key == "chai1_device":
                    continue
                if self.config.get(config_key) is not None:
                    val = self.config[config_key]
                    if isinstance(val, bool) and val:
                        model_command.append(f"--{cli_opt}")
                    elif not isinstance(val, bool):
                        model_command.extend([f"--{cli_opt}", str(val)])
            
            chai_device_arg = self.config.get("chai1_device")
            if gpu_id is not None and chai_device_arg is None:
                chai_device_arg = "cuda:0"
            
            if chai_device_arg:
                model_command.extend(["--device", chai_device_arg])

        else:
            return -1, "", f"Unsupported model_name: {model_name}"

        full_cmd = []
        if model_name == "alphafold3":
            full_cmd = ["singularity", "run"]
            if gpu_id is not None:
                full_cmd.append("--nv")
            for host_path, container_path_spec in binds.items():
                if os.path.exists(host_path):
                    full_cmd.extend(["--bind", f"{host_path}:{container_path_spec}"])
                else:
                    logger.warning(f"Host path for AF3 binding does not exist, skipping bind: {host_path}")
            full_cmd.append(sif_path)
            full_cmd.extend(model_command)
        
        else:
            full_cmd = ["singularity", "exec"]
            if gpu_id is not None:
                full_cmd.append("--nv")
            for host_path, container_path_spec in binds.items():
                if os.path.exists(host_path):
                    full_cmd.extend(["--bind", f"{host_path}:{container_path_spec}"])
                else:
                    logger.warning(f"Host path for {model_name} binding does not exist, skipping bind: {host_path}")
            full_cmd.append(sif_path)
            full_cmd.extend(model_command)

        logger.info(f"Final assembled command for {model_name}: {' '.join(full_cmd)}")

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