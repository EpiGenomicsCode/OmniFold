import argparse
import os
import logging
import sys

from .orchestrator import Orchestrator

def setup_logging(log_level_str: str):
    """Configures logging for the application."""
    numeric_level = getattr(logging, log_level_str.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f'Invalid log level: {log_level_str}')
    
    # Basic configuration ensuring all modules use this
    # Get the root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(numeric_level)

    # Remove any existing handlers to avoid duplicate logs if re-running in same session
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    # Add a new handler
    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    root_logger.addHandler(handler)
    
    # Get a logger for this module
    logger = logging.getLogger(__name__)
    logger.info(f"Logging configured at level: {log_level_str.upper()}")
    return logger


def main():
    parser = argparse.ArgumentParser(
        description="HPC Protein Ensemble Prediction CLI. \n"
                    "Runs protein structure prediction using an ensemble of models (AlphaFold3 and Boltz-1).",
        formatter_class=argparse.RawTextHelpFormatter  # Allows for better formatting of help text
    )

    # --- Input and Output Arguments ---
    parser.add_argument(
        "--input_file",
        required=True,
        type=str,
        help="Path to the input file (FASTA, AlphaFold3 JSON, or Boltz YAML)."
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        type=str,
        help="Path to the directory where results will be saved."
    )

    # --- Singularity Container Paths ---
    container_group = parser.add_argument_group('Singularity Container Paths')
    container_group.add_argument(
        "--alphafold3_sif_path",
        required=True,
        type=str,
        help="Path to the AlphaFold3 Singularity image (.sif file)."
    )
    container_group.add_argument(
        "--boltz1_sif_path",
        required=True,
        type=str,
        help="Path to the Boltz-1 Singularity image (.sif file)."
    )

    # --- Model Specific Paths ---
    model_paths_group = parser.add_argument_group('Model Specific Paths')
    model_paths_group.add_argument(
        "--alphafold3_model_weights_dir",
        required=True, # As per user: "they have to provide the model weights otherwise we cant run af3"
        type=str,
        help="Path to the directory containing AlphaFold3 model parameters/weights."
    )
    model_paths_group.add_argument(
        "--alphafold3_database_dir",
        type=str,
        help="Path to the root directory for AlphaFold3 databases (used for template search by AF3 itself, "
             "or if AF3 MSA pipeline is run directly for full data processing). "
             "MSAManager might use a more specific subset if it runs AF3 MSA only."
    )

    msa_group = parser.add_argument_group('MSA Generation Configuration')
    msa_group.add_argument(
        "--msa_method",
        choices=["alphafold3", "colabfold"],
        default="alphafold3",
        help="Method for MSA generation if not provided in input (default: alphafold3). "
             "If 'colabfold', local MMseqs2 execution is assumed unless --no_msa is also specified "
             "and Boltz is expected to use its own API/server."
    )
    msa_group.add_argument(
        "--no_msa",
        action="store_true",
        help="Skip MSA generation by this tool. Input must contain alignments, or models run MSA-free."
    )
    msa_group.add_argument(
        "--allow_msa_fallback",
        action="store_true",
        help="If the preferred MSA method (e.g., alphafold3) fails, attempt to use the other (e.g., colabfold local) as a fallback."
    )
    msa_group.add_argument(
        "--colabfold_msa_server_url",
        type=str,
        help="(Optional) URL for an external ColabFold MMseqs2 API server. If provided, and Boltz is configured to use it, "
             "this URL might be passed to the Boltz container. Not used by MSAManager's local colabfold/MMseqs2 run."
    )

    # --- Execution Control Arguments ---
    exec_group = parser.add_argument_group('Execution Control')
    exec_group.add_argument(
        "--sequential",
        action="store_true",
        help="Force sequential execution of models, even if multiple GPUs are available."
    )
    exec_group.add_argument(
        "--default_seed",
        type=int,
        default=42,
        help="Default random seed for stochastic parts of the pipeline (default: 42)."
    )
    exec_group.add_argument(
        "--log_level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        default="INFO",
        help="Set the logging level (default: INFO)."
    )

    # --- AlphaFold3 Specific Parameters ---
    af3_specific_group = parser.add_argument_group('AlphaFold3 Specific Parameters')
    af3_specific_group.add_argument(
        "--af3_num_recycles",
        type=int,
        default=10,
        help="Number of recycles for AlphaFold3 (default: 10)."
    )
    af3_specific_group.add_argument(
        "--af3_num_diffusion_samples",
        type=int,
        default=5,
        help="Number of diffusion samples for AlphaFold3 (default: 5)."
    )
    af3_specific_group.add_argument(
        "--af3_num_seeds",
        type=int,
        default=None,
        help="Number of seeds for AlphaFold3. If set, input JSON must provide a single seed. (default: None, uses seeds from input JSON)."
    )
    af3_specific_group.add_argument(
        "--af3_save_embeddings",
        action="store_true",
        help="Save final trunk single and pair embeddings in AlphaFold3 output (default: False)."
    )
    af3_specific_group.add_argument(
        "--af3_max_template_date",
        type=str,
        default="2021-09-30",
        help="Maximum template release date for AlphaFold3 (YYYY-MM-DD, default: 2021-09-30)."
    )
    af3_specific_group.add_argument(
        "--af3_conformer_max_iterations",
        type=int,
        default=None,
        help="Max iterations for RDKit conformer search in AlphaFold3 (default: RDKit default)."
    )
    af3_specific_group.add_argument(
        "--af3_buckets",
        type=str,
        nargs='+',
        default=None,
        help="(Optional) Token sizes for caching compilations in AlphaFold3. If not provided, AF3 defaults are used."
    )

    # --- Boltz-1 Specific Parameters ---
    boltz_specific_group = parser.add_argument_group('Boltz-1 Specific Parameters')
    boltz_specific_group.add_argument(
        "--boltz_recycling_steps",
        type=int,
        default=3,
        help="Number of recycling steps for Boltz-1 (default: 3)."
    )
    boltz_specific_group.add_argument(
        "--boltz_sampling_steps",
        type=int,
        default=200,
        help="Number of sampling steps for Boltz-1 (default: 200)."
    )
    boltz_specific_group.add_argument(
        "--boltz_diffusion_samples",
        type=int,
        default=1,
        help="Number of diffusion samples for Boltz-1 (default: 1)."
    )
    boltz_specific_group.add_argument(
        "--boltz_step_scale",
        type=float,
        default=1.638,
        help="Step size for diffusion process sampling (recommended between 1 and 2, default: 1.638)."
    )
    boltz_specific_group.add_argument(
        "--boltz_no_potentials",
        action="store_true",
        help="Disable potentials for steering in Boltz-1 (default: False)."
    )
    boltz_specific_group.add_argument(
        "--boltz_write_full_pae",
        action="store_true",
        help="Write full PAE (Predicted Aligned Error) to output (default: False)."
    )
    boltz_specific_group.add_argument(
        "--boltz_write_full_pde",
        action="store_true",
        help="Write full PDE (Predicted Distance Error) to output (default: False)."
    )
    boltz_specific_group.add_argument(
        "--boltz_output_format",
        choices=["pdb", "mmcif"],
        default="mmcif",
        help="Output format for Boltz-1 predictions (default: mmcif)."
    )
    
    args = parser.parse_args()
    logger = setup_logging(args.log_level)

    # Create the output directory if it doesn't exist
    try:
        os.makedirs(args.output_dir, exist_ok=True)
        logger.info(f"Ensured output directory exists: {args.output_dir}")
    except OSError as e:
        logger.error(f"Failed to create output directory {args.output_dir}: {e}")
        sys.exit(1) # Exit if we can\'t create the output directory

    # --- Prepare Configuration for Orchestrator ---
    # This config will be passed to all modules.
    config = {
        "input_file": os.path.abspath(args.input_file),
        "output_dir": os.path.abspath(args.output_dir),
        
        "alphafold3_sif_path": os.path.abspath(args.alphafold3_sif_path),
        "boltz1_sif_path": os.path.abspath(args.boltz1_sif_path),
        
        "alphafold3_model_weights_dir": os.path.abspath(args.alphafold3_model_weights_dir),
        "alphafold3_database_dir": os.path.abspath(args.alphafold3_database_dir) if args.alphafold3_database_dir else None,

        "msa_method_preference": args.msa_method,
        "skip_msa_generation": args.no_msa,
        "allow_msa_fallback": args.allow_msa_fallback,
        "colabfold_msa_server_url": args.colabfold_msa_server_url,

        "jackhmmer_binary_path": "jackhmmer",
        "nhmmer_binary_path": "nhmmer",
        "hmmsearch_binary_path": "hmmsearch",
        "hmmbuild_binary_path": "hmmbuild",
        "hmmalign_binary_path": "hmmalign",
        
        "run_sequentially": args.sequential,
        "default_seed": args.default_seed,
        "log_level": args.log_level.upper(),
        
        "alphafold_database_root_path": os.path.abspath(args.alphafold3_database_dir) if args.alphafold3_database_dir else None,
        "alphafold_model_params_path": os.path.abspath(args.alphafold3_model_weights_dir),

        # AlphaFold3 specific execution parameters
        "af3_num_recycles": args.af3_num_recycles,
        "af3_num_diffusion_samples": args.af3_num_diffusion_samples,
        "af3_num_seeds": args.af3_num_seeds,
        "af3_save_embeddings": args.af3_save_embeddings,
        "af3_max_template_date": args.af3_max_template_date,
        "af3_conformer_max_iterations": args.af3_conformer_max_iterations,

        # Boltz-1 specific execution parameters
        "boltz_recycling_steps": args.boltz_recycling_steps,
        "boltz_sampling_steps": args.boltz_sampling_steps,
        "boltz_diffusion_samples": args.boltz_diffusion_samples,
        "boltz_step_scale": args.boltz_step_scale,
        "boltz_no_potentials": args.boltz_no_potentials,
        "boltz_write_full_pae": args.boltz_write_full_pae,
        "boltz_write_full_pde": args.boltz_write_full_pde,
        "boltz_output_format": args.boltz_output_format,

        # GPU detection/assignment  needs to go here.
    }
    
    # Conditionally add af3_buckets to config
    if args.af3_buckets is not None:
        config["af3_buckets"] = args.af3_buckets
    
    # Add a check for required paths within the config that Orchestrator/submodules need
    required_paths = [
        "alphafold3_sif_path", 
        "boltz1_sif_path", 
        "alphafold3_model_weights_dir", 
        # alphafold3_database_dir is optional for prediction, required only for AF3 MSA gen
    ]
    missing_paths = [p for p in required_paths if not config.get(p) or not os.path.exists(config[p])]
    if missing_paths:
        for path_key in missing_paths:
            if not config.get(path_key):
                 logger.error(f"Missing required configuration parameter: --{path_key.replace('_sif_path','-sif-path').replace('_dir','-dir')}") # Attempt to match arg name
            else:
                 logger.error(f"Required path does not exist: {config[path_key]} (from --{path_key.replace('_sif_path','-sif-path').replace('_dir','-dir')})")
        sys.exit(1)

    logger.info("Configuration prepared. Initializing Orchestrator.")
    
    try:
        orchestrator = Orchestrator(config=config, output_dir=config["output_dir"])
        orchestrator.run_pipeline(input_file=config["input_file"])
        logger.info("Pipeline execution finished.")
    except Exception as e:
        logger.error(f"Pipeline execution failed: {e}", exc_info=True) # Log traceback
        sys.exit(1)


if __name__ == "__main__":
    main() 
