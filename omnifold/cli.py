import argparse
import os
import logging
import sys
from pathlib import Path

from .orchestrator import Orchestrator


def setup_logging(log_level_str: str):
    """
    Configures logging for the application.
    
    Basic configuration ensuring all modules use this by configuring the root logger.
    Removes any existing handlers to avoid duplicate logs if re-running in same session.
    """
    numeric_level = getattr(logging, log_level_str.upper(), None)
    if not isinstance(numeric_level, int):
        raise ValueError(f'Invalid log level: {log_level_str}')
    
    root_logger = logging.getLogger()
    root_logger.setLevel(numeric_level)

    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)

    handler = logging.StreamHandler(sys.stdout)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    root_logger.addHandler(handler)
    
    logger = logging.getLogger(__name__)
    logger.info(f"Logging configured at level: {log_level_str.upper()}")
    return logger


def main():
    """
    OmniFold: An HPC Protein Ensemble Prediction CLI.
    Runs protein structure prediction using an ensemble of models (AlphaFold3, Boltz, and Chai-1).
    """
    parser = argparse.ArgumentParser(
        description="OmniFold: An HPC Protein Ensemble Prediction CLI that runs protein structure prediction "
                    "using an ensemble of models (AlphaFold3, Boltz, and Chai-1).",
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument(
        "--input_file",
        required=False,  # Changed to False as it's not always required
        type=str,
        help="Path to the input file (FASTA, AlphaFold3 JSON, or Boltz YAML)."
    )
    parser.add_argument(
        "--output_dir",
        required=True,
        type=str,
        help="Path to the directory where results will be saved."
    )

    pipeline_group = parser.add_argument_group('Pipeline Control')
    mut_group = pipeline_group.add_mutually_exclusive_group()
    mut_group.add_argument(
        "--msa_only",
        action="store_true",
        help="Run only the MSA and config generation phase. The pipeline will stop before model execution. "
             "The specified --output_dir can then be used as --input_dir for a --gpu_only run."
    )
    mut_group.add_argument(
        "--gpu_only",
        action="store_true",
        help="Run only the model execution (GPU) phase. Provide an --output_dir that contains the results of a previous --msa_only run."
    )
    mut_group.add_argument(
        "--report_only",
        action="store_true",
        help="Generate/repair the HTML report only. Provide --output_dir pointing to a completed job directory."
    )

    container_group = parser.add_argument_group('Singularity Container Paths')
    container_group.add_argument(
        "--alphafold3_sif_path",
        type=str,
        default=None,
        help="Path to the AlphaFold3 Singularity image (.sif file). Required if running AlphaFold3."
    )
    container_group.add_argument(
        "--boltz1_sif_path",
        type=str,
        default=None,
        help="Path to the Boltz Singularity image (.sif file). Required if running Boltz."
    )
    container_group.add_argument(
        "--chai1_sif_path",
        type=str, 
        default=None,
        help="Path to the Chai-1 Singularity image (.sif file). If provided, Chai-1 will be run."
    )

    model_paths_group = parser.add_argument_group('Model Specific Paths')
    model_paths_group.add_argument(
        "--alphafold3_model_weights_dir",
        type=str,
        default=None,
        help="Path to the directory containing AlphaFold3 model parameters/weights. Required if running AlphaFold3."
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

    exec_group = parser.add_argument_group('Execution Control')
    exec_group.add_argument(
        "--no_report",
        action="store_true",
        help="Skip HTML report generation even if models succeed. Useful for debugging or when running on constrained filesystems."
    )

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

    boltz_specific_group = parser.add_argument_group('Boltz Specific Parameters')
    boltz_specific_group.add_argument(
        "--boltz_recycling_steps",
        type=int,
        help="Number of recycling steps for Boltz (default: 3)."
    )
    boltz_specific_group.add_argument(
        "--boltz_sampling_steps",
        type=int,
        help="Number of sampling steps for Boltz (default: 200)."
    )
    boltz_specific_group.add_argument(
        "--boltz_diffusion_samples",
        type=int,
        help="Number of diffusion samples for Boltz (default: 1)."
    )
    boltz_specific_group.add_argument(
        "--boltz_step_scale",
        type=float,
        help="Step size for diffusion process sampling (recommended between 1 and 2, default: 1.638)."
    )
    boltz_specific_group.add_argument(
        "--boltz_no_potentials",
        action="store_true",
        help="Disable potentials for steering in Boltz (default: False)."
    )
    boltz_specific_group.add_argument(
        "--boltz_write_full_pae",
        action="store_true",
        default=True,
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
        help="Output format for Boltz predictions (default: mmcif)."
    )

    chai1_specific_group = parser.add_argument_group('Chai-1 Specific Parameters')
    chai1_specific_group.add_argument(
        "--chai1_use_msa_server",
        action="store_true",
        help="Use MSA server for Chai-1 (default: False)."
    )
    chai1_specific_group.add_argument(
        "--chai1_use_templates_server",
        action="store_true",
        help="Use templates server for Chai-1 (default: False)."
    )
    chai1_specific_group.add_argument(
        "--chai1_disable_esm_embeddings",
        action="store_true",
        help="Disable ESM embeddings for Chai-1 (default: ESM embeddings are enabled)."
    )
    chai1_specific_group.add_argument(
        "--chai1_msa_server_url",
        type=str,
        default="https://api.colabfold.com",
        help="URL for ColabFold MSA server used by Chai-1 (default: https://api.colabfold.com)."
    )
    chai1_specific_group.add_argument(
        "--chai1_msa_directory",
        type=str,
        default=None,
        help="Path to precomputed MSAs for Chai-1 (optional)."
    )
    chai1_specific_group.add_argument(
        "--chai1_constraint_path",
        type=str,
        default=None,
        help="Path to constraints file for Chai-1 (optional)."
    )
    chai1_specific_group.add_argument(
        "--chai1_template_hits_path",
        type=str,
        default=None,
        help="Path to template hits file (e.g., .m8) for Chai-1 (optional)."
    )
    chai1_specific_group.add_argument(
        "--chai1_recycle_msa_subsample",
        type=int, 
        default=None, # Default for a parameter that might not be set
        help="Recycle MSA subsample size for Chai-1 (e.g., 64). Default: Chai-1 internal."
    )
    chai1_specific_group.add_argument(
        "--chai1_num_trunk_recycles",
        type=int, 
        default=None, 
        help="Number of trunk recycles for Chai-1. Default: Chai-1 internal."
    )
    chai1_specific_group.add_argument(
        "--chai1_num_diffn_timesteps", 
        type=int, 
        default=None, 
        help="Number of diffusion timesteps for Chai-1. Default: Chai-1 internal."
    )
    chai1_specific_group.add_argument(
        "--chai1_num_diffn_samples", 
        type=int, 
        default=None, 
        help="Number of diffusion samples for Chai-1. Default: Chai-1 internal."
    )
    chai1_specific_group.add_argument(
        "--chai1_num_trunk_samples", 
        type=int, 
        default=None, 
        help="Number of trunk samples for Chai-1. Default: Chai-1 internal."
    )
    chai1_specific_group.add_argument(
        "--chai1_seed", 
        type=int, 
        default=None, 
        help="Random seed for Chai-1. Default: Chai-1 internal."
    )
    chai1_specific_group.add_argument(
        "--chai1_device",
        type=str,
        default=None, # Corresponds to 'cuda:0' in chai1.py if None
        help="Device for Chai-1 (e.g., 'cuda:0', 'cpu'; default handled by Chai-1: 'cuda:0')."
    )
    chai1_specific_group.add_argument(
        "--chai1_disable_low_memory",
        action="store_true",
        help="Disable low memory mode for Chai-1 (default: low memory mode is enabled)."
    )
    
    # Store defaults before parsing
    defaults = {}
    for action in parser._actions:
        # Filter out help actions and other non-argument actions
        if action.dest != "help" and hasattr(action, 'default'):
            defaults[action.dest] = action.default

    args = parser.parse_args()
    logger = setup_logging(args.log_level)

    # Shortcut for report-only mode
    if args.report_only:
        from omnifold.html_report.generate_report import run_report_generation
        try:
            run_report_generation(Path(os.path.abspath(args.output_dir)))
            logger.info("Report generation completed successfully.")
            return
        except Exception as e:
            logger.error(f"Report generation failed: {e}", exc_info=True)
            sys.exit(1)

    # --- Validation for new flags ---
    if args.report_only:
        if not args.output_dir:
            parser.error("--output_dir is required for --report_only mode.")
    elif args.gpu_only:
        if not args.output_dir:
            parser.error("--output_dir pointing to a previous --msa_only run is required for --gpu_only mode.")
        if not os.path.isdir(args.output_dir):
            parser.error(f"--output_dir specified for --gpu_only is not a valid directory: {args.output_dir}")
        if not os.path.exists(os.path.join(args.output_dir, "omnifold_job.json")):
             parser.error(f"omnifold_job.json not found in {args.output_dir}. Please ensure the directory is a valid output from an --msa_only run.")
        if args.input_file:
            logger.warning("--input_file is ignored in --gpu_only mode as it uses the previous MSA run's output directory.")
    elif args.msa_only:
        if not args.input_file:
            parser.error("--msa_only requires --input_file to be specified.")
        if not args.output_dir:
            parser.error("--output_dir is required for --msa_only run.")
    else: # Full run
        if not args.input_file:
            parser.error("--input_file is required for a full pipeline run.")
        if not args.output_dir:
            parser.error("--output_dir is required for a full run.")
    # ---

    config = vars(args) # Convert parsed args to a dictionary

    # --- Add _is_user_specified flags for Chai-1 --- 
    # These correspond to keys in optional_chai_params in runner.py (excluding 'device')
    chai1_config_keys_map = {
        "chai1_recycle_msa_subsample": "recycle-msa-subsample",
        "chai1_num_trunk_recycles": "num-trunk-recycles",
        "chai1_num_diffn_timesteps": "num-diffn-timesteps",
        "chai1_num_diffn_samples": "num-diffn-samples",
        "chai1_num_trunk_samples": "num-trunk-samples",
        "chai1_seed": "seed",
        "chai1_use_templates_server": "use-templates-server" # This is an action=store_true flag
    }

    for arg_dest, _ in chai1_config_keys_map.items():
        parsed_value = getattr(args, arg_dest, None)
        default_value = defaults.get(arg_dest)
        action = next((act for act in parser._actions if act.dest == arg_dest), None)

        user_specified = False
        if isinstance(action, argparse._StoreTrueAction):
            # For store_true, it's user-specified if the value is True (meaning flag was present)
            # and its default was False (or None, though typically False for store_true).
            # If default is True, it's user-specified if value is False (meaning flag was NOT present, 
            # but this scenario is less common for positive flags).
            # Simpler: if value is not the default for a boolean action.
            if parsed_value != default_value:
                user_specified = True
        elif isinstance(action, argparse._StoreFalseAction):
            if parsed_value != default_value:
                user_specified = True
        elif parsed_value is not None and parsed_value != default_value:
            # For other types, if it's not None and not the default, it was user-specified.
            # This also handles cases where default is None and user provides a value.
            user_specified = True
        elif default_value is None and parsed_value is not None:
            # Explicitly handles case where default is None and user provides any value.
            user_specified = True
        
        config[f"{arg_dest}_is_user_specified"] = user_specified
        if user_specified:
            logger.debug(f"Chai-1 arg '{arg_dest}' was user-specified with value: {parsed_value}")
    # ---

    config = {
        "input_file": os.path.abspath(args.input_file) if args.input_file else None,
        "output_dir": os.path.abspath(args.output_dir),
        
        "alphafold3_sif_path": os.path.abspath(args.alphafold3_sif_path) if args.alphafold3_sif_path else None,
        "boltz1_sif_path": os.path.abspath(args.boltz1_sif_path) if args.boltz1_sif_path else None,
        "chai1_sif_path": os.path.abspath(args.chai1_sif_path) if args.chai1_sif_path else None,
        
        "alphafold3_model_weights_dir": os.path.abspath(args.alphafold3_model_weights_dir) if args.alphafold3_model_weights_dir else None,
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
        "generate_report": not args.no_report,
        "log_level": args.log_level.upper(),
        
        "alphafold_database_root_path": os.path.abspath(args.alphafold3_database_dir) if args.alphafold3_database_dir else None,
        "alphafold_model_params_path": os.path.abspath(args.alphafold3_model_weights_dir) if args.alphafold3_model_weights_dir else None,

        "af3_num_recycles": args.af3_num_recycles,
        "af3_num_diffusion_samples": args.af3_num_diffusion_samples,
        "af3_num_seeds": args.af3_num_seeds,
        "af3_save_embeddings": args.af3_save_embeddings,
        "af3_max_template_date": args.af3_max_template_date,
        "af3_conformer_max_iterations": args.af3_conformer_max_iterations,

        "boltz_recycling_steps": args.boltz_recycling_steps,
        "boltz_sampling_steps": args.boltz_sampling_steps,
        "boltz_diffusion_samples": args.boltz_diffusion_samples,
        "boltz_diffusion_samples_is_default": parser.get_default("boltz_diffusion_samples") == args.boltz_diffusion_samples,
        "boltz_step_scale": args.boltz_step_scale,
        "boltz_no_potentials": args.boltz_no_potentials,
        "boltz_write_full_pae": args.boltz_write_full_pae,
        "boltz_write_full_pde": args.boltz_write_full_pde,
        "boltz_output_format": args.boltz_output_format,

        "chai1_use_msa_server": args.chai1_use_msa_server,
        "chai1_use_templates_server": args.chai1_use_templates_server,
        "chai1_use_esm_embeddings": not args.chai1_disable_esm_embeddings,
        "chai1_msa_server_url": args.chai1_msa_server_url,
        "chai1_msa_directory": os.path.abspath(args.chai1_msa_directory) if args.chai1_msa_directory else None,
        "chai1_constraint_path": os.path.abspath(args.chai1_constraint_path) if args.chai1_constraint_path else None,
        "chai1_template_hits_path": os.path.abspath(args.chai1_template_hits_path) if args.chai1_template_hits_path else None,
        "chai1_recycle_msa_subsample": args.chai1_recycle_msa_subsample,
        "chai1_num_trunk_recycles": args.chai1_num_trunk_recycles,
        "chai1_num_diffn_timesteps": args.chai1_num_diffn_timesteps,
        "chai1_num_diffn_samples": args.chai1_num_diffn_samples,
        "chai1_num_trunk_samples": args.chai1_num_trunk_samples,
        "chai1_seed": args.chai1_seed,
        "chai1_device": args.chai1_device,
        "chai1_low_memory": not args.chai1_disable_low_memory,
    }
    
    if args.af3_buckets is not None:
        config["af3_buckets"] = args.af3_buckets
    
    # --- Path and Model Sanity Checks ---
    msa_only_colabfold = args.msa_only and args.msa_method == "colabfold"

    sifs_provided = [
        config["alphafold3_sif_path"],
        config["boltz1_sif_path"],
        config["chai1_sif_path"]
    ]

    # Enforce at least one SIF **unless** we are in an MSA-only run that will use ColabFold (no container needed)
    if not any(sifs_provided) and not msa_only_colabfold:
        logger.error("No Singularity image file (.sif) provided. "
                     "Please provide at least one SIF path for AlphaFold3, Boltz, or Chai-1 "
                     "using --alphafold3_sif_path, --boltz1_sif_path, or --chai1_sif_path.")
        sys.exit(1)

    error_messages = []

    # Check AlphaFold3 paths only if a SIF is needed / provided
    if config["alphafold3_sif_path"]:
        if not os.path.exists(config["alphafold3_sif_path"]):
            error_messages.append(f"AlphaFold3 SIF path does not exist: {config['alphafold3_sif_path']} (from --alphafold3_sif_path)")
        if not config["alphafold3_model_weights_dir"]:
            error_messages.append("--alphafold3_model_weights_dir is required when --alphafold3_sif_path is provided.")
        elif not os.path.exists(config["alphafold3_model_weights_dir"]):
            error_messages.append(f"AlphaFold3 model weights directory does not exist: {config['alphafold3_model_weights_dir']} (from --alphafold3_model_weights_dir)")
    else:
        # In MSA-only mode with alphafold3 pipeline we MUST have the SIF and weights
        if args.msa_only and args.msa_method == "alphafold3":
            error_messages.append("--msa_only with msa_method=alphafold3 requires --alphafold3_sif_path to be provided.")

    # Check Boltz path if SIF is provided and we are not in msa_only_colabfold mode
    if config["boltz1_sif_path"] and not os.path.exists(config["boltz1_sif_path"]):
        error_messages.append(f"Boltz SIF path does not exist: {config['boltz1_sif_path']} (from --boltz1_sif_path)")

    # Check Chai-1 path if SIF is provided and not in msa_only_colabfold mode
    if config["chai1_sif_path"] and not os.path.exists(config["chai1_sif_path"]):
        error_messages.append(f"Chai-1 SIF path does not exist: {config['chai1_sif_path']} (from --chai1_sif_path)")

    if error_messages:
        for msg in error_messages:
            logger.error(msg)
        sys.exit(1)

    logger.info("Configuration prepared. Initializing Orchestrator.")
    
    try:
        orchestrator = Orchestrator(config=config, output_dir=config["output_dir"])
        input_path = config["output_dir"] if args.gpu_only else config["input_file"]
        orchestrator.run_pipeline(
            input_path=input_path, 
            msa_only=args.msa_only, 
            gpu_only=args.gpu_only
        )
        logger.info("Pipeline execution finished.")
    except Exception as e:
        logger.error(f"Pipeline execution failed: {e}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main() 
