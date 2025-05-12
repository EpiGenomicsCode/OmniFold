import json
import os
import subprocess
import logging
from typing import Any, Dict, Optional, Tuple

# Configure logging
logger = logging.getLogger(__name__)


class MSAManager:
    """
    Manages the generation of Multiple Sequence Alignments (MSAs) using
    either the AlphaFold 3 pipeline or a ColabFold-like (MMseqs2) pipeline.
    """

    def __init__(self, config: Dict[str, Any], job_input: Dict[str, Any], output_dir: str):
        """
        Initializes the MSA_Manager.

        Args:
            config: A dictionary containing configuration for MSA tools,
                    container paths, database paths, etc.
            job_input: Standardized input data from the InputHandler.
                       Expected to contain sequences and a flag indicating
                       if MSA is already provided.
            output_dir: The main output directory for the job.
        """
        self.config = config
        self.job_input = job_input
        self.output_dir = output_dir
        self.msa_tmp_dir = os.path.join(self.output_dir, "msa_intermediate")
        os.makedirs(self.msa_tmp_dir, exist_ok=True)
        logger.info("MSA_Manager initialized.")

    def generate_msa(self, msa_method: str = "alphafold3") -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
        """
        Generates MSA based on the chosen method.

        Args:
            msa_method: The method to use for MSA generation.
                        Can be "alphafold3" or "colabfold".

        Returns:
            A tuple containing:
            - msa_results: Dictionary containing MSA data (e.g., path to AF3 data JSON or A3M file content).
                           Specific structure depends on the method.
            - error_message: A string containing an error message if MSA generation failed, otherwise None.
        """
        if self.job_input.get("alignment_provided"):
            logger.info("MSA data found in input. Skipping MSA generation.")
            # TODO: Process provided alignment if necessary (e.g., format conversion)
            # For now, assume it's directly usable or handled by ConfigGenerator
            return self.job_input.get("alignment_data"), None

        logger.info(f"MSA generation requested using method: {msa_method}")

        if msa_method == "alphafold3":
            return self._run_alphafold3_msa_pipeline()
        elif msa_method == "colabfold":
            return self._run_colabfold_msa_pipeline()
        else:
            error_msg = f"Unsupported MSA method: {msa_method}"
            logger.error(error_msg)
            return None, error_msg

    def _run_alphafold3_msa_pipeline(self) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
        """
        Runs the AlphaFold 3 data pipeline to generate MSAs.
        This involves running the AF3 singularity container with flags to only
        perform data preparation.
        """
        logger.info("Starting AlphaFold 3 MSA pipeline.")
        af3_container_path = self.config.get("alphafold3_container_sif")
        if not af3_container_path or not os.path.exists(af3_container_path):
            return None, "AlphaFold 3 container SIF path not configured or not found."

        # Prepare a minimal AF3 JSON input for the data pipeline
        # This typically contains just the sequences.
        # The structure of job_input.sequences is assumed to be a list of dicts,
        # e.g., [{"id": "A", "sequence": "SEQ", "type": "protein"}, ...]
        
        # Create a unique name for this MSA run to avoid conflicts if multiple chains
        # or if the main job_input.name is too generic.
        input_name = self.job_input.get("name", "msa_job")
        if not self.job_input.get("sequences"):
             return None, "No sequences found in job_input for AF3 MSA pipeline."

        # For simplicity, we'll assume a single sequence or the first sequence for now.
        # A more robust implementation would handle multiple chains for complexes,
        # potentially creating separate AF3 inputs or a combined one.
        
        first_chain = self.job_input["sequences"][0]
        chain_id = first_chain.get("id", "A")
        sequence = first_chain.get("sequence")
        molecule_type = first_chain.get("type", "protein").lower() # protein, rna, dna

        if not sequence:
            return None, f"Sequence for chain {chain_id} is missing."

        minimal_af3_input = {
            "name": f"{input_name}_{chain_id}_msa_gen",
            "sequences": [
                {molecule_type: {"id": chain_id, "sequence": sequence}}
            ],
            # Ensure dialect and version are compatible with AF3 data pipeline
            "dialect": "alphafold3",
            "version": 3 
        }
        
        # Potentially add other required fields if AF3 data pipeline needs them,
        # e.g., modelSeeds (even if just a dummy one for data pipeline).
        # From run_alphafold.py, modelSeeds is a top-level list.
        minimal_af3_input["modelSeeds"] = [self.config.get("default_seed", 42)]


        input_json_path = os.path.join(self.msa_tmp_dir, f"{minimal_af3_input['name']}_input.json")
        with open(input_json_path, 'w') as f:
            json.dump(minimal_af3_input, f)
        
        logger.info(f"Minimal AF3 input for MSA generation written to: {input_json_path}")

        # Define output directory for this specific MSA run
        af3_msa_output_dir = os.path.join(self.msa_tmp_dir, f"{minimal_af3_input['name']}_output")
        os.makedirs(af3_msa_output_dir, exist_ok=True)

        # Construct the singularity command
        # Based on PRD: singularity exec --bind /db/path:/db --bind $PWD/tmp_out:/app/output \
        # af3_container.sif af3 run --json_path input.json --output_dir /app/output --norun_inference
        # We need to adapt paths and flags.
        # The run_alphafold.py script has flags like --run_data_pipeline=True and --run_inference=False
        
        cmd = [
            "singularity", "exec",
            "--nv", # As per PRD, for GPU resources, though MSA might be CPU
        ]

        # Database bindings - these need to come from self.config
        # Example: self.config.get("alphafold_database_paths") -> dict of db_name: host_path
        # The container might expect them at specific locations or via env vars.
        # For now, let's assume a general data dir binding if provided
        # run_alphafold.py uses DB_DIR which can be multiple paths.
        # It also defines specific database paths like _UNIREF90_DATABASE_PATH
        
        # Simplified binding for critical paths; a full setup would bind all necessary DBs.
        # This needs to be robustly configured.
        db_dir_host = self.config.get("alphafold_database_root_path") # e.g., /path/to/all_dbs
        if db_dir_host and os.path.exists(db_dir_host):
             # This assumes the AF3 container knows to look inside /db for its constituent DBs
             # or that environment variables inside the container are set accordingly.
            cmd.extend(["-B", f"{db_dir_host}:/db"]) # A common convention
        else:
            logger.warning("AlphaFold database root path not configured or not found. MSA quality may be affected.")
            # Potentially return error if essential DBs are missing for AF3 MSA

        # Binding for input JSON and output directory
        cmd.extend(["-B", f"{os.path.abspath(input_json_path)}:/app/input.json:ro"])
        cmd.extend(["-B", f"{os.path.abspath(af3_msa_output_dir)}:/app/output"])
        
        cmd.append(af3_container_path)
        
        # Command inside singularity: run_alphafold.py or af3 tool
        # Assuming run_alphafold.py is the entry point or accessible
        cmd.extend([
            "python3", "run_alphafold.py", # Or the correct entrypoint for af3
            "--json_path", "/app/input.json",
            "--output_dir", "/app/output",
            "--model_dir", self.config.get("alphafold_model_params_path", "/app/models"), # Path inside or bound
            "--run_data_pipeline=True",
            "--run_inference=False", # Crucial for MSA only
            # Add paths to binaries like jackhmmer, hhsearch etc. if container needs them
            # These should be configured in self.config and passed if not in container's PATH
            "--jackhmmer_binary_path", self.config.get("jackhmmer_binary_path", "jackhmmer"),
            "--nhmmer_binary_path", self.config.get("nhmmer_binary_path", "nhmmer"),
            "--hmmsearch_binary_path", self.config.get("hmmsearch_binary_path", "hmmsearch"),
            "--hmmbuild_binary_path", self.config.get("hmmbuild_binary_path", "hmmbuild"),
            # Add database paths as per run_alphafold.py flags, pointing to /db/...
            # This is complex due to many DBs. A simpler way might be if AF3 container
            # is pre-configured to look for DBs in /db and we bind the root DB dir.
            # Or, set environment variables inside singularity if AF3 respects them.
            # For now, omitting explicit DB paths flags assuming /db binding is sufficient
            # and AF3 knows how to find them, or these are defaults in AF3 run script.
            # Referencing run_alphafold.py, it has many --*_database_path flags.
            # These would need to be populated from self.config and point to paths *inside* the container.
            # Example: --uniref90_database_path /db/uniref90_2022_05.fa (if uniref90 is in the root of db_dir_host)
        ])
        
        # Add DB flags if properly configured
        # This part is crucial and complex. Needs careful mapping from host paths to container paths
        # and ensuring run_alphafold.py gets the correct arguments.
        # Example:
        # if db_dir_host: # And specific DBs are configured
        #    cmd.extend(["--uniref90_database_path", f"/db/{os.path.basename(self.config.get('uniref90_path_on_host'))}"])
        #    ... and so on for all DBs listed in run_alphafold.py ...
        # This requires self.config to have paths like 'uniref90_path_on_host'.

        logger.info(f"Executing AlphaFold 3 MSA command: {' '.join(cmd)}")
        try:
            process = subprocess.run(cmd, capture_output=True, text=True, check=True)
            logger.info("AlphaFold 3 MSA pipeline completed successfully.")
            logger.debug(f"AF3 MSA stdout: {process.stdout}")
            
            # Expected output is <name>_data.json in the af3_msa_output_dir
            # The name is minimal_af3_input['name']
            output_data_json_path = os.path.join(af3_msa_output_dir, f"{minimal_af3_input['name']}_data.json")
            
            if os.path.exists(output_data_json_path):
                logger.info(f"AlphaFold 3 MSA data JSON found at: {output_data_json_path}")
                # The result for now is the path to this file.
                # ConfigGenerator will parse it.
                return {"af3_data_json_path": output_data_json_path}, None
            else:
                error_msg = f"AlphaFold 3 MSA output data JSON not found at {output_data_json_path}."
                logger.error(error_msg)
                logger.error(f"AF3 MSA stderr: {process.stderr}")
                return None, error_msg

        except subprocess.CalledProcessError as e:
            error_msg = f"AlphaFold 3 MSA pipeline failed. Error: {e.stderr}"
            logger.error(error_msg)
            return None, error_msg
        except FileNotFoundError:
            error_msg = "Singularity command not found. Ensure Singularity is installed and in PATH."
            logger.error(error_msg)
            return None, error_msg

    def _run_colabfold_msa_pipeline(self) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
        """
        Runs a ColabFold-like (MMseqs2) MSA pipeline.
        This could involve calling a local MMseqs2 pipeline or Boltz's
        functionality if it supports local MSA generation to A3M.
        """
        logger.info("Starting ColabFold (MMseqs2) MSA pipeline.")
        # This is a placeholder. Implementation depends on how MMseqs2/Boltz is invoked.
        # PRD: "use a ColabFold MMseqs2-based search to quickly gather MSAs.
        # This could involve calling a local MMseqs2 pipeline (as in ColabFold)
        # or using Boltz's ability to query an MMseqs2 API."
        # "on an offline cluster, an equivalent local search will be performed"
        # Output should be an A3M file.

        # Example: If we have a script that runs local MMseqs2
        # mmseqs_runner_script = self.config.get("mmseqs_runner_script")
        # input_fasta_path = ... (create FASTA from self.job_input.sequences)
        # output_a3m_path = os.path.join(self.msa_tmp_dir, f"{self.job_input.get('name', 'msa_job')}.a3m")
        # cmd = [mmseqs_runner_script, "-i", input_fasta_path, "-o", output_a3m_path,
        #        "--db-path", self.config.get("mmseqs_database_path")]
        #
        # try:
        #    subprocess.run(cmd, check=True, capture_output=True, text=True)
        #    logger.info(f"MMseqs2 MSA generation successful. A3M at: {output_a3m_path}")
        #    return {"a3m_file_path": output_a3m_path}, None
        # except subprocess.CalledProcessError as e:
        #    error_msg = f"ColabFold/MMseqs2 MSA pipeline failed. Error: {e.stderr}"
        #    logger.error(error_msg)
        #    return None, error_msg
        
        # For now, returning not implemented
        logger.warning("ColabFold/MMseqs2 MSA pipeline is not fully implemented yet.")
        return None, "ColabFold/MMseqs2 MSA pipeline not implemented."

    def get_msa_results(self) -> Tuple[Optional[Dict[str, Any]], Optional[str]]:
        """
        Orchestrates MSA generation based on user preference or defaults.
        This method could eventually allow choosing between AF3 MSA or ColabFold MSA
        based on a flag from the main CLI.
        """
        # Defaulting to AlphaFold3 MSA for now as per PRD preference for AF3 pipeline.
        # This choice could be a parameter or configured.
        # cli.py may parse --msa-method [af3|mmseqs]
        msa_method_preference = self.config.get("msa_method_preference", "alphafold3")
        
        if "alignment_data" in self.job_input and self.job_input.get("alignment_provided"):
            logger.info("Using pre-existing alignment data from input.")
            return self.job_input["alignment_data"], None

        logger.info(f"Attempting MSA generation with preferred method: {msa_method_preference}")
        msa_results, error = self.generate_msa(msa_method=msa_method_preference)

        if error:
            logger.error(f"MSA generation with {msa_method_preference} failed: {error}")
            # Potential fallback strategy: if AF3's search fails, perhaps try the MMseqs2 pipeline.
            # PRD: "A potential fallback strategy: if AF3's search fails, perhaps try the MMseqs2 pipeline"
            if msa_method_preference == "alphafold3" and self.config.get("allow_msa_fallback", False):
                logger.info("Attempting fallback MSA generation with colabfold.")
                msa_results, error = self.generate_msa(msa_method="colabfold")
                if error:
                    logger.error(f"Fallback MSA generation with colabfold also failed: {error}")
                    return None, f"All MSA generation attempts failed. Main: {msa_method_preference}, Fallback: colabfold. Last error: {error}"
            elif msa_method_preference == "colabfold" and self.config.get("allow_msa_fallback", False):
                logger.info("Attempting fallback MSA generation with alphafold3.")
                msa_results, error = self.generate_msa(msa_method="alphafold3")
                if error:
                    logger.error(f"Fallback MSA generation with alphafold3 also failed: {error}")
                    return None, f"All MSA generation attempts failed. Main: {msa_method_preference}, Fallback: alphafold3. Last error: {error}"
            else:
                return None, f"MSA generation with {msa_method_preference} failed: {error}. No fallback attempted or fallback also failed."
        
        if msa_results:
            logger.info(f"MSA generation successful using {msa_method_preference if not error else 'fallback method'}.") # Small logic fix here
            return msa_results, None
        else:
            # This case should ideally be covered by the error conditions above.
            return None, "MSA generation did not produce results for an unknown reason."


if __name__ == '__main__':
    # Example Usage (requires a dummy config and job_input)
    logging.basicConfig(level=logging.INFO)
    
    # Create dummy directories and files for testing
    os.makedirs("/tmp/protein_ensemble_pred_test/output/msa_intermediate", exist_ok=True)
    os.makedirs("/tmp/protein_ensemble_pred_test/output/msa_intermediate/msa_job_A_msa_gen_output", exist_ok=True)
    
    # Dummy AF3 container (replace with actual path or skip if not testing execution)
    # DUMMY_AF3_SIF = "/path/to/your/dummy_alphafold3.sif" # Needs to exist
    # with open(DUMMY_AF3_SIF, "w") as f: f.write("#!/bin/bash\necho 'Fake AF3 SIF'")
    # os.chmod(DUMMY_AF3_SIF, 0o755)

    # Dummy run_alphafold.py script that the container might call
    # This is highly dependent on the actual container setup.
    # For a real test, this script would need to simulate run_alphafold.py behavior
    # and create the expected *_data.json file.

    # For this example, we'll assume the container execution might fail or skip it
    # by not providing a valid SIF path, so we test other logic.

    dummy_config = {
        "alphafold3_container_sif": "non_existent_dummy_alphafold3.sif", # Set to a real SIF for actual run
        "alphafold_database_root_path": "/tmp/dummy_dbs", # Needs to exist if used
        "alphafold_model_params_path": "/tmp/dummy_models", # Needs to exist if used
        "jackhmmer_binary_path": "jackhmmer", # Mock or ensure in PATH
        "nhmmer_binary_path": "nhmmer",
        "hmmsearch_binary_path": "hmmsearch",
        "hmmbuild_binary_path": "hmmbuild",
        "default_seed": 1234,
        "msa_method_preference": "alphafold3", # or "colabfold"
        "allow_msa_fallback": True,
        # Add other paths like uniref90_path_on_host etc. if fully testing AF3 cmd
    }
    os.makedirs(dummy_config["alphafold_database_root_path"], exist_ok=True)
    os.makedirs(dummy_config["alphafold_model_params_path"], exist_ok=True)

    dummy_job_input = {
        "name": "test_protein",
        "sequences": [
            {"id": "A", "sequence": "MADEINHEAVEN", "type": "protein"}
        ],
        "alignment_provided": False, # Set to True to test skipping
        # "alignment_data": {"a3m_file_path": "/path/to/existing.a3m"} # If alignment_provided is True
    }

    msa_manager = MSAManager(
        config=dummy_config,
        job_input=dummy_job_input,
        output_dir="/tmp/protein_ensemble_pred_test/output"
    )

    msa_data, error_msg = msa_manager.get_msa_results()

    if error_msg:
        logger.error(f"MSA Generation failed: {error_msg}")
    else:
        logger.info(f"MSA Generation successful. Results: {msa_data}")
        # Example: if AF3 MSA ran, msa_data might be {'af3_data_json_path': '...'}
        # If ColabFold ran, msa_data might be {'a3m_file_path': '...'}

    # To simulate a successful AF3 data pipeline run for subsequent components:
    # Create a dummy *_data.json if the SIF was non_existent
    if dummy_config["alphafold3_container_sif"] == "non_existent_dummy_alphafold3.sif" and not error_msg:
        # This part is tricky because we don't know the exact output name without running the _run_alphafold3_msa_pipeline
        # Let's assume it would have been:
        # /tmp/protein_ensemble_pred_test/output/msa_intermediate/test_protein_A_msa_gen_output/test_protein_A_msa_gen_data.json
        
        # Simplified: if the function returned a path, we'd use that.
        # For testing, you might manually create this dummy file:
        dummy_data_json_content = {
            "name": "test_protein_A_msa_gen",
            "sequences": [{"protein": {"id": "A", "sequence": "MADEINHEAVEN"}}],
            "msas": [{"protein_A": ["seq1", "seq2"]}], # Simplified MSA structure
            "templates": []
        }
        # Construct the expected path if the msa_data contains af3_data_json_path
        if msa_data and "af3_data_json_path" in msa_data:
             # This would only happen if the dummy SIF somehow "worked" or if the path was returned
             # despite the SIF not existing (which shouldn't happen with current error checks).
             # For robust offline testing of components *after* MSA manager,
             # one would typically mock the MSA manager's output or prepare test data files.
             pass
        elif not msa_data and not error_msg: # If it "succeeded" but returned None (e.g. ColabFold not implemented)
            logger.info("MSA generation path might have succeeded but returned no data (e.g., not implemented branch).")


    # Clean up dummy files/dirs if necessary
    # shutil.rmtree("/tmp/protein_ensemble_pred_test")
    # if os.path.exists(DUMMY_AF3_SIF): os.remove(DUMMY_AF3_SIF)
    # shutil.rmtree(dummy_config["alphafold_database_root_path"])
    # shutil.rmtree(dummy_config["alphafold_model_params_path"]) 