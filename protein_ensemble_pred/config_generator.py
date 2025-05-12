import json
import os
import random
from typing import List, Optional, Dict, Any
import yaml
from pathlib import Path

from .util.definitions import JobInput, SequenceInfo, SequenceType
from .af3_models import (
    Af3Input, Protein, ProteinChain, RNA, RNAChain, DNA, DNAChain, Ligand, LigandMolecule,
    MolId, ProtSeq, RNASeq, DNASeq # Type aliases for validation
)

class ConfigGenerator:
    def __init__(self):
        """Initializes the ConfigGenerator."""
        pass

    def generate_af3_json_from_job_input(
        self, 
        job_input: JobInput,
        output_dir: Path,
        # msa_results: Optional[Dict[str, str]] = None # Future: from MSAManager {chain_id: msa_path}
    ) -> Path:
        """
        Generates an AlphaFold3 JSON input file from a JobInput object.
        Prioritizes MSA paths from job_input.input_msa_paths if available.
        Later, will integrate with msa_results from MSAManager.

        Args:
            job_input: The JobInput object from InputHandler.
            output_dir: Directory to save the generated JSON file.
            # msa_results: Optional dictionary mapping chain_id to its generated MSA file path.

        Returns:
            Path to the generated JSON file.
        Raises:
            ValueError: If JobInput is empty or contains no sequences, or no valid sequences processed.
            IOError: If there's an error writing the JSON file.
        """
        if not job_input or not job_input.sequences:
            raise ValueError("JobInput is empty or contains no sequences.")

        output_dir.mkdir(parents=True, exist_ok=True)
        json_file_path = output_dir / f"{job_input.name_stem}_af3_generated.json"

        try:
            af3_sequences = []
            for seq_info in job_input.sequences:
                chain_id = seq_info.chain_id
                msa_path_for_chain = None

                # Priority 1: MSA path from the input AF3 JSON itself
                if job_input.input_msa_paths and chain_id in job_input.input_msa_paths:
                    msa_path_for_chain = job_input.input_msa_paths[chain_id]
                
                # Priority 2: MSA path from MSAManager (Future)
                # elif msa_results and chain_id in msa_results:
                #     msa_path_for_chain = msa_results[chain_id]
                
                common_chain_args = {
                    "id": chain_id,
                    "sequence": seq_info.sequence # Will be cast to ProtSeq, RNASeq, DNASeq by Pydantic
                }
                if msa_path_for_chain: # Only add if a path is available
                    common_chain_args["unpairedMsaPath"] = msa_path_for_chain
                    # Note: pairedMsaPath and templates are not handled yet from JobInput

                if seq_info.molecule_type == "protein":
                    protein_chain = ProteinChain(**common_chain_args)
                    af3_sequences.append(Protein(protein=protein_chain))
                elif seq_info.molecule_type == "rna":
                    rna_chain = RNAChain(**common_chain_args)
                    af3_sequences.append(RNA(rna=rna_chain))
                elif seq_info.molecule_type == "dna":
                    dna_chain = DNAChain(**common_chain_args)
                    af3_sequences.append(DNA(dna=dna_chain))
                elif seq_info.molecule_type == "ligand_ccd":
                    ligand_mol = LigandMolecule(id=chain_id, ccdCodes=[seq_info.sequence])
                    af3_sequences.append(Ligand(ligand=ligand_mol))
                elif seq_info.molecule_type == "ligand_smiles":
                    ligand_mol = LigandMolecule(id=chain_id, smiles=seq_info.sequence)
                    af3_sequences.append(Ligand(ligand=ligand_mol))
                else: # unknown
                    print(f"Warning: Skipping sequence {seq_info.original_name} (ID: {chain_id}) due to unknown type: {seq_info.molecule_type}")
                    continue
            
            if not af3_sequences:
                raise ValueError("No valid sequences could be processed for AF3 JSON.")

            model_seeds = [random.randint(1, 100000)]

            af3_input_data = Af3Input(
                name=job_input.name_stem,
                modelSeeds=model_seeds,
                sequences=af3_sequences,
                # bondedAtomPairs will be None unless specifically added from somewhere else
            )

            json_string = af3_input_data.model_dump_json(indent=2, by_alias=True, exclude_none=True)
            
            with open(json_file_path, "w") as f:
                f.write(json_string)
            
            print(f"AlphaFold3 JSON config generated at: {json_file_path}")

            # If the original input was an AF3 JSON and it was essentially just copied/re-serialized 
            # (e.g. it already contained MSA paths and we just used them),
            # we might not need a separate "_generated.json" file in that specific scenario.
            # For now, always creating a new "_generated.json" for clarity of what this function produces.
            # The Orchestrator can later decide if it needs to use an original input JSON directly.
            if job_input.raw_input_type == "af3_json":
                original_json_path = output_dir / f"{job_input.name_stem}.json" # Assuming original might be needed
                # This part is tricky: if the input JSON was data-rich (e.g. _data.json from AF3 pipeline)
                # we might want to use it directly. The current logic re-builds it. 
                # PRD: "AF3 JSON ... the handler will extract the sequences (and any provided alignments or templates) 
                # for use by the other model's pipeline as well."
                # PRD: "If the input is already in one model's advanced format...the handler will extract sequences...for use by other model's pipeline"
                # This implies our generated JSON from this function IS the standardized one for AF3 run.
                pass # Current behavior is to always generate a new one.

            return json_file_path

        except ValueError as ve: # Catch specific ValueError from pre-checks
            raise ve 
        except Exception as e: # Broad catch for Pydantic errors, file IO, etc.
            raise IOError(f"Error generating AlphaFold3 JSON for {job_input.name_stem} at {json_file_path}: {e}")

    def generate_boltz_yaml_from_job_input(
        self,
        job_input: JobInput,
        output_dir: Path
    ) -> Path:
        """
        Generates a Boltz-1 YAML configuration file from a JobInput object.
        The YAML will be saved in the specified output_dir.
        Returns the path to the generated YAML file.
        """
        boltz_config = {"version": 1, "sequences": []}
        sequences_list = []

        for seq_info in job_input.sequences:
            entity = {}
            if seq_info.molecule_type == "protein":
                protein_data = {
                    "id": seq_info.chain_id,
                    "sequence": seq_info.sequence
                }
                # Check if an MSA path was provided in the input for this specific chain
                if seq_info.chain_id in job_input.input_msa_paths:
                    protein_data["msa"] = job_input.input_msa_paths[seq_info.chain_id]
                # else: Boltz will attempt to generate MSA if not provided or run single sequence
                entity["protein"] = protein_data
            elif seq_info.molecule_type == "rna":
                entity["rna"] = {
                    "id": seq_info.chain_id,
                    "sequence": seq_info.sequence
                }
            elif seq_info.molecule_type == "dna":
                entity["dna"] = {
                    "id": seq_info.chain_id,
                    "sequence": seq_info.sequence
                }
            elif seq_info.molecule_type == "ligand_smiles":
                entity["ligand"] = {
                    "id": seq_info.chain_id,
                    "smiles": seq_info.sequence
                }
            elif seq_info.molecule_type == "ligand_ccd":
                entity["ligand"] = {
                    "id": seq_info.chain_id,
                    "ccd": seq_info.sequence
                }
            elif seq_info.molecule_type == "unknown":
                # PRD: "Unknown" sequence type handled with a warning by ConfigGenerator
                # For Boltz, it's unclear how it handles "unknown". 
                # We will skip adding it to the Boltz YAML and log a warning.
                print(f"Warning: Sequence '{seq_info.original_name}' (chain {seq_info.chain_id}) has an unknown molecule type and will be skipped for Boltz YAML generation.")
                continue # Skip this sequence for Boltz
            else:
                # Should not happen if SequenceType is strictly followed
                print(f"Warning: Unhandled molecule type '{seq_info.molecule_type}' for sequence '{seq_info.original_name}'. Skipping for Boltz YAML.")
                continue
            
            sequences_list.append(entity)

        if not sequences_list:
            raise ValueError("No suitable sequences found in JobInput to generate Boltz YAML.")

        boltz_config["sequences"] = sequences_list

        if job_input.constraints:
            boltz_config["constraints"] = job_input.constraints
        
        # Ensure output directory exists
        output_dir.mkdir(parents=True, exist_ok=True)
        
        output_filename = f"{job_input.name_stem}_boltz_generated.yaml"
        output_filepath = output_dir / output_filename

        try:
            with open(output_filepath, 'w') as f:
                yaml.dump(boltz_config, f, sort_keys=False, default_flow_style=False, indent=2)
        except yaml.YAMLError as e:
            raise IOError(f"Error writing Boltz YAML file '{output_filepath}': {e}")
        except Exception as e:
            raise IOError(f"An unexpected error occurred while writing Boltz YAML file '{output_filepath}': {e}")

        return output_filepath


# Example Usage (commented out)
# if __name__ == '__main__':
#     from input_handler import InputHandler
#     test_handler = InputHandler()
#     # Test with FASTA
#     # ... (fasta test setup) ...
#     # job_fasta = test_handler.parse_input_file("test.fasta")
#     # if job_fasta:
#     #     cg = ConfigGenerator()
#     #     cg.generate_af3_json_from_job_input(job_fasta, "./test_output")

#     # Test with an input AF3 JSON that has MSA paths
#     dummy_af3_input_content = {
#         "name": "complex_with_msa",
#         "modelSeeds": [456],
#         "sequences": [
#             {"protein": {"id": "P1", "sequence": "MPEPTIDE", "unpairedMsaPath": "msas/p1_msa.a3m"}},
#             {"rna": {"id": "R1", "sequence": "ACGUACGU", "unpairedMsaPath": "msas/r1_msa.sto"}} # Example diff format
#         ],
#         "dialect": "alphafold3",
#         "version": 1
#     }
#     input_json_path = "_test_input_af3.json"
#     with open(input_json_path, "w") as f: json.dump(dummy_af3_input_content, f)
    
#     job_from_json = test_handler.parse_input_file(input_json_path)
#     if job_from_json:
#         print(f"Parsed from JSON: {job_from_json.name_stem}, raw_type: {job_from_json.raw_input_type}")
#         print(f"Input MSA paths found: {job_from_json.input_msa_paths}")
#         cg = ConfigGenerator()
#         generated_json_path = cg.generate_af3_json_from_job_input(job_from_json, "./test_output_json_input")
#         if generated_json_path and os.path.exists(generated_json_path):
#             print(f"Generated JSON: {generated_json_path}")
#             with open(generated_json_path, 'r') as f_read: print(f_read.read()) # Verify MSA paths are included
#             # os.remove(generated_json_path)
#             # if not os.listdir("./test_output_json_input"): os.rmdir("./test_output_json_input")
#     os.remove(input_json_path) 