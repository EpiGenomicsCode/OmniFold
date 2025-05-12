import json
import os
import random
import logging # Added logging
from typing import List, Optional, Dict, Any
import yaml
from pathlib import Path

from .util.definitions import JobInput, SequenceInfo, SequenceType
from .af3_models import (
    Af3Input, Protein, ProteinChain, RNA, RNAChain, DNA, DNAChain, Ligand, LigandMolecule,
    MolId, ProtSeq, RNASeq, DNASeq # Type aliases for validation
)

logger = logging.getLogger(__name__) # Setup logger

class ConfigGenerator:
    def __init__(self):
        """Initializes the ConfigGenerator."""
        pass

    def generate_configs(self, job_input: JobInput, output_dir: Path) -> Optional[Dict[str, str]]:
        """
        Generates configuration files for both AlphaFold3 and Boltz-1 based on the job input.

        Args:
            job_input: The JobInput object containing sequences and potentially MSA info.
            output_dir: The base output directory for the entire job.

        Returns:
            A dictionary containing paths to the generated config files (e.g., 
            {"af3_config_path": "...", "boltz_config_path": "..."}), or None if generation fails.
        """
        config_output_dir = output_dir / "configs"
        config_output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Generating model configs in: {config_output_dir}")
        
        config_paths = {}
        try:
            # Generate AF3 config if needed (always generate for now, runner uses it)
            af3_config_path = self._generate_af3_json_from_job_input(job_input, config_output_dir)
            if af3_config_path:
                config_paths["af3_config_path"] = str(af3_config_path)
            else:
                 logger.warning("Failed to generate AlphaFold3 JSON config.") # Don't fail pipeline if only one fails?

            # Generate Boltz config if needed
            # Check if there are sequences suitable for Boltz first?
            if any(s.molecule_type != 'unknown' for s in job_input.sequences):
                boltz_config_path = self._generate_boltz_yaml_from_job_input(job_input, config_output_dir)
                if boltz_config_path:
                    config_paths["boltz_config_path"] = str(boltz_config_path)
                else:
                    logger.warning("Failed to generate Boltz YAML config.")
            else:
                logger.info("No sequences suitable for Boltz found, skipping YAML generation.")

            if not config_paths:
                 logger.error("Failed to generate any model configuration files.")
                 return None
                 
            return config_paths

        except Exception as e:
            logger.error(f"Error during config generation: {e}", exc_info=True)
            return None

    def _generate_af3_json_from_job_input(
        self, 
        job_input: JobInput,
        config_output_dir: Path,
    ) -> Optional[Path]:
        """
        Generates an AlphaFold3 JSON input file from a JobInput object.
        Prioritizes MSA paths from job_input.protein_id_to_a3m_path if available.

        Args:
            job_input: The JobInput object from InputHandler.
            config_output_dir: Directory to save the generated JSON file.

        Returns:
            Path to the generated JSON file or None on failure.
        """
        if not job_input or not job_input.sequences:
            logger.error("Cannot generate AF3 JSON: JobInput is empty or contains no sequences.")
            return None

        json_file_path = config_output_dir / f"{job_input.name_stem}_af3_generated.json"
        protein_msa_paths = job_input.get("protein_id_to_a3m_path", {}) # Use the dict from Orchestrator
        input_msa_paths = job_input.input_msa_paths # Fallback from original input JSON

        try:
            af3_sequences = []
            processed_ids = set()
            for seq_info in job_input.sequences:
                chain_id = seq_info.chain_id
                if chain_id in processed_ids:
                    logger.warning(f"Duplicate chain ID '{chain_id}' found in job_input sequences. Check InputHandler logic. Skipping duplicate.")
                    continue
                
                msa_path_for_chain = None

                # Priority 1: MSA path from Orchestrator (extracted/generated)
                if protein_msa_paths and chain_id in protein_msa_paths:
                    msa_path_for_chain = protein_msa_paths[chain_id]
                    logger.debug(f"Using extracted/generated MSA path for AF3 chain {chain_id}: {msa_path_for_chain}")
                # Priority 2: MSA path from the input AF3 JSON itself (less likely needed here now)
                elif input_msa_paths and chain_id in input_msa_paths:
                    msa_path_for_chain = input_msa_paths[chain_id]
                    logger.debug(f"Using input MSA path for AF3 chain {chain_id}: {msa_path_for_chain}")
                
                common_chain_args = {
                    "id": chain_id,
                    "sequence": seq_info.sequence
                }
                
                entity_to_add = None
                if seq_info.molecule_type == "protein":
                    protein_chain_args = common_chain_args.copy()
                    if msa_path_for_chain: 
                        protein_chain_args["unpairedMsaPath"] = str(Path(msa_path_for_chain).resolve()) # Use absolute path
                        # AF3 requires pairedMsa field to be present if unpairedMsaPath is, even if empty.
                        # Setting to empty string signifies we want to use the unpairedMsa but run pairedMSA-free.
                        # If no msa_path_for_chain, let AF3 pipeline handle it (or run MSA-free if pipeline disabled)
                        # We assume pairedMsa is handled by the AF3 data pipeline if run, or not needed otherwise.
                        protein_chain_args["pairedMsaPath"] = None # Explicitly null if not using our MSA?
                        protein_chain_args["pairedMsa"] = "" # Crucial: Set pairedMsa to empty string if unpairedMsaPath is provided
                    else:
                        # Let AF3 pipeline handle it if enabled, or run MSA free
                        pass 
                    
                    protein_chain = ProteinChain(**protein_chain_args)
                    entity_to_add = Protein(protein=protein_chain)
                elif seq_info.molecule_type == "rna":
                    rna_chain_args = common_chain_args.copy()
                    if msa_path_for_chain: 
                         rna_chain_args["unpairedMsaPath"] = str(Path(msa_path_for_chain).resolve())
                    rna_chain = RNAChain(**rna_chain_args)
                    entity_to_add = RNA(rna=rna_chain)
                elif seq_info.molecule_type == "dna":
                    dna_chain_args = common_chain_args.copy()
                    # Add MSA path if needed for DNA in future?
                    dna_chain = DNAChain(**dna_chain_args)
                    entity_to_add = DNA(dna=dna_chain)
                elif seq_info.molecule_type == "ligand_ccd":
                    ligand_mol = LigandMolecule(id=chain_id, ccdCodes=[seq_info.sequence])
                    entity_to_add = Ligand(ligand=ligand_mol)
                elif seq_info.molecule_type == "ligand_smiles":
                    ligand_mol = LigandMolecule(id=chain_id, smiles=seq_info.sequence)
                    entity_to_add = Ligand(ligand=ligand_mol)
                else: # unknown
                    logger.warning(f"Skipping sequence {seq_info.original_name} (ID: {chain_id}) for AF3 JSON due to unknown type: {seq_info.molecule_type}")
                    continue
                
                af3_sequences.append(entity_to_add)
                processed_ids.add(chain_id)
            
            if not af3_sequences:
                logger.error("No valid sequences could be processed for AF3 JSON.")
                return None

            # Use seeds from input if provided, otherwise generate one
            model_seeds = job_input.get("model_seeds") or [random.randint(1, 100000)]

            af3_input_data = Af3Input(
                name=job_input.name_stem,
                modelSeeds=model_seeds,
                sequences=af3_sequences,
                bondedAtomPairs=job_input.get("bonded_atom_pairs") 
                # userCCD/userCCDPath would need to be handled if needed
            )

            json_string = af3_input_data.model_dump_json(indent=2, by_alias=True, exclude_none=True)
            
            with open(json_file_path, "w") as f:
                f.write(json_string)
            
            logger.info(f"AlphaFold3 JSON config generated at: {json_file_path}")
            return json_file_path

        except ValueError as ve: 
            logger.error(f"Validation error generating AF3 JSON: {ve}")
            return None
        except Exception as e: 
            logger.error(f"Error generating AlphaFold3 JSON for {job_input.name_stem} at {json_file_path}: {e}", exc_info=True)
            # Clean up potentially partially written file?
            # if json_file_path.exists(): json_file_path.unlink()
            return None

    def _generate_boltz_yaml_from_job_input(
        self,
        job_input: JobInput,
        config_output_dir: Path
    ) -> Optional[Path]:
        """
        Generates a Boltz-1 YAML configuration file from a JobInput object.
        Prioritizes MSA paths from job_input["protein_id_to_a3m_path"].
        The YAML will be saved in the specified config_output_dir.
        Returns the path to the generated YAML file or None on failure.
        """
        boltz_config = {"version": 1, "sequences": []}
        sequences_list = []
        protein_msa_paths = job_input.get("protein_id_to_a3m_path", {}) # Use the dict from Orchestrator
        input_msa_paths = job_input.input_msa_paths # Fallback
        processed_ids = set()

        for seq_info in job_input.sequences:
            chain_id = seq_info.chain_id
            if chain_id in processed_ids:
                 logger.warning(f"Duplicate chain ID '{chain_id}' found in job_input sequences for Boltz YAML. Skipping duplicate.")
                 continue

            entity = {}
            msa_path_for_chain = None
            
            # Determine MSA path for this chain
            if protein_msa_paths and chain_id in protein_msa_paths:
                msa_path_for_chain = protein_msa_paths[chain_id]
                logger.debug(f"Using extracted/generated MSA path for Boltz chain {chain_id}: {msa_path_for_chain}")
            elif input_msa_paths and chain_id in input_msa_paths: # Fallback
                msa_path_for_chain = input_msa_paths[chain_id]
                logger.debug(f"Using input MSA path for Boltz chain {chain_id}: {msa_path_for_chain}")

            if seq_info.molecule_type == "protein":
                protein_data = {
                    "id": chain_id,
                    "sequence": seq_info.sequence
                }
                if msa_path_for_chain:
                    protein_data["msa"] = str(Path(msa_path_for_chain).resolve()) # Use absolute path
                else:
                    logger.info(f"No MSA path found for protein {chain_id}. Boltz will run MSA-free or generate its own if configured.")
                entity["protein"] = protein_data
            elif seq_info.molecule_type == "rna":
                entity["rna"] = {"id": chain_id, "sequence": seq_info.sequence}
                # Add MSA path for RNA if needed/supported by Boltz
            elif seq_info.molecule_type == "dna":
                entity["dna"] = {"id": chain_id, "sequence": seq_info.sequence}
            elif seq_info.molecule_type == "ligand_smiles":
                entity["ligand"] = {"id": chain_id, "smiles": seq_info.sequence}
            elif seq_info.molecule_type == "ligand_ccd":
                entity["ligand"] = {"id": chain_id, "ccd": seq_info.sequence}
            elif seq_info.molecule_type == "unknown":
                logger.warning(f"Sequence '{seq_info.original_name}' (chain {chain_id}) has unknown type. Skipping for Boltz YAML.")
                continue
            else:
                logger.warning(f"Unhandled molecule type '{seq_info.molecule_type}' for '{seq_info.original_name}'. Skipping for Boltz YAML.")
                continue
            
            sequences_list.append(entity)
            processed_ids.add(chain_id)

        if not sequences_list:
            logger.error("No suitable sequences found in JobInput to generate Boltz YAML.")
            return None

        boltz_config["sequences"] = sequences_list

        if job_input.get("constraints"):
            boltz_config["constraints"] = job_input["constraints"]
        
        output_filename = f"{job_input.name_stem}_boltz_generated.yaml"
        output_filepath = config_output_dir / output_filename

        try:
            with open(output_filepath, 'w') as f:
                yaml.dump(boltz_config, f, sort_keys=False, default_flow_style=False, indent=2)
            logger.info(f"Boltz YAML config generated at: {output_filepath}")
            return output_filepath
        except yaml.YAMLError as e:
            logger.error(f"Error writing Boltz YAML file '{output_filepath}': {e}")
            return None
        except Exception as e:
            logger.error(f"An unexpected error occurred while writing Boltz YAML file '{output_filepath}': {e}", exc_info=True)
            # if output_filepath.exists(): output_filepath.unlink()
            return None
