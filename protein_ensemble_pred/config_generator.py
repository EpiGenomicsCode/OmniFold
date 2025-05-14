import json
import os
import random
import logging # Added logging
from typing import List, Optional, Dict, Any
import yaml
from pathlib import Path

from .util.definitions import JobInput, SequenceInfo, SequenceType
from .util.msa_utils import is_a3m_singleton # Import the new helper
from .af3_models import (
    Af3Input, Protein, ProteinChain, RNA, RNAChain, DNA, DNAChain, Ligand, LigandMolecule,
    MolId, ProtSeq, RNASeq, DNASeq # Type aliases for validation
)

logger = logging.getLogger(__name__) # Setup logger

class ConfigGenerator:
    def __init__(self):
        """Initializes the ConfigGenerator."""
        pass

    def generate_configs(self, job_input: JobInput, output_dir: Path, cli_config: dict) -> Optional[Dict[str, str]]:
        """
        Generates configuration files for both AlphaFold3 and Boltz-1.
        Respects original user-provided configs if they exist and are complete.
        """
        config_output_dir = output_dir / "configs"
        config_output_dir.mkdir(parents=True, exist_ok=True)
        logger.info(f"Generating model configs in: {config_output_dir}")
        
        config_paths = {}
        af3_config_generated_internally = False
        boltz_config_generated_internally = False

        try:
            # --- AlphaFold3 Configuration ---
            if job_input.original_af3_config_path and job_input.has_msa:
                config_paths["af3_config_path"] = job_input.original_af3_config_path
                logger.info(f"Using user-provided AF3 JSON (with MSAs) directly for AF3 inference: {config_paths['af3_config_path']}")
            
            elif job_input.af3_data_json: # MSAs were generated (from any input type), result is in af3_data_json
                logger.info(f"AF3 MSAs are available from pipeline output: {job_input.af3_data_json}")
                logger.info("Creating final AF3 inference JSON by merging MSA data with original/default seeds and parameters.")
                inference_json_name = f"{job_input.name_stem}_af3_inference_final.json"
                final_af3_path = self._create_af3_inference_json_with_merged_data(
                    job_input=job_input, # This has the original/default seeds, bonded_atom_pairs etc.
                    msa_data_source_json_path=job_input.af3_data_json, # This has the MSAs from the data pipeline
                    config_output_dir=config_output_dir,
                    target_filename=inference_json_name
                )
                if final_af3_path:
                    config_paths["af3_config_path"] = str(final_af3_path)
                    af3_config_generated_internally = True # We still consider this internal generation
                else:
                    logger.error("Failed to create merged AF3 JSON for inference.")
            
            elif not job_input.is_boltz_config: # No original complete AF3, no AF3 MSA pipeline output, but not Boltz input
                                              # This implies FASTA/other input + ColabFold MSAs, so we need to generate AF3 JSON from A3Ms
                logger.info("Attempting to generate new AlphaFold3 JSON config for inference using A3M paths.")
                af3_inference_filename = f"{job_input.name_stem}_af3_inference_from_a3m.json"
                # _generate_af3_json_from_job_input uses job_input.model_seeds (which has CLI default if applicable)
                # and job_input.protein_id_to_a3m_path
                af3_config_path = self._generate_af3_json_from_job_input(job_input, config_output_dir, af3_inference_filename)
                if af3_config_path:
                    config_paths["af3_config_path"] = str(af3_config_path)
                    af3_config_generated_internally = True
                else:
                    logger.warning("Failed to generate AlphaFold3 JSON config for inference from A3M paths.")
            else:
                logger.info("Skipping AlphaFold3 config generation (input is Boltz YAML or no suitable AF3 data source for inference).")

            # --- Boltz-1 Configuration ---
            if job_input.original_boltz_config_path:
                # User provided Boltz YAML. Use it directly.
                # Per user: "if the configuration of a particular model is given we shouldn't be messing with that input file"
                # So, we don't adjust n_preds here based on AF3 seeds if original Boltz config is used.
                config_paths["boltz_config_path"] = job_input.original_boltz_config_path
                logger.info(f"Using user-provided Boltz YAML directly for Boltz inference: {config_paths['boltz_config_path']}")
            elif any(s.molecule_type != 'unknown' for s in job_input.sequences): # Check if there are any processable sequences for Boltz
                logger.info("Attempting to generate new Boltz YAML config for inference.")
                boltz_inference_filename = f"{job_input.name_stem}_boltz_inference_generated.yaml"
                # Pass cli_config here for n_preds logic
                boltz_config_path = self._generate_boltz_yaml_from_job_input(job_input, config_output_dir, boltz_inference_filename, cli_config)
                if boltz_config_path:
                    config_paths["boltz_config_path"] = str(boltz_config_path)
                    boltz_config_generated_internally = True
                else:
                    logger.warning("Failed to generate Boltz YAML config for inference.")
            else:
                logger.info("No sequences suitable for Boltz found, or original Boltz config provided. Skipping new Boltz YAML generation.")
                
            if not config_paths:
                 logger.error("Failed to identify or generate any model configuration files.")
                 return None
                 
            logger.info(f"Final model config paths: {config_paths}")
            return config_paths

        except Exception as e:
            logger.error(f"Error during config generation: {e}", exc_info=True)
            return None

    def _create_af3_inference_json_with_merged_data(
        self,
        job_input: JobInput, # Contains target seeds, bonded pairs, etc.
        msa_data_source_json_path: str, # Path to *_data.json from MSA pipeline (has MSAs)
        config_output_dir: Path,
        target_filename: str
    ) -> Optional[Path]:
        """Creates a new AF3 JSON for inference by taking MSAs from msa_data_source_json_path
           and other parameters (seeds, name, bondedAtoms) from job_input."""
        try:
            with open(msa_data_source_json_path, 'r') as f:
                msa_data = json.load(f)
            
            # Extract sequence entities (which include MSA fields like unpairedMsa) from the MSA data source
            # We need to ensure these are in the Pydantic model structure if we pass them to Af3Input
            # For simplicity, we rebuild the sequence part using the structure of msa_data,
            # assuming it matches what Af3Input expects for its sequences list.
            final_sequences_for_pydantic = []
            for entity_spec_msa_source in msa_data.get("sequences", []):
                # We assume entity_spec_msa_source is already in the correct format 
                # { "protein": { "id": ..., "sequence": ..., "unpairedMsa": ... } } or similar for RNA/DNA
                # No easy way to map job_input.sequences back to these if chain IDs changed or something.
                # Safest is to trust the structure from msa_data_source, which AF3 produced.
                # The main thing is we will override top-level fields like modelSeeds and name.
                final_sequences_for_pydantic.append(entity_spec_msa_source)

            if not final_sequences_for_pydantic:
                logger.error(f"No sequence data found in MSA source JSON: {msa_data_source_json_path}")
                return None

            # Determine seeds: Priority to job_input.model_seeds (from original input or CLI default)
            seeds_for_inference = job_input.model_seeds
            if seeds_for_inference is None:
                # Fallback to seeds from the MSA data source JSON itself, then random
                seeds_for_inference = msa_data.get("modelSeeds")
                if seeds_for_inference is None:
                    seeds_for_inference = [random.randint(1, 100000)]
                    logger.info(f"Using random seed for merged AF3 inference JSON as no other seeds found: {seeds_for_inference}")
                else:
                    logger.info(f"Using seeds from MSA data source for merged AF3 inference JSON: {seeds_for_inference}")
            else:
                logger.info(f"Using seeds from job_input (original/CLI) for merged AF3 inference JSON: {seeds_for_inference}")

            af3_input_for_inference = Af3Input(
                name=job_input.name_stem, # Use original/intended name_stem
                modelSeeds=seeds_for_inference,
                sequences=final_sequences_for_pydantic, # These have MSAs embedded
                bondedAtomPairs=job_input.bonded_atom_pairs # From original job_input
                # dialect, version, etc., could also be set if needed.
            )

            output_json_path = config_output_dir / target_filename
            with open(output_json_path, "w") as f:
                f.write(af3_input_for_inference.model_dump_json(indent=2, by_alias=True, exclude_none=True))
            logger.info(f"Generated merged AF3 JSON for inference: {output_json_path}")
            return output_json_path

        except Exception as e:
            logger.error(f"Failed to create merged AF3 inference JSON from {msa_data_source_json_path}: {e}", exc_info=True)
            return None

    def _generate_af3_json_from_job_input(
        self, 
        job_input: JobInput,
        config_output_dir: Path,
        filename: str # Explicit filename parameter
    ) -> Optional[Path]:
        """
        Generates an AlphaFold3 JSON input file from a JobInput object.
        Used for generating AF3 input for MSA pipeline or for inference if no direct _data.json is available.
        """
        if not job_input or not job_input.sequences:
            logger.error("Cannot generate AF3 JSON: JobInput is empty or contains no sequences.")
            return None

        json_file_path = config_output_dir / filename # Use provided filename
        protein_msa_paths = job_input.protein_id_to_a3m_path # Use the dict from Orchestrator
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
                        protein_chain_args["unpairedMsaPath"] = str(Path(msa_path_for_chain).resolve()) 
                        protein_chain_args["pairedMsa"] = f">query\n{seq_info.sequence}\n"
                    else:
                        # If no external MSA is provided, do not include MSA-related fields.
                        # AF3 will run MSA-free or use its internal pipeline if those flags are set (not our case for inference).
                        # For pure inference with no precomputed MSA and no pipeline, it would likely also fail.
                        # However, our logic ensures msa_path_for_chain IS populated if MSAs were generated/extracted.
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

            # Use seeds from job_input (captured from original AF3 JSON or set by CLI default via Orchestrator if applicable)
            # If job_input.model_seeds is None (e.g. FASTA input), generate a random seed.
            # This ensures AF3 always gets seeds. Orchestrator/CLI handles default seed logic.
            current_model_seeds = job_input.model_seeds
            if current_model_seeds is None:
                 # Fallback to a single random seed if no seeds are specified anywhere for this job_input.
                 # For the temp AF3 JSON for MSA, job_input.model_seeds is explicitly None.
                 # A random seed here is fine as the MSA pipeline is deterministic for fixed inputs.
                 current_model_seeds = [random.randint(1, 100000)]
                 logger.info(f"Generating random seed for AF3 JSON ({filename}): {current_model_seeds}")

            af3_input_data = Af3Input(
                name=job_input.name_stem,
                modelSeeds=current_model_seeds, # Use seeds from job_input
                sequences=af3_sequences,
                bondedAtomPairs=job_input.bonded_atom_pairs 
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
        config_output_dir: Path,
        filename: str, # Explicit filename parameter
        cli_config: dict # Added cli_config for n_preds logic
    ) -> Optional[Path]:
        """
        Generates a Boltz-1 YAML configuration file.
        Adjusts number of predictions based on AF3 seeds if CLI default for Boltz samples is used.
        """
        yaml_file_path = config_output_dir / filename # Use provided filename
        boltz_config = {
            "name": job_input.name_stem, # Add top-level name
            "version": 1, 
            "segments": [] # Changed from "sequences" to "segments" for Boltz
        }
        # Default sampling settings, can be overridden
        boltz_config["sampling"] = {
            "n_preds": cli_config.get("boltz_diffusion_samples", 1), # Get from cli_config
            "recycling_steps": cli_config.get("boltz_recycling_steps", 3),
            "sampling_steps": cli_config.get("boltz_sampling_steps", 200),
            "step_scale": cli_config.get("boltz_step_scale", 1.638)
        }
        boltz_config["potentials"] = {"use_potentials": not cli_config.get("boltz_no_potentials", False)}
        boltz_config["output"] = {
            "format": cli_config.get("boltz_output_format", "mmcif"),
            "write_full_pae": cli_config.get("boltz_write_full_pae", False),
            "write_full_pde": cli_config.get("boltz_write_full_pde", False)
        }


        # Adjust n_preds based on AF3 input seeds if Boltz n_preds is default
        if cli_config.get("boltz_diffusion_samples_is_default", False) and \
           job_input.num_model_seeds_from_input is not None and \
           job_input.num_model_seeds_from_input > 0:
            boltz_config["sampling"]["n_preds"] = job_input.num_model_seeds_from_input
            logger.info(f"Boltz 'n_preds' set to {job_input.num_model_seeds_from_input} based on AF3 input modelSeeds, as Boltz CLI default was used.")

        if job_input.constraints:
            boltz_config["constraints"] = job_input.constraints
            logger.info(f"Added constraints to Boltz YAML from job input.")

        segments_list = []
        protein_msa_paths = job_input.protein_id_to_a3m_path # Use the dict from Orchestrator
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
                    # "id": chain_id, # Boltz uses "name" within segment for identifier
                    "name": chain_id, 
                    "sequence": seq_info.sequence
                }
                if msa_path_for_chain:
                    if is_a3m_singleton(msa_path_for_chain, seq_info.sequence):
                        protein_data["msa"] = "empty"
                        logger.info(f"A3M for protein {chain_id} is a singleton. Setting Boltz msa to 'empty'.")
                    else:
                        protein_data["msa"] = str(Path(msa_path_for_chain).resolve()) # Use absolute path
                else:
                    logger.info(f"No MSA path found for protein {chain_id}. Setting Boltz msa to 'empty'.")
                    protein_data["msa"] = "empty" # Default to empty if no MSA path
                entity["protein"] = protein_data
            elif seq_info.molecule_type == "rna":
                # entity["rna"] = {"id": chain_id, "sequence": seq_info.sequence}
                # Boltz YAML structure for RNA/DNA needs confirmation, assume similar to protein for now if supported
                segments_list.append({"name": chain_id, "sequence": seq_info.sequence, "type": "rna"})
            elif seq_info.molecule_type == "dna":
                # entity["dna"] = {"id": chain_id, "sequence": seq_info.sequence}
                segments_list.append({"name": chain_id, "sequence": seq_info.sequence, "type": "dna"})
            elif seq_info.molecule_type == "ligand_smiles":
                # entity["ligand"] = {"id": chain_id, "smiles": seq_info.sequence}
                 segments_list.append({"name": chain_id, "smiles": seq_info.sequence, "type": "ligand"})
            elif seq_info.molecule_type == "ligand_ccd":
                # entity["ligand"] = {"id": chain_id, "ccd": seq_info.sequence}
                segments_list.append({"name": chain_id, "ccd": seq_info.sequence, "type": "ligand"}) # Boltz might prefer SMILES
            elif seq_info.molecule_type == "unknown":
                logger.warning(f"Sequence '{seq_info.original_name}' (chain {chain_id}) has unknown type. Skipping for Boltz YAML.")
                continue
            else:
                logger.warning(f"Unhandled molecule type '{seq_info.molecule_type}' for '{seq_info.original_name}'. Skipping for Boltz YAML.")
                continue
            
            segments_list.append(entity)
            processed_ids.add(chain_id)

        if not segments_list:
            logger.warning(f"No suitable segments generated for Boltz YAML for job {job_input.name_stem}")
            # Still might create a file if only general settings were important, but likely an issue.
            # For now, let it proceed to write an empty segments list if that's the case.

        boltz_config["segments"] = segments_list
        
        try:
            with open(yaml_file_path, 'w') as f:
                yaml.dump(boltz_config, f, sort_keys=False, Dumper=yaml.SafeDumper) # Use SafeDumper
            logger.info(f"Boltz YAML config generated at: {yaml_file_path}")
            return yaml_file_path
        except Exception as e:
            logger.error(f"Error writing Boltz YAML to {yaml_file_path}: {e}", exc_info=True)
            return None
