import json
import os
import random
import logging
from typing import List, Optional, Dict, Any
import yaml
from pathlib import Path

from .util.definitions import JobInput, SequenceInfo, SequenceType
from .util.msa_utils import is_a3m_singleton
from .af3_models import (
    Af3Input, Protein, ProteinChain, RNA, RNAChain, DNA, DNAChain, Ligand, LigandMolecule,
    MolId, ProtSeq, RNASeq, DNASeq
)

logger = logging.getLogger(__name__)

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
            if job_input.original_af3_config_path and job_input.has_msa:
                config_paths["af3_config_path"] = job_input.original_af3_config_path
                logger.info(f"Using user-provided AF3 JSON (with MSAs) directly for AF3 inference: {config_paths['af3_config_path']}")
            
            elif job_input.af3_data_json:
                logger.info(f"AF3 MSAs are available from pipeline output: {job_input.af3_data_json}")
                logger.info("Creating final AF3 inference JSON by merging MSA data with original/default seeds and parameters.")
                inference_json_name = f"{job_input.name_stem}_af3_inference_final.json"
                final_af3_path = self._create_af3_inference_json_with_merged_data(
                    job_input=job_input,
                    msa_data_source_json_path=job_input.af3_data_json,
                    config_output_dir=config_output_dir,
                    target_filename=inference_json_name
                )
                if final_af3_path:
                    config_paths["af3_config_path"] = str(final_af3_path)
                    af3_config_generated_internally = True
                else:
                    logger.error("Failed to create merged AF3 JSON for inference.")
            
            elif not job_input.is_boltz_config:
                logger.info("Attempting to generate new AlphaFold3 JSON config for inference using A3M paths.")
                af3_inference_filename = f"{job_input.name_stem}_af3_inference_from_a3m.json"

                af3_config_path = self._generate_af3_json_from_job_input(job_input, config_output_dir, af3_inference_filename)
                if af3_config_path:
                    config_paths["af3_config_path"] = str(af3_config_path)
                    af3_config_generated_internally = True
                else:
                    logger.warning("Failed to generate AlphaFold3 JSON config for inference from A3M paths.")
            else:
                logger.info("Skipping AlphaFold3 config generation (input is Boltz YAML or no suitable AF3 data source for inference).")

            if job_input.protein_id_to_a3m_path and any(s.molecule_type == 'protein' for s in job_input.sequences):
                logger.info("MSAs were generated by the pipeline. Generating a new Boltz YAML to include these MSA paths.")
                boltz_inference_filename = f"{job_input.name_stem}_boltz_inference_generated_with_msas.yaml"
                boltz_config_path = self._generate_boltz_yaml_from_job_input(job_input, config_output_dir, boltz_inference_filename, cli_config)
                if boltz_config_path:
                    config_paths["boltz_config_path"] = str(boltz_config_path)
                    boltz_config_generated_internally = True
                else:
                    logger.error("CRITICAL: Failed to generate Boltz YAML with newly created MSA paths. Boltz execution will likely fail.")
            
            elif job_input.original_boltz_config_path:
                config_paths["boltz_config_path"] = job_input.original_boltz_config_path
                logger.info(f"Using user-provided Boltz YAML directly: {config_paths['boltz_config_path']}. MSAs were not (re)generated by this pipeline run for this input, or no protein sequences present.")
            
            elif any(s.molecule_type != 'unknown' for s in job_input.sequences):
                logger.info("Attempting to generate new Boltz YAML config for inference (no prior MSAs from pipeline for these sequences and no original Boltz config).")
                boltz_inference_filename = f"{job_input.name_stem}_boltz_inference_generated.yaml"
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
        job_input: JobInput,
        msa_data_source_json_path: str,
        config_output_dir: Path,
        target_filename: str
    ) -> Optional[Path]:
        """Creates a new AF3 JSON for inference by taking MSAs from msa_data_source_json_path
           and other parameters (seeds, name, bondedAtoms) from job_input."""
        try:
            with open(msa_data_source_json_path, 'r') as f:
                msa_data = json.load(f)
            
            final_sequences_for_pydantic = []
            for entity_spec_msa_source in msa_data.get("sequences", []):
                final_sequences_for_pydantic.append(entity_spec_msa_source)

            if not final_sequences_for_pydantic:
                logger.error(f"No sequence data found in MSA source JSON: {msa_data_source_json_path}")
                return None

            seeds_for_inference = job_input.model_seeds
            if seeds_for_inference is None:
                seeds_for_inference = msa_data.get("modelSeeds")
                if seeds_for_inference is None:
                    seeds_for_inference = [random.randint(1, 100000)]
                    logger.info(f"Using random seed for merged AF3 inference JSON as no other seeds found: {seeds_for_inference}")
                else:
                    logger.info(f"Using seeds from MSA data source for merged AF3 inference JSON: {seeds_for_inference}")
            else:
                logger.info(f"Using seeds from job_input (original/CLI) for merged AF3 inference JSON: {seeds_for_inference}")

            msa_version = msa_data.get("version")
            msa_dialect = msa_data.get("dialect", "alphafold3") # Default to alphafold3 if not present

            if msa_version is None:
                logger.error(f"'version' field is missing in the MSA data source JSON: {msa_data_source_json_path}. Cannot create valid Af3Input.")
                return None

            af3_input_for_inference = Af3Input(
                name=job_input.name_stem,
                modelSeeds=seeds_for_inference,
                sequences=final_sequences_for_pydantic,
                bondedAtomPairs=job_input.bonded_atom_pairs,
                version=msa_version, # Pass the extracted version
                dialect=msa_dialect  # Pass the extracted or default dialect
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
        filename: str
    ) -> Optional[Path]:
        """
        Generates an AlphaFold3 JSON input file from a JobInput object.
        Used for generating AF3 input for MSA pipeline or for inference if no direct _data.json is available.
        """
        if not job_input or not job_input.sequences:
            logger.error("Cannot generate AF3 JSON: JobInput is empty or contains no sequences.")
            return None

        json_file_path = config_output_dir / filename
        protein_msa_paths = job_input.protein_id_to_a3m_path
        input_msa_paths = job_input.input_msa_paths

        try:
            af3_sequences = []
            processed_ids = set()
            for seq_info in job_input.sequences:
                chain_id = seq_info.chain_id
                if chain_id in processed_ids:
                    logger.warning(f"Duplicate chain ID '{chain_id}' found in job_input sequences. Check InputHandler logic. Skipping duplicate.")
                    continue
                
                common_chain_args = {
                    "id": chain_id,
                    "sequence": seq_info.sequence
                }
                
                entity_to_add = None
                if seq_info.molecule_type == "protein":
                    protein_chain_args = common_chain_args.copy()
                    
                    # Correctly check for and use the generated MSA paths
                    if isinstance(protein_msa_paths, dict) and "unpaired" in protein_msa_paths:
                        unpaired_path = protein_msa_paths.get("unpaired", {}).get(chain_id)
                        paired_path = protein_msa_paths.get("paired", {}).get(chain_id)

                        # Paths in the JSON must be relative to the config file's location
                        if unpaired_path:
                            protein_chain_args["unpairedMsaPath"] = os.path.relpath(unpaired_path, config_output_dir)
                        if paired_path:
                            protein_chain_args["pairedMsaPath"] = os.path.relpath(paired_path, config_output_dir)
                        
                        if unpaired_path and not paired_path:
                            protein_chain_args["pairedMsa"] = ""
                    
                    # When providing MSAs, we must also explicitly disable templates if we don't have any
                    protein_chain_args["templates"] = []
                    
                    protein_chain = ProteinChain(**protein_chain_args)
                    entity_to_add = Protein(protein=protein_chain)
                elif seq_info.molecule_type == "rna":
                    # Simplified logic for RNA as an example, assuming only unpaired for now
                    rna_chain_args = common_chain_args.copy()
                    unpaired_path = protein_msa_paths.get("unpaired", {}).get(chain_id)
                    if unpaired_path: 
                         rna_chain_args["unpairedMsaPath"] = str(Path(unpaired_path).resolve())
                    rna_chain = RNAChain(**rna_chain_args)
                    entity_to_add = RNA(rna=rna_chain)
                elif seq_info.molecule_type == "dna":
                    dna_chain_args = common_chain_args.copy()
                    dna_chain = DNAChain(**dna_chain_args)
                    entity_to_add = DNA(dna=dna_chain)
                elif seq_info.molecule_type == "ligand_ccd":
                    ligand_mol = LigandMolecule(id=chain_id, ccdCodes=[seq_info.sequence])
                    entity_to_add = Ligand(ligand=ligand_mol)
                elif seq_info.molecule_type == "ligand_smiles":
                    ligand_mol = LigandMolecule(id=chain_id, smiles=seq_info.sequence)
                    entity_to_add = Ligand(ligand=ligand_mol)
                else:
                    logger.warning(f"Skipping sequence {seq_info.original_name} (ID: {chain_id}) for AF3 JSON due to unknown type: {seq_info.molecule_type}")
                    continue
                
                af3_sequences.append(entity_to_add)
                processed_ids.add(chain_id)
            
            if not af3_sequences:
                logger.error("No valid sequences could be processed for AF3 JSON.")
                return None

            current_model_seeds = job_input.model_seeds
            if current_model_seeds is None:
                 current_model_seeds = [random.randint(1, 100000)]
                 logger.info(f"Generating random seed for AF3 JSON ({filename}): {current_model_seeds}")

            af3_input_data = Af3Input(
                name=job_input.name_stem,
                modelSeeds=current_model_seeds,
                sequences=af3_sequences,
                bondedAtomPairs=job_input.bonded_atom_pairs,
                version=3,  # Default version as integer
                dialect="alphafold3"  # Default dialect
            )

            json_string = af3_input_data.model_dump_json(indent=2, by_alias=True, exclude_none=True)
            
            with open(json_file_path, "w") as f:
                f.write(json_string)
            
            logger.info(f"AlphaFold3 JSON config generated at: {json_file_path}")
            logger.info(f"--- AlphaFold3 Input JSON Content ---\n{json_string}\n------------------------------------")
            return json_file_path

        except ValueError as ve: 
            logger.error(f"Validation error generating AF3 JSON: {ve}")
            return None
        except Exception as e: 
            logger.error(f"Error generating AlphaFold3 JSON for {job_input.name_stem} at {json_file_path}: {e}", exc_info=True)
            return None

    def _generate_boltz_yaml_from_job_input(
        self,
        job_input: JobInput,
        config_output_dir: Path,
        filename: str,
        cli_config: dict
    ) -> Optional[Path]:
        """
        Generates a Boltz-1 YAML configuration file.
        Adjusts number of predictions based on AF3 seeds if CLI default for Boltz samples is used.
        """
        yaml_file_path = config_output_dir / filename
        boltz_config = {
            "name": job_input.name_stem,
            "version": 1, 
            "sequences": [] 
        }

        boltz_config["sampling"] = {
            "n_preds": cli_config.get("boltz_diffusion_samples", 1),
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

        # This logic to set n_preds based on AF3 seeds should be more robust
        # It assumes the default was used for Boltz based on a flag from the CLI parser.
        if cli_config.get("boltz_diffusion_samples_is_default", False) and \
           job_input.num_model_seeds_from_input is not None and \
           job_input.num_model_seeds_from_input > 0:
            boltz_config["sampling"]["n_preds"] = job_input.num_model_seeds_from_input
            logger.info(f"Boltz 'n_preds' set to {job_input.num_model_seeds_from_input} based on AF3 input modelSeeds, as Boltz CLI default was used.")

        if job_input.constraints:
            boltz_config["constraints"] = job_input.constraints
            logger.info(f"Added constraints to Boltz YAML from job input.")

        segments_list = []
        protein_msa_paths = job_input.protein_id_to_a3m_path # This will now be a dict of dicts
        input_msa_paths = job_input.input_msa_paths
        processed_ids = set()
        sequence_to_msa_path_map: Dict[str, str] = {}

        # --- New: Prioritize Boltz CSV MSAs if they were generated ---
        if job_input.boltz_csv_msa_dir and Path(job_input.boltz_csv_msa_dir).is_dir():
            logger.info(f"Using pre-generated Boltz CSV MSAs from: {job_input.boltz_csv_msa_dir}")
            job_root_path = config_output_dir.parent
            for seq_info in job_input.sequences:
                if seq_info.molecule_type == "protein":
                    # If we've already stored an MSA path for this exact sequence, no need to recompute –
                    # duplicate chains with identical sequences MUST use the same MSA file (Boltz requirement).
                    if seq_info.sequence in sequence_to_msa_path_map:
                        continue
                    host_csv_path = Path(job_input.boltz_csv_msa_dir) / f"{seq_info.chain_id}.csv"
                    if host_csv_path.is_file() and host_csv_path.stat().st_size > 10: # Check for non-trivial file
                        relative_csv_dir = Path(job_input.boltz_csv_msa_dir).relative_to(job_root_path)
                        container_csv_path = str(Path("/data/job_output") / relative_csv_dir / f"{seq_info.chain_id}.csv")
                        sequence_to_msa_path_map[seq_info.sequence] = container_csv_path
                        logger.info(f"Mapped sequence for chain {seq_info.chain_id} to CSV MSA: {container_csv_path}")
                    else:
                        logger.warning(
                            f"CSV file for chain {seq_info.chain_id} not found or is empty at {host_csv_path}. "
                            "Will default to 'empty' MSA for this chain and all identical sequences."
                        )
                        # Only set to empty if this sequence hasn't been mapped before.
                        sequence_to_msa_path_map.setdefault(seq_info.sequence, "empty")
        
        for seq_info in job_input.sequences:
            chain_id = seq_info.chain_id
            if chain_id in processed_ids:
                 logger.warning(f"Duplicate chain ID '{chain_id}' found in job_input sequences for Boltz YAML. Skipping duplicate.")
                 continue

            entity = {}
            
            if seq_info.molecule_type == "protein":
                # Retrieve the shared MSA path for this sequence (set above). If not present, default to empty.
                msa_path = sequence_to_msa_path_map.get(seq_info.sequence, "empty")
                if msa_path == "empty":
                    logger.info(f"Setting MSA to 'empty' for Boltz protein {chain_id} (no valid CSV found).")

                entity = {
                    "protein": {
                        "id": chain_id,
                        "sequence": seq_info.sequence,
                        "msa": msa_path
                    }
                }
            elif seq_info.molecule_type == "rna":
                entity = {"id": chain_id, "sequence": seq_info.sequence, "type": "rna"} 
            elif seq_info.molecule_type == "dna":
                entity = {"id": chain_id, "sequence": seq_info.sequence, "type": "dna"}
            elif seq_info.molecule_type == "ligand_smiles":
                entity = { "ligand": { "id": chain_id, "smiles": seq_info.sequence } }
            elif seq_info.molecule_type == "ligand_ccd":
                entity = { "ligand": { "id": chain_id, "ccd": seq_info.sequence } }
            else:
                logger.warning(f"Unhandled molecule type '{seq_info.molecule_type}' for '{seq_info.original_name}'. Skipping for Boltz YAML.")
                continue

            if entity:
                segments_list.append(entity)
            processed_ids.add(chain_id)

        if not segments_list:
            logger.warning(f"No suitable segments generated for Boltz YAML for job {job_input.name_stem}")

        boltz_config["sequences"] = segments_list 
        
        try:
            with open(yaml_file_path, 'w') as f:
                yaml.dump(boltz_config, f, sort_keys=False, Dumper=yaml.SafeDumper)
            
            # Log the content of the generated YAML
            with open(yaml_file_path, 'r') as f:
                yaml_content = f.read()
            logger.info(f"Boltz YAML config generated at: {yaml_file_path}")
            logger.info(f"--- Boltz Input YAML Content ---\n{yaml_content.strip()}\n---------------------------------")
            return yaml_file_path
        except Exception as e:
            logger.error(f"Error writing Boltz YAML to {yaml_file_path}: {e}", exc_info=True)
            return None
