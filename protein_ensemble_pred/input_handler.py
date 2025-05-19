import os
import json 
import pyfastx
from typing import List, Optional, Dict, Union, Any
from pathlib import Path
import yaml
import logging

from .util.definitions import JobInput, SequenceInfo, idgen, as_entity, SequenceType
from .af3_models import Af3Input as Af3PydanticModel  
from .af3_models import Protein, RNA, DNA, Ligand  

logger = logging.getLogger(__name__)

class InputHandler:
    def __init__(self):
        """Initializes the InputHandler."""
        pass

    def parse_input(self, input_filepath: str) -> Optional[JobInput]:
        """
        Parses the input file and returns a JobInput object.

        Args:
            input_filepath: Path to the input file.

        Returns:
            A JobInput object if parsing is successful, None otherwise.
        """
        input_path = Path(input_filepath)
        if not input_path.is_file():
            logger.error(f"Input file not found: {input_filepath}")
            return None

        file_extension = input_path.suffix.lower()
        name_stem = input_path.stem

        try:
            if file_extension in [".fasta", ".fa", ".fna", ".faa"]:
                logger.info(f"Parsing FASTA file: {input_filepath}")
                return self._parse_fasta(input_path, name_stem)
            elif file_extension == ".json":
                logger.info(f"Attempting to parse as AlphaFold3 JSON: {input_filepath}")
                return self._parse_af3_json(input_path, name_stem)
            elif file_extension in [".yaml", ".yml"]:
                logger.info(f"Attempting to parse as Boltz YAML: {input_filepath}")
                return self._parse_boltz_yaml(input_path, name_stem)
            else:
                logger.error(f"Unsupported file extension: {file_extension}. Please use FASTA, AF3 JSON, or Boltz YAML.")
                return None
        except Exception as e:
            logger.error(f"Failed to parse input file {input_filepath}: {e}", exc_info=True)
            return None

    def _parse_fasta(self, file_path: Path, name_stem: str) -> Optional[JobInput]:
        """
        Parses a FASTA file.
        """
        sequences_info: List[SequenceInfo] = []
        chain_id_generator = idgen()
        current_sequence_parts: List[str] = []
        current_original_name: Optional[str] = None

        try:
            with open(file_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    if line.startswith('>'):
                        if current_original_name is not None and current_sequence_parts:
                            full_sequence = "".join(current_sequence_parts)
                            chain_id = next(chain_id_generator)
                            sequences_info.append(as_entity(full_sequence, chain_id, current_original_name))
                        

                        current_original_name = line[1:].strip()  
                        current_sequence_parts = []
                    elif current_original_name is not None: 
                        cleaned_line = ''.join(filter(str.isalpha, line))
                        current_sequence_parts.append(cleaned_line)
            
            if current_original_name is not None and current_sequence_parts:
                full_sequence = "".join(current_sequence_parts)
                chain_id = next(chain_id_generator)
                sequences_info.append(as_entity(full_sequence, chain_id, current_original_name))

            if not sequences_info:
                logger.error(f"No valid sequences found in FASTA file: {file_path}")
                return None

            return JobInput(
                name_stem=name_stem,
                sequences=sequences_info,
                raw_input_type="fasta",
                output_dir="" 
            )

        except StopIteration:
            logger.error(f"Exceeded maximum number of chain IDs (ZZ). Too many sequences in FASTA: {file_path}")
            return None
        except Exception as e:
            logger.error(f"Error parsing FASTA file {file_path}: {e}", exc_info=True)
            return None

    def _parse_af3_json(self, file_path: Path, name_stem: str) -> Optional[JobInput]:
        """
        Parses an AlphaFold3 JSON file.
        Extracts sequences, IDs, types, MSA paths, seeds, and bonds.
        Performs a content-based check for MSAs.
        """
        sequences_info: List[SequenceInfo] = []
        model_seeds_from_json: Optional[List[int]] = None
        num_model_seeds_val: Optional[int] = None
        bonded_atom_pairs_val: Optional[List[Any]] = None
        has_msa_content = False

        try:
            with open(file_path, 'r') as f:
                data = json.load(f)
            
            try:
                af3_input_model = Af3PydanticModel(**data)
                data = af3_input_model.model_dump(by_alias=True) 
            except Exception as pydantic_error:
                logger.warning(f"Input JSON {file_path} failed Pydantic validation for Af3Input: {pydantic_error}. Attempting to parse leniently.")

            job_name_from_json = data.get("name", name_stem)
            
            if "modelSeeds" in data and isinstance(data["modelSeeds"], list) and data["modelSeeds"]:
                model_seeds_from_json = [int(s) for s in data["modelSeeds"] if isinstance(s, (int, str)) and str(s).isdigit()]
                if model_seeds_from_json: 
                    num_model_seeds_val = len(model_seeds_from_json)
                else: 
                    model_seeds_from_json = None 
                    num_model_seeds_val = None
            
            bonded_atom_pairs_val = data.get("bondedAtomPairs")

            chain_id_generator = idgen()
            processed_ids = set()

            for entity_spec in data.get("sequences", []):
                if not isinstance(entity_spec, dict) or len(entity_spec) != 1:
                    logger.warning(f"Skipping invalid entity spec in AF3 JSON: {entity_spec}")
                    continue
                
                entity_type_key, entity_data = list(entity_spec.items())[0]
                if not isinstance(entity_data, dict):
                    logger.warning(f"Skipping invalid entity data for type '{entity_type_key}': {entity_data}")
                    continue

                raw_ids = entity_data.get("id")
                chain_ids: List[str]
                if isinstance(raw_ids, str): chain_ids = [raw_ids]
                elif isinstance(raw_ids, list): chain_ids = [str(i) for i in raw_ids if isinstance(i, str)]
                else:
                    logger.warning(f"Entity type '{entity_type_key}' missing or invalid 'id'. Assigning automatically.")
                    try: chain_ids = [next(chain_id_generator)]; logger.warning(f"Assigned automatic ID: {chain_ids[0]}")
                    except StopIteration: logger.error(f"Ran out of chain IDs: {file_path}"); return None
                
                seq: Optional[str] = entity_data.get("sequence")
                seq_type_val: SequenceType
                original_name_val = f"{entity_type_key}_{'_'.join(chain_ids)}"

                entity_has_msa = False
                if entity_type_key == "protein": seq_type_val = "protein"
                elif entity_type_key == "rna": seq_type_val = "rna"
                elif entity_type_key == "dna": seq_type_val = "dna"
                elif entity_type_key == "ligand":
                    if "ccdCodes" in entity_data and entity_data["ccdCodes"]: seq = entity_data["ccdCodes"][0]; seq_type_val = "ligand_ccd"
                    elif "smiles" in entity_data: seq = entity_data["smiles"]; seq_type_val = "ligand_smiles"
                    else: logger.warning(f"Ligand entity invalid: {entity_data}. Skipping."); continue
                else: logger.warning(f"Unsupported entity type '{entity_type_key}'. Skipping."); continue

                if seq is None and entity_type_key != "ligand":
                    logger.warning(f"Entity '{entity_type_key}' ID(s) '{chain_ids}' missing sequence. Skipping.")
                    continue

                if seq_type_val in ["protein", "rna", "dna"]:
                    if entity_data.get("unpairedMsaPath") and Path(file_path.parent / entity_data.get("unpairedMsaPath")).is_file(): # Check if path is valid
                        entity_has_msa = True
                    if not entity_has_msa and entity_data.get("pairedMsaPath") and Path(file_path.parent / entity_data.get("pairedMsaPath")).is_file(): # Check if path is valid
                        entity_has_msa = True
                    
                    unpaired_msa_val = entity_data.get("unpairedMsa")
                    if not entity_has_msa and isinstance(unpaired_msa_val, str) and unpaired_msa_val.strip(): # Non-empty string
                        entity_has_msa = True

                    paired_msa_val = entity_data.get("pairedMsa")
                    if not entity_has_msa and isinstance(paired_msa_val, str) and paired_msa_val.strip(): # Non-empty string
                        entity_has_msa = True

                if entity_has_msa:
                    has_msa_content = True 

                for chain_id in chain_ids:
                    if chain_id in processed_ids: logger.warning(f"Duplicate chain ID '{chain_id}'. Skipping."); continue
                    sequences_info.append(SequenceInfo(original_name=original_name_val, sequence=str(seq), molecule_type=seq_type_val, chain_id=chain_id))
                    processed_ids.add(chain_id)
            
            if not sequences_info: logger.error(f"No valid sequences in AF3 JSON: {file_path}"); return None

            job_input = JobInput(
                name_stem=job_name_from_json,
                sequences=sequences_info,
                raw_input_type="af3_json",
                original_af3_config_path=str(file_path.resolve()), 
                model_seeds=model_seeds_from_json,
                num_model_seeds_from_input=num_model_seeds_val,
                bonded_atom_pairs=bonded_atom_pairs_val,
                has_msa=has_msa_content, 
                output_dir="" 
            )

            if has_msa_content:
                job_input.af3_data_json = str(file_path.resolve()) 
                logger.info(f"Detected MSA content within input AF3 JSON: {file_path}. This file will be used as af3_data_json.")
            else:
                logger.info(f"No MSA content detected within input AF3 JSON: {file_path}. MSA generation may be required.")

            return job_input

        except json.JSONDecodeError as e:
            logger.error(f"Invalid JSON format in {file_path}: {e}", exc_info=True)
            return None
        except StopIteration: 
             logger.error(f"Exceeded maximum number of chain IDs while parsing AF3 JSON: {file_path}")
             return None
        except Exception as e:
            logger.error(f"Error parsing AF3 JSON file {file_path}: {e}", exc_info=True)
            return None

    def _parse_boltz_yaml(self, file_path: Path, name_stem: str) -> Optional[JobInput]:
        sequences_info: List[SequenceInfo] = []
        has_msa_content = False
        constraints_val: Optional[List[Dict]] = None


        try:
            with open(file_path, 'r') as f:
                data = yaml.safe_load(f)

            job_name_from_yaml = data.get("name", name_stem)
            
            if "constraints" in data and data["constraints"]:
                constraints_val = data["constraints"]

            yaml_sequences = data.get("segments", data.get("sequences", []))
            if not isinstance(yaml_sequences, list):
                logger.error(f"Boltz YAML {file_path} 'segments' or 'sequences' key is not a list or is missing.")
                return None

            chain_id_gen = idgen()

            for i, seg_data in enumerate(yaml_sequences):
                if not isinstance(seg_data, dict):
                    logger.warning(f"Skipping non-dictionary segment in Boltz YAML: {seg_data}"); continue

                entity_type = next(iter(seg_data.keys()))
                entity_data = seg_data[entity_type]

                if not isinstance(entity_data, dict):
                    logger.warning(f"Invalid entity data in Boltz YAML: {entity_data}"); continue

                if entity_type == "protein":
                    seq_str = entity_data.get("sequence")
                    molecule_type_val: SequenceType = "protein" # Default for protein
                elif entity_type == "ligand":
                    if "smiles" in entity_data:
                        seq_str = entity_data.get("smiles")
                        molecule_type_val = "ligand_smiles"
                    elif "ccd" in entity_data:
                        seq_str = entity_data.get("ccd")
                        molecule_type_val = "ligand_ccd"
                    else:
                        logger.warning(f"Ligand entity in Boltz YAML is missing both 'smiles' and 'ccd': {entity_data}. Skipping.")
                        continue
                else:
                    logger.warning(f"Unsupported entity type in Boltz YAML: {entity_type}"); continue

                chain_id = entity_data.get("id")
                if not chain_id:
                    try:
                        chain_id = next(chain_id_gen)
                        logger.info(f"Generated chain ID {chain_id} for Boltz sequence index {i}")
                    except StopIteration:
                        logger.error(f"Ran out of chain IDs for Boltz YAML: {file_path}")
                        return None
                
                original_name = chain_id 

                if not seq_str:
                    logger.warning(f"Segment '{chain_id}' in Boltz YAML is missing 'sequence', 'smiles', or 'ccd'. Skipping.")
                    continue
                
                if entity_type == "ligand":
                    seq_info_obj = SequenceInfo(
                        original_name=original_name,
                        sequence=str(seq_str),
                        molecule_type=molecule_type_val, 
                        chain_id=str(chain_id)
                    )
                elif entity_type == "protein": 
                     seq_info_obj = SequenceInfo(
                        original_name=original_name,
                        sequence=str(seq_str),
                        molecule_type="protein", 
                        chain_id=str(chain_id)
                    )
                else: 
                    logger.warning(f"Unexpected entity type '{entity_type}' fell through to as_entity. Check logic.")
                    seq_info_obj = as_entity(str(seq_str), str(chain_id), original_name)

                sequences_info.append(seq_info_obj)

                msa_entry = entity_data.get("msa")
                if msa_entry: 
                    has_msa_content = True 

            if not sequences_info:
                logger.error(f"No valid sequences extracted from Boltz YAML: {file_path}")
                return None

            return JobInput(
                name_stem=job_name_from_yaml,
                sequences=sequences_info,
                raw_input_type="boltz_yaml",
                original_boltz_config_path=str(file_path.resolve()),
                has_msa=has_msa_content,
                is_boltz_config=True, 
                constraints=constraints_val,
                output_dir="" 
            )

        except yaml.YAMLError as e:
            logger.error(f"Invalid YAML format in {file_path}: {e}", exc_info=True)
            return None
        except StopIteration: 
             logger.error(f"Exceeded maximum number of chain IDs while parsing Boltz YAML: {file_path}")
             return None
        except Exception as e:
            logger.error(f"Error parsing Boltz YAML file {file_path}: {e}", exc_info=True)
            return None
