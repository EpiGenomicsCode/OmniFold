import os
import json # For loading JSON file content
import pyfastx
from typing import List, Optional, Dict, Union
from pathlib import Path
import yaml # Added import
import logging

from .util.definitions import JobInput, SequenceInfo, idgen, as_entity, SequenceType
from .af3_models import Af3Input as Af3PydanticModel # Rename to avoid confusion with function
from .af3_models import Protein, RNA, DNA, Ligand # To identify sequence item types

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
                 # Add placeholder or call specific parser if Boltz YAML is distinct
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
                        # Process previous sequence if exists
                        if current_original_name is not None and current_sequence_parts:
                            full_sequence = "".join(current_sequence_parts)
                            chain_id = next(chain_id_generator)
                            sequences_info.append(as_entity(full_sequence, chain_id, current_original_name))
                        
                        # Start new sequence
                        current_original_name = line[1:].strip() # Name is everything after '>'
                        current_sequence_parts = []
                    elif current_original_name is not None: # Only append if we are inside a sequence block
                        # Filter out non-alphabet characters (safer than just upper)
                        cleaned_line = ''.join(filter(str.isalpha, line))
                        current_sequence_parts.append(cleaned_line)
            
            # Process the last sequence in the file
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
                raw_input_type="fasta"
                # input_msa_paths and constraints will be empty for FASTA
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
        """
        sequences_info: List[SequenceInfo] = []
        input_msa_paths: Dict[str, str] = {}
        model_seeds: Optional[List[int]] = None
        bonded_atom_pairs: Optional[List[Any]] = None
        has_msa_flag = False # Flag if any MSA path is found

        try:
            with open(file_path, 'r') as f:
                data = json.load(f)

            # Basic validation using Pydantic model if desired (optional but recommended)
            try:
                # Pass Path object to Pydantic model if it expects it
                af3_input_model = Af3Input(**data)
                # Use validated data from the model
                data = af3_input_model.model_dump(by_alias=True) 
            except Exception as pydantic_error: # Catch validation errors specifically
                 logger.warning(f"Input JSON {file_path} failed Pydantic validation for Af3Input: {pydantic_error}. Attempting to parse leniently.")
                 # Continue with raw dictionary `data`, but parsing might fail later

            job_name_from_json = data.get("name", name_stem)
            model_seeds = data.get("modelSeeds")
            bonded_atom_pairs = data.get("bondedAtomPairs")

            chain_id_generator = idgen() # Needed if IDs are missing?
            processed_ids = set()

            for entity in data.get("sequences", []):
                if not isinstance(entity, dict) or len(entity) != 1:
                    logger.warning(f"Skipping invalid entity in AF3 JSON sequences: {entity}")
                    continue
                
                entity_type, entity_data = list(entity.items())[0]
                if not isinstance(entity_data, dict):
                     logger.warning(f"Skipping invalid entity data for type '{entity_type}': {entity_data}")
                     continue

                raw_ids = entity_data.get("id")
                chain_ids: List[str]
                if isinstance(raw_ids, str):
                    chain_ids = [raw_ids]
                elif isinstance(raw_ids, list):
                    chain_ids = [str(i) for i in raw_ids if isinstance(i, str)]
                else:
                    logger.warning(f"Entity type '{entity_type}' missing or invalid 'id'. Assigning automatically.")
                    # This case might indicate an issue, AF3 format expects IDs.
                    # Assigning one for robustness, but logging warning.
                    try:
                        chain_ids = [next(chain_id_generator)]
                        logger.warning(f"Assigned automatic ID: {chain_ids[0]}")
                    except StopIteration:
                         logger.error(f"Ran out of chain IDs while parsing AF3 JSON: {file_path}")
                         return None
                
                # Determine sequence and type
                seq: Optional[str] = None
                seq_type: SequenceType = "unknown"
                msa_path: Optional[str] = None
                original_name = f"{entity_type}_{'_'.join(chain_ids)}" # Construct a name

                if entity_type == "protein":
                    seq = entity_data.get("sequence")
                    seq_type = "protein"
                    msa_path = entity_data.get("unpairedMsaPath")
                elif entity_type == "rna":
                    seq = entity_data.get("sequence")
                    seq_type = "rna"
                    # msa_path = entity_data.get("unpairedMsaPath") # If DNA MSAs become supported
                elif entity_type == "dna":
                    seq = entity_data.get("sequence")
                    seq_type = "dna"
                    # msa_path = entity_data.get("unpairedMsaPath") # If DNA MSAs become supported
                elif entity_type == "ligand":
                    if "ccdCodes" in entity_data and isinstance(entity_data["ccdCodes"], list) and entity_data["ccdCodes"]:
                        seq = entity_data["ccdCodes"][0] # Take the first CCD code as the sequence identifier
                        seq_type = "ligand_ccd"
                    elif "smiles" in entity_data and isinstance(entity_data["smiles"], str):
                        seq = entity_data["smiles"]
                        seq_type = "ligand_smiles"
                    else:
                         logger.warning(f"Ligand entity has invalid/missing ccdCodes or smiles: {entity_data}. Skipping.")
                         continue
                else:
                    logger.warning(f"Unsupported entity type '{entity_type}' in AF3 JSON. Skipping.")
                    continue

                if seq is None:
                     logger.warning(f"Entity type '{entity_type}' with ID(s) '{chain_ids}' is missing sequence data. Skipping.")
                     continue
                
                # Add sequence info for each ID (handles homomers)
                for chain_id in chain_ids:
                    if chain_id in processed_ids:
                         logger.warning(f"Duplicate chain ID '{chain_id}' encountered in AF3 JSON. Check input format. Skipping duplicate.")
                         continue
                    
                    # Use as_entity for basic validation/normalization of sequence based on type?
                    # seq_info = as_entity(seq, chain_id, original_name)
                    # if seq_info.molecule_type != seq_type and seq_type != "unknown":
                    #      logger.warning(f"Guessed type {seq_info.molecule_type} differs from JSON type {seq_type} for ID {chain_id}")
                    # For now, trust the JSON type directly:
                    seq_info = SequenceInfo(original_name=original_name, sequence=seq, molecule_type=seq_type, chain_id=chain_id)
                    sequences_info.append(seq_info)
                    
                    if msa_path and seq_type in ["protein", "rna"]:
                        # Ensure path is absolute or relative to JSON location
                        abs_msa_path = Path(msa_path) 
                        if not abs_msa_path.is_absolute():
                            abs_msa_path = file_path.parent / msa_path
                        input_msa_paths[chain_id] = str(abs_msa_path.resolve())
                        has_msa_flag = True
                    processed_ids.add(chain_id)
            
            if not sequences_info:
                logger.error(f"No valid sequences could be extracted from AF3 JSON: {file_path}")
                return None

            job_input = JobInput(
                name_stem=job_name_from_json,
                sequences=sequences_info,
                raw_input_type="af3_json",
                input_msa_paths=input_msa_paths,
                constraints=None, # AF3 JSON doesn't have constraints like Boltz
                # Store seeds and bonds if they exist
                **({"model_seeds": model_seeds} if model_seeds else {}),
                **({"bonded_atom_pairs": bonded_atom_pairs} if bonded_atom_pairs else {})
            )
            # Set has_msa flag based on whether paths were found
            if has_msa_flag:
                job_input.__dict__["has_msa"] = True # Add flag if not default
            # Check if the input file looks like a full AF3 _data.json (contains MSA strings)
            # This is a simple heuristic: check if ANY protein has a non-empty unpairedMsa string
            if any(
                isinstance(ent, dict) and ent.get("protein", {}).get("unpairedMsa") 
                for ent in data.get("sequences", [])
            ):
                 job_input.__dict__["af3_data_json"] = str(file_path.resolve()) 
                 job_input.__dict__["has_msa"] = True # Has MSA *content*
                 logger.info("Input JSON appears to contain MSA content (unpairedMsa fields found). Setting af3_data_json path.")

            return job_input

        except json.JSONDecodeError as e:
            logger.error(f"Invalid JSON format in file {file_path}: {e}")
            return None
        except Exception as e:
            logger.error(f"Error parsing AlphaFold3 JSON file {file_path}: {e}", exc_info=True)
            return None

    def _parse_boltz_yaml(self, file_path: Path, name_stem: str) -> Optional[JobInput]:
        """
        Parses a Boltz YAML file.
        Extracts sequences, IDs, types, MSA paths, and constraints.
        """
        sequences_info: List[SequenceInfo] = []
        input_msa_paths: Dict[str, str] = {}
        constraints: Optional[List[Dict[str, Any]]] = None
        has_msa_flag = False

        try:
            with open(file_path, 'r') as f:
                data = yaml.safe_load(f)

            if not isinstance(data, dict):
                logger.error(f"Invalid YAML format: Root is not a dictionary in {file_path}")
                return None

            constraints = data.get("constraints")
            if constraints is not None and not isinstance(constraints, list):
                logger.warning(f"Ignoring invalid 'constraints' format (expected list) in {file_path}")
                constraints = None

            processed_ids = set()
            for entity in data.get("sequences", []):
                if not isinstance(entity, dict) or len(entity) != 1:
                    logger.warning(f"Skipping invalid entity in Boltz YAML sequences: {entity}")
                    continue
                
                entity_type, entity_data = list(entity.items())[0]
                if not isinstance(entity_data, dict):
                     logger.warning(f"Skipping invalid entity data for type '{entity_type}': {entity_data}")
                     continue
                
                chain_id = entity_data.get("id")
                if not chain_id or not isinstance(chain_id, str):
                    logger.warning(f"Entity type '{entity_type}' missing or invalid 'id'. Skipping.")
                    continue
                
                if chain_id in processed_ids:
                    logger.warning(f"Duplicate chain ID '{chain_id}' encountered in Boltz YAML. Skipping duplicate.")
                    continue

                seq: Optional[str] = None
                seq_type: SequenceType = "unknown"
                msa_path: Optional[str] = None
                original_name = f"{entity_type}_{chain_id}"

                if entity_type == "protein":
                    seq = entity_data.get("sequence")
                    seq_type = "protein"
                    msa_path = entity_data.get("msa") # Boltz uses "msa" key
                elif entity_type == "rna":
                    seq = entity_data.get("sequence")
                    seq_type = "rna"
                    # msa_path = entity_data.get("msa") # If Boltz supports RNA MSA
                elif entity_type == "dna":
                    seq = entity_data.get("sequence")
                    seq_type = "dna"
                elif entity_type == "ligand":
                    if "ccd" in entity_data and isinstance(entity_data["ccd"], str):
                        seq = entity_data["ccd"]
                        seq_type = "ligand_ccd"
                    elif "smiles" in entity_data and isinstance(entity_data["smiles"], str):
                        seq = entity_data["smiles"]
                        seq_type = "ligand_smiles"
                    else:
                         logger.warning(f"Ligand entity has invalid/missing ccd or smiles: {entity_data}. Skipping.")
                         continue
                else:
                    logger.warning(f"Unsupported entity type '{entity_type}' in Boltz YAML. Skipping.")
                    continue

                if seq is None:
                     logger.warning(f"Entity type '{entity_type}' with ID '{chain_id}' is missing sequence data. Skipping.")
                     continue

                seq_info = SequenceInfo(original_name=original_name, sequence=seq, molecule_type=seq_type, chain_id=chain_id)
                sequences_info.append(seq_info)
                
                if msa_path and seq_type == "protein":
                    abs_msa_path = Path(msa_path)
                    if not abs_msa_path.is_absolute():
                        abs_msa_path = file_path.parent / msa_path
                    input_msa_paths[chain_id] = str(abs_msa_path.resolve())
                    has_msa_flag = True
                processed_ids.add(chain_id)

            if not sequences_info:
                logger.error(f"No valid sequences could be extracted from Boltz YAML: {file_path}")
                return None

            job_input = JobInput(
                name_stem=name_stem,
                sequences=sequences_info,
                raw_input_type="boltz_yaml",
                input_msa_paths=input_msa_paths,
                constraints=constraints
            )
            if has_msa_flag:
                job_input.__dict__["has_msa"] = True # Add flag if not default
                
            # Boltz YAML is less likely to contain full MSA *content* like AF3 data.json
            # So we don't set af3_data_json path here.
            
            return job_input

        except yaml.YAMLError as e:
            logger.error(f"Invalid YAML format in file {file_path}: {e}")
            return None
        except Exception as e:
            logger.error(f"Error parsing Boltz YAML file {file_path}: {e}", exc_info=True)
            return None

# Example Usage (commented out)
# if __name__ == '__main__':
#     handler = InputHandler()
#     # Test FASTA parsing
#     dummy_fasta_content = (
#         ">seq1 protein\n"
#         "MPEPTIDEASE\n"
#         ">seq2 rna\n"
#         "ACGUACGU\n"
#     )
#     with open("dummy.fasta", "w") as f: f.write(dummy_fasta_content)
#     job_data_fasta = handler.parse_input_file("dummy.fasta")
#     if job_data_fasta: print(f"FASTA Parsed: {job_data_fasta.name_stem}, {len(job_data_fasta.sequences)} seqs, type: {job_data_fasta.raw_input_type}")
#     os.remove("dummy.fasta")

#     # Test AF3 JSON parsing
#     dummy_af3_json_content = {
#         "name": "test_complex",
#         "modelSeeds": [123],
#         "sequences": [
#             {"protein": {"id": "A", "sequence": "MPEPTIDE"}},
#             {"protein": {"id": "B", "sequence": "ASEQW", "unpairedMsaPath": "./msa_B.a3m"}},
#             {"rna": {"id": "C", "sequence": "ACGU"}},
#             {"ligand": {"id": "X", "smiles": "CC(=O)O"}}
#         ],
#         "dialect": "alphafold3",
#         "version": 1
#     }
#     with open("dummy.json", "w") as f: json.dump(dummy_af3_json_content, f)
#     job_data_json = handler.parse_input_file("dummy.json")
#     if job_data_json:
#         print(f"AF3 JSON Parsed: {job_data_json.name_stem}, {len(job_data_json.sequences)} seqs, type: {job_data_json.raw_input_type}")
#         for seq in job_data_json.sequences:
#             print(f"  ID: {seq.chain_id}, Type: {seq.molecule_type}, Seq: {seq.sequence[:10]}...")
#         if job_data_json.input_msa_paths:
#             print(f"  Input MSA Paths: {job_data_json.input_msa_paths}")
#     os.remove("dummy.json") 