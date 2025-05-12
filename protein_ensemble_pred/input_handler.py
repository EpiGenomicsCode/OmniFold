import os
import json # For loading JSON file content
import pyfastx
from typing import List, Optional, Dict, Union
from pathlib import Path
import yaml # Added import

from .util.definitions import JobInput, SequenceInfo, idgen, as_entity, SequenceType
from .af3_models import Af3Input as Af3PydanticModel # Rename to avoid confusion with function
from .af3_models import Protein, RNA, DNA, Ligand # To identify sequence item types

class InputHandler:
    def __init__(self):
        """Initializes the InputHandler."""
        pass

    def parse_input_file(self, filepath_str: str) -> JobInput:
        """
        Parses the input file (FASTA, AF3 JSON, or Boltz YAML) and returns a JobInput object.
        """
        filepath = Path(filepath_str)
        if not filepath.exists():
            raise FileNotFoundError(f"Input file not found: {filepath}")

        name_stem = filepath.stem
        file_extension = filepath.suffix.lower()

        if file_extension in [".fasta", ".fas", ".fa", ".fna", ".ffn", ".faa", ".frn"]:
            return self._parse_fasta_file(filepath, name_stem)
        elif file_extension == ".json":
            # Could be AF3 JSON. Add more checks if other JSON types are expected.
            return self._parse_af3_json_file(filepath, name_stem)
        elif file_extension in [".yaml", ".yml"]:
            return self._parse_boltz_yaml_file(filepath, name_stem)
        else:
            raise ValueError(
                f"Unsupported file extension: '{file_extension}'. "
                "Please provide a FASTA (.fasta, .fa, etc.), AF3 JSON (.json), or Boltz YAML (.yaml, .yml) file."
            )

    def _parse_fasta_file(self, filepath: Path, name_stem: str) -> JobInput:
        sequences_info: List[SequenceInfo] = []
        input_msa_paths: Dict[str, str] = {}
        chain_id_generator = idgen()

        try:
            for record in pyfastx.Fasta(str(filepath)):
                original_name = record.name
                sequence = record.seq
                
                current_chain_id: Optional[str] = None
                molecule_type_explicitly_defined = False
                final_molecule_type: SequenceType = "unknown" # Default

                if "|" in original_name:
                    parts = original_name.split("|", 2)
                    parsed_chain_id = parts[0].strip()
                    # Use parsed_chain_id if non-empty, otherwise will fall to idgen later if needed
                    current_chain_id = parsed_chain_id if parsed_chain_id else None 
                    
                    entity_type_str = parts[1].lower().strip() if len(parts) > 1 else ""
                    msa_path_str = parts[2].strip() if len(parts) > 2 and parts[2].strip() else None

                    temp_molecule_type: Optional[SequenceType] = None
                    if entity_type_str == "protein":
                        temp_molecule_type = "protein"
                    elif entity_type_str == "rna":
                        temp_molecule_type = "rna"
                    elif entity_type_str == "dna":
                        temp_molecule_type = "dna"
                    elif entity_type_str == "ccd":
                        temp_molecule_type = "ligand_ccd"
                    elif entity_type_str == "smiles":
                        temp_molecule_type = "ligand_smiles"
                    
                    if temp_molecule_type:
                        final_molecule_type = temp_molecule_type
                        molecule_type_explicitly_defined = True
                        if not current_chain_id: # If parsed_chain_id was empty
                            current_chain_id = next(chain_id_generator)
                        
                        if final_molecule_type == "protein" and msa_path_str and msa_path_str.lower() != "empty":
                            input_msa_paths[current_chain_id] = msa_path_str
                    # else: explicit type not recognized, will fall through to as_entity
                
                if not molecule_type_explicitly_defined:
                    # Determine chain_id for as_entity: use parsed if available & valid, else generate
                    chain_id_for_guessing = current_chain_id if current_chain_id else next(chain_id_generator)
                    
                    # Perform type guessing
                    guessed_seq_info = as_entity(sequence, chain_id_for_guessing, original_name)
                    final_molecule_type = guessed_seq_info.molecule_type
                    current_chain_id = guessed_seq_info.chain_id # as_entity provides the definitive chain_id used
                    # Note: MSA path from header is NOT used if we fell back to type guessing, 
                    # as the explicit 'protein' type was not declared with it.

                # Ensure current_chain_id is set (should be by now)
                if not current_chain_id:
                     current_chain_id = next(chain_id_generator) # Should be very rare, a final fallback

                sequences_info.append(SequenceInfo(
                    original_name=original_name,
                    sequence=sequence, # as_entity might adjust sequence (e.g. case), ensure we use the right one
                                       # Current as_entity returns original sequence for ligand_smiles/unknown
                                       # and upper-cased for protein/rna/dna/ligand_ccd.
                                       # For consistency, let's ensure sequence added is what as_entity would have determined
                                       # or the explicit sequence if type was explicit.
                    molecule_type=final_molecule_type,
                    chain_id=current_chain_id
                ))
            # Correcting the sequence for SequenceInfo based on final_molecule_type determination logic
            # This part needs to be inside the loop, or sequences_info needs to be adjusted post-loop.
            # For simplicity, let's refine what goes into SequenceInfo directly.

        except Exception as e:
            raise ValueError(f"Error parsing FASTA file '{filepath}': {e}")

        # Refinement for sequence string based on molecule type after the loop (if necessary)
        # Or, more cleanly, ensure the sequence passed to SequenceInfo constructor is the final one.
        # The current logic for as_entity already returns the sequence in its final form (e.g., uppercase).
        # If explicit type, the original sequence is used. This might be an inconsistency.
        # Let's assume `as_entity` is the source of truth for sequence format if type guessing occurs.
        # If type is explicit, we use the raw sequence. This seems acceptable.

        final_sequences_info = []
        for i, s_info_draft in enumerate(sequences_info):
            # If type was guessed, as_entity already formatted the sequence.
            # If type was explicit, s_info_draft.sequence is the raw sequence.
            # To ensure consistency (e.g. for ligand_ccd being uppercase as per as_entity):
            final_seq_str = s_info_draft.sequence
            if s_info_draft.molecule_type in ["protein", "rna", "dna", "ligand_ccd"]:
                 final_seq_str = s_info_draft.sequence.upper()
            # For smiles and unknown, keep original casing as `as_entity` does.

            final_sequences_info.append(SequenceInfo(
                original_name=s_info_draft.original_name,
                sequence=final_seq_str, # Use consistently formatted sequence
                molecule_type=s_info_draft.molecule_type,
                chain_id=s_info_draft.chain_id,
                molecule_type_confidence=s_info_draft.molecule_type_confidence if hasattr(s_info_draft, 'molecule_type_confidence') else 1.0
            ))
        sequences_info = final_sequences_info


        if not sequences_info:
            raise ValueError(f"No sequences found in FASTA file: {filepath}")

        return JobInput(
            name_stem=name_stem,
            sequences=sequences_info,
            raw_input_type="fasta",
            input_msa_paths=input_msa_paths
        )

    def _parse_af3_json_file(self, filepath: Path, name_stem: str) -> JobInput:
        """
        Parses an AlphaFold3 JSON file into a JobInput object.
        Extracts sequences, chain IDs, and any provided MSA paths.
        """
        try:
            with open(filepath, 'r') as f:
                data = json.load(f)
            
            # Validate and parse using Pydantic model
            parsed_af3_input = Af3PydanticModel.model_validate(data)

            job_sequences: List[SequenceInfo] = []
            input_msa_paths: Dict[str, str] = {}

            chain_id_gen_for_ligands = idgen() # In case ligands in JSON don't have IDs we can use

            for seq_item_wrapper in parsed_af3_input.sequences:
                seq_info = None
                # The wrapper can be Protein, RNA, DNA, or Ligand
                # We need to get the actual chain object (ProteinChain, RNAChain, etc.)
                actual_chain_obj = None
                molecule_type_literal: Optional[SequenceType] = None
                original_name_from_json = name_stem # Default, can be improved if chain has own name

                if isinstance(seq_item_wrapper, Protein):
                    actual_chain_obj = seq_item_wrapper.protein
                    molecule_type_literal = "protein"
                elif isinstance(seq_item_wrapper, RNA):
                    actual_chain_obj = seq_item_wrapper.rna
                    molecule_type_literal = "rna"
                elif isinstance(seq_item_wrapper, DNA):
                    actual_chain_obj = seq_item_wrapper.dna
                    molecule_type_literal = "dna"
                elif isinstance(seq_item_wrapper, Ligand):
                    actual_chain_obj = seq_item_wrapper.ligand
                    # Determine if ligand is CCD or SMILES based on which field is present
                    if actual_chain_obj.ccdCodes:
                        molecule_type_literal = "ligand_ccd"
                    elif actual_chain_obj.smiles:
                        molecule_type_literal = "ligand_smiles"
                
                if actual_chain_obj and molecule_type_literal:
                    # Handle single ID or list of IDs (for homomers)
                    ids_to_process: List[str]
                    if isinstance(actual_chain_obj.id, str):
                        ids_to_process = [actual_chain_obj.id]
                    else: # It's a list of IDs
                        ids_to_process = actual_chain_obj.id
                    
                    # For molecules with actual sequences
                    sequence_str = ""
                    if hasattr(actual_chain_obj, 'sequence') and actual_chain_obj.sequence:
                        sequence_str = actual_chain_obj.sequence
                    elif molecule_type_literal == "ligand_ccd" and actual_chain_obj.ccdCodes:
                        sequence_str = actual_chain_obj.ccdCodes[0] # Taking the first CCD code as sequence
                    elif molecule_type_literal == "ligand_smiles" and actual_chain_obj.smiles:
                        sequence_str = actual_chain_obj.smiles

                    for chain_id_str in ids_to_process:
                        seq_info = SequenceInfo(
                            original_name=original_name_from_json, # Could try to get more specific name if available
                            sequence=sequence_str,
                            molecule_type=molecule_type_literal,
                            chain_id=chain_id_str
                        )
                        job_sequences.append(seq_info)

                        # Check for MSA paths (for Protein, RNA, DNA)
                        if hasattr(actual_chain_obj, 'unpairedMsaPath') and actual_chain_obj.unpairedMsaPath:
                            input_msa_paths[chain_id_str] = actual_chain_obj.unpairedMsaPath
                        # Could also check for pairedMsaPath if relevant for JobInput
                else:
                    # This branch should ideally not be hit if Pydantic models are comprehensive
                    # and input conforms to one of the wrapped types (Protein, RNA, DNA, Ligand).
                    # If it is hit, it means seq_item_wrapper was not an instance of these.
                    raise ValueError(f"Encountered an unexpected sequence item type in AF3 JSON '{filepath}'. Expected Protein, RNA, DNA, or Ligand wrapper.")

            if not job_sequences:
                # This case implies the 'sequences' list in JSON was empty or all items were unparseable.
                raise ValueError(f"No valid sequences extracted from AF3 JSON '{filepath}'.")

            return JobInput(
                name_stem=name_stem,
                sequences=job_sequences,
                raw_input_type="af3_json",
                input_msa_paths=input_msa_paths if input_msa_paths else {}
            )

        except FileNotFoundError: # Specific exception for file not found
            raise FileNotFoundError(f"AF3 JSON file not found: {filepath}")
        except json.JSONDecodeError as e:
            raise ValueError(f"Error decoding AF3 JSON file '{filepath}': Invalid JSON format. {e}")
        except Exception as e: # Covers Pydantic validation errors, other unexpected issues
            # Consider logging the original exception e for debugging, or re-raising if specific handling is not needed
            raise ValueError(f"Error parsing AlphaFold3 JSON file '{filepath}': {e}")
            # Removed return None, execution will stop due to exception

    def _parse_boltz_yaml_file(self, filepath: Path, name_stem: str) -> JobInput:
        """Parses a Boltz YAML input file."""
        sequences_info: List[SequenceInfo] = []
        input_msa_paths: Dict[str, str] = {}
        constraints: List[Dict[str, Any]] | None = None

        try:
            with open(filepath, 'r') as f:
                data = yaml.safe_load(f)
        except yaml.YAMLError as e:
            raise ValueError(f"Error parsing Boltz YAML file '{filepath}': Invalid YAML format. {e}")
        except FileNotFoundError:
            raise FileNotFoundError(f"Boltz YAML file not found: {filepath}")
        except Exception as e:
            raise ValueError(f"An unexpected error occurred while reading Boltz YAML file '{filepath}': {e}")


        if not isinstance(data, dict):
            raise ValueError(f"Invalid Boltz YAML content in '{filepath}': Expected a dictionary at the root.")

        parsed_sequences = data.get("sequences", [])
        if not isinstance(parsed_sequences, list):
            raise ValueError(f"Invalid Boltz YAML '{filepath}': 'sequences' field must be a list.")

        for i, item in enumerate(parsed_sequences):
            if not isinstance(item, dict) or len(item) != 1:
                raise ValueError(f"Invalid Boltz YAML '{filepath}': Each item in 'sequences' list (index {i}) must be a dictionary with a single key (e.g., 'protein', 'ligand').")

            entity_type_key = list(item.keys())[0] # e.g., "protein", "rna", "dna", "ligand"
            entity_data = item[entity_type_key]

            if not isinstance(entity_data, dict):
                raise ValueError(f"Invalid Boltz YAML '{filepath}': Value for '{entity_type_key}' (index {i}) must be a dictionary.")

            chain_ids = entity_data.get("id")
            if isinstance(chain_ids, str):
                chain_ids = [chain_ids] # Normalize to list
            elif not isinstance(chain_ids, list) or not all(isinstance(cid, str) for cid in chain_ids):
                raise ValueError(f"Invalid Boltz YAML '{filepath}': 'id' for '{entity_type_key}' (index {i}) must be a string or list of strings.")

            sequence_str = entity_data.get("sequence")
            smiles_str = entity_data.get("smiles")
            ccd_code = entity_data.get("ccd")
            msa_path = entity_data.get("msa")

            molecule_type: SequenceType = "unknown"
            actual_sequence_data = ""

            if entity_type_key == "protein":
                molecule_type = "protein"
                if not isinstance(sequence_str, str):
                    raise ValueError(f"Invalid Boltz YAML '{filepath}': 'sequence' missing or not a string for protein (index {i}).")
                actual_sequence_data = sequence_str
            elif entity_type_key == "rna":
                molecule_type = "rna"
                if not isinstance(sequence_str, str):
                    raise ValueError(f"Invalid Boltz YAML '{filepath}': 'sequence' missing or not a string for RNA (index {i}).")
                actual_sequence_data = sequence_str
            elif entity_type_key == "dna":
                molecule_type = "dna"
                if not isinstance(sequence_str, str):
                    raise ValueError(f"Invalid Boltz YAML '{filepath}': 'sequence' missing or not a string for DNA (index {i}).")
                actual_sequence_data = sequence_str
            elif entity_type_key == "ligand":
                if smiles_str and isinstance(smiles_str, str):
                    molecule_type = "ligand_smiles"
                    actual_sequence_data = smiles_str
                elif ccd_code and isinstance(ccd_code, str):
                    molecule_type = "ligand_ccd"
                    actual_sequence_data = ccd_code
                else:
                    raise ValueError(f"Invalid Boltz YAML '{filepath}': Ligand (index {i}) must have a 'smiles' or 'ccd' string.")
            else:
                raise ValueError(f"Unsupported entity type '{entity_type_key}' in Boltz YAML (index {i}) from file '{filepath}'.")

            for chain_id_str in chain_ids:
                sequences_info.append(SequenceInfo(
                    original_name=chain_id_str, # Use chain_id as original_name
                    sequence=actual_sequence_data,
                    molecule_type=molecule_type,
                    chain_id=chain_id_str
                ))
                if molecule_type == "protein" and msa_path and isinstance(msa_path, str) and msa_path.lower() != "empty":
                    input_msa_paths[chain_id_str] = msa_path
        
        parsed_constraints = data.get("constraints")
        if parsed_constraints is not None:
            if not isinstance(parsed_constraints, list):
                 raise ValueError(f"Invalid Boltz YAML '{filepath}': 'constraints' field, if present, must be a list.")
            constraints = parsed_constraints


        if not sequences_info:
            raise ValueError(f"No sequences found in Boltz YAML file: {filepath}")

        return JobInput(
            name_stem=name_stem,
            sequences=sequences_info,
            raw_input_type="boltz_yaml",
            input_msa_paths=input_msa_paths,
            constraints=constraints
        )

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