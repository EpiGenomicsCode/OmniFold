import os
import json # For loading JSON file content
import pyfastx
from typing import List, Optional, Dict, Union

from .util.definitions import JobInput, SequenceInfo, idgen, as_entity, SequenceType
from .af3_models import Af3Input as Af3PydanticModel # Rename to avoid confusion with function
from .af3_models import Protein, RNA, DNA, Ligand # To identify sequence item types

class InputHandler:
    def __init__(self):
        """Initializes the InputHandler."""
        pass

    def parse_input_file(self, file_path: str) -> Optional[JobInput]:
        """
        Parses the input file, determines its type, and processes accordingly.
        """
        if not os.path.exists(file_path):
            print(f"Error: Input file not found at {file_path}")
            return None
            
        file_ext = os.path.splitext(file_path)[1].lower()

        if file_ext in ['.fasta', '.fa', '.fna', '.fsa']:
            return self._parse_fasta_file(file_path)
        elif file_ext == '.json':
            return self._parse_af3_json_file(file_path)
        # elif file_ext in ['.yaml', '.yml']:
        #     # Placeholder for Boltz YAML input
        #     print(f"Boltz YAML input parsing not yet implemented for {file_path}")
        #     return None
        else:
            print(f"Error: Unrecognized file type for {file_path}")
            return None

    def _parse_fasta_file(self, fasta_path: str) -> Optional[JobInput]:
        """
        Parses a FASTA file and converts its content into a JobInput object.
        Uses pyfastx for parsing and utility functions for type guessing and ID generation.
        """
        if not os.path.exists(fasta_path):
            print(f"Error: FASTA file not found at {fasta_path}")
            return None

        try:
            # Extract a name stem from the FASTA file path
            name_stem = os.path.splitext(os.path.basename(fasta_path))[0]
            
            parsed_sequences: List[SequenceInfo] = []
            chain_id_generator = idgen() # Initialize chain ID generator

            for record in pyfastx.Fasta(fasta_path):
                name = record.name
                sequence = record.seq
                
                current_chain_id = next(chain_id_generator)
                
                # Guess sequence type and create SequenceInfo object
                # The as_entity function is responsible for creating the SequenceInfo object
                seq_info = as_entity(seq=sequence, chain_id=current_chain_id, original_name=name)
                parsed_sequences.append(seq_info)
            
            if not parsed_sequences:
                print(f"Warning: No sequences found or parsed in {fasta_path}")
                return None

            job_input = JobInput(
                name_stem=name_stem,
                sequences=parsed_sequences,
                raw_input_type="fasta"
            )
            return job_input

        except Exception as e:
            print(f"Error parsing FASTA file {fasta_path}: {e}")
            return None

    def _parse_af3_json_file(self, json_path: str) -> Optional[JobInput]:
        """
        Parses an AlphaFold3 JSON file into a JobInput object.
        Extracts sequences, chain IDs, and any provided MSA paths.
        """
        try:
            with open(json_path, 'r') as f:
                data = json.load(f)
            
            # Validate and parse using Pydantic model
            parsed_af3_input = Af3PydanticModel.model_validate(data)

            name_stem = parsed_af3_input.name
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
                    print(f"Warning: Encountered an unexpected sequence item type in {json_path}")

            if not job_sequences:
                print(f"Warning: No sequences extracted from AF3 JSON {json_path}")
                return None

            return JobInput(
                name_stem=name_stem,
                sequences=job_sequences,
                raw_input_type="af3_json",
                input_msa_paths=input_msa_paths if input_msa_paths else None
            )

        except Exception as e:
            # Covers JSON decoding errors, Pydantic validation errors, etc.
            print(f"Error parsing AlphaFold3 JSON file {json_path}: {e}")
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