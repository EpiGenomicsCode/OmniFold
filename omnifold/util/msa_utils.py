import json
import logging
import os
from typing import Dict, List, Optional, Union

logger = logging.getLogger(__name__)

def extract_all_protein_a3ms_from_af3_json(json_path: str, output_dir: str) -> Optional[Dict[str, str]]:
    """
    Extracts the unpaired MSA (A3M format) for all protein sequences
    found in an AlphaFold3 data.json file and writes them to individual files
    in the specified output directory.

    Handles both single IDs and lists of IDs (for homomers), assuming the same
    MSA applies to all IDs in a list.

    Args:
        json_path: Path to the input AlphaFold3 data.json file.
        output_dir: Directory where the extracted A3M files (e.g., msa_A.a3m) should be written.

    Returns:
        A dictionary mapping each protein ID (str) to the absolute path (str)
        of its extracted A3M file, or None if a critical error occurs.
        Partial results might be written even if warnings occur for some proteins.
    """
    logger.info(f"Attempting to extract all protein A3Ms from AF3 JSON: {json_path} into {output_dir}")
    extracted_paths: Dict[str, str] = {}
    has_critical_error = False

    try:
        with open(json_path, 'r') as f:
            data = json.load(f)

        if not isinstance(data, dict) or "sequences" not in data:
            logger.error(f"Invalid AF3 JSON format in {json_path}: Missing 'sequences' key.")
            return None 

        os.makedirs(output_dir, exist_ok=True)

        sequences = data.get("sequences", [])
        if not sequences:
             logger.warning(f"No sequences found in {json_path}.")
             return {} 

        protein_found = False
        for entity_index, entity in enumerate(sequences):
            if isinstance(entity, dict) and "protein" in entity:
                protein_found = True
                protein_data = entity["protein"]
                
                if not isinstance(protein_data, dict):
                    logger.warning(f"Protein entry at index {entity_index} is not a dictionary. Skipping.")
                    continue

                protein_ids_raw = protein_data.get('id')
                msa_content = protein_data.get("unpairedMsa")

                if not protein_ids_raw:
                    logger.warning(f"Protein entry at index {entity_index} missing 'id'. Skipping.")
                    continue
                    
                if msa_content is None:
                    logger.info(f"Protein ID(s) {protein_ids_raw} has no 'unpairedMsa' field. Skipping extraction for this entity.")
                    continue 

                if not isinstance(msa_content, str):
                    logger.warning(f"Protein ID(s) {protein_ids_raw} has 'unpairedMsa' but it's not a string. Skipping extraction for this entity.")
                    continue
                    
                if not msa_content.strip():
                     logger.warning(f"Protein ID(s) {protein_ids_raw} has an empty or whitespace-only 'unpairedMsa'. Writing empty file(s)." )

                protein_ids: List[str]
                if isinstance(protein_ids_raw, str):
                    protein_ids = [protein_ids_raw]
                elif isinstance(protein_ids_raw, list):
                    protein_ids = [str(pid) for pid in protein_ids_raw if isinstance(pid, str)]
                    if len(protein_ids) != len(protein_ids_raw):
                        logger.warning(f"Non-string IDs found in list for protein entity {entity_index}: {protein_ids_raw}. Skipping non-string IDs.")
                else:
                    logger.warning(f"Protein ID for entity {entity_index} has unexpected type: {type(protein_ids_raw)}. Skipping.")
                    continue
                
                if not protein_ids:
                    logger.warning(f"No valid string IDs found for protein entity {entity_index}. Skipping.")
                    continue

                for protein_id in protein_ids:
                    output_a3m_filename = f"msa_{protein_id}.a3m"
                    output_a3m_path = os.path.abspath(os.path.join(output_dir, output_a3m_filename))
                    
                    try:
                        with open(output_a3m_path, 'w') as f:
                            f.write(msa_content)
                        logger.info(f"Successfully extracted A3M for Protein ID '{protein_id}' to: {output_a3m_path}")
                        if protein_id in extracted_paths:
                             logger.warning(f"Duplicate Protein ID '{protein_id}' encountered. Overwriting previous A3M path.")
                        extracted_paths[protein_id] = output_a3m_path
                    except IOError as e:
                        logger.error(f"Error writing extracted A3M file for ID '{protein_id}' to {output_a3m_path}: {e}")
                        has_critical_error = True 
                        break 
                if has_critical_error: break 
                        
        if not protein_found:
             logger.warning(f"No protein entities found in the 'sequences' list in {json_path}.")

        if has_critical_error:
             logger.error("Critical error occurred during A3M extraction. Returning None.")
             return None
        else:
            logger.info(f"Finished extracting A3Ms. Found {len(extracted_paths)} MSAs.")
            return extracted_paths

    except FileNotFoundError:
        logger.error(f"AF3 JSON file not found: {json_path}")
        return None
    except json.JSONDecodeError:
        logger.error(f"Error decoding JSON from file: {json_path}")
        return None
    except IOError as e:
        logger.error(f"Error creating output directory or reading input file {json_path}: {e}")
        return None
    except Exception as e:
        logger.error(f"An unexpected error occurred during A3M extraction: {e}", exc_info=True)
        return None 

def is_a3m_singleton(a3m_file_path: str, query_sequence: str) -> bool:
    """
    Checks if an A3M file contains only a single sequence (the query) and no other alignments.

    Args:
        a3m_file_path: Path to the A3M file.
        query_sequence: The expected query sequence (without gaps).

    Returns:
        True if the A3M is a singleton (only query), False otherwise.
    """
    try:
        with open(a3m_file_path, 'r') as f:
            lines = [line.strip() for line in f if line.strip()] 
        
        if not lines: 
            return False

        num_sequences = 0
        first_sequence_lines = []
        in_first_sequence = False

        for line in lines:
            if line.startswith('>'):
                num_sequences += 1
                if num_sequences > 1:
                    return False 
                in_first_sequence = True
            elif in_first_sequence:
                first_sequence_lines.append(line)

        if num_sequences != 1:
            return False 

        first_sequence_from_file = "".join(first_sequence_lines).replace("-", "").replace(".", "")
        
        clean_query_sequence = query_sequence.replace("-", "").replace(".", "")
        
        return first_sequence_from_file.upper() == clean_query_sequence.upper()

    except FileNotFoundError:
        logger.warning(f"A3M file not found for singleton check: {a3m_file_path}")
        return False 
    except Exception as e:
        logger.error(f"Error checking A3M file {a3m_file_path} for singleton status: {e}")
        return False 