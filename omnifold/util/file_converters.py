import json
import requests
import time
import functools
from pathlib import Path
import logging
from typing import Optional, Tuple, Dict, List, Any, Set
from .definitions import JobInput

logger = logging.getLogger(__name__)

# Create a global session for connection pooling and better performance
_session = requests.Session()

#  glycan ring codes verified to be handled correctly by Chai
GLYCAN_RING_CODES: Set[str] = {
    "A2G", "A2M", "ARA", "BMA", "FRU",
    "FUC", "GAL", "GLC", "KDN", "MAN",
    "NAG", "NDG", "RIB", "SIA", "XYS"
}

# Cache for CCD SMILES strings to avoid redundant API calls
# Increased to 1024 to better handle large batch runs
@functools.lru_cache(maxsize=1024)
def fetch_ccd_smiles(ccd_code: str) -> Optional[str]:
    """
    Fetch the canonical SMILES string for a CCD code from the RCSB PDB API.
    Includes retry logic, proper error handling, and checking multiple descriptor fields.
    
    Args:
        ccd_code: The Chemical Component Dictionary code (e.g., "HEM")
        
    Returns:
        The canonical SMILES string or None if not found/error
    """
    # Ensure CCD code is uppercase (PDB codes are standardized as uppercase)
    ccd_code = ccd_code.upper()
    
    # Retry parameters
    max_retries = 3
    
    # URL for the RCSB API
    url = f"https://data.rcsb.org/rest/v1/core/chemcomp/{ccd_code}"
    
    for attempt in range(max_retries):
        try:
            # Make the request with timeout using the session object
            response = _session.get(url, timeout=15)
            
            # Raise an exception for HTTP errors
            response.raise_for_status()
            
            # Parse the response
            data = response.json()
            
            # Check multiple descriptor fields in order of preference
            # Using the correct field names from the API
            descriptors = data.get('rcsb_chem_comp_descriptor', {})
            
            # Try different SMILES fields in order of preference (with correct API field names)
            for field in ['isomeric_smiles', 'smiles', 'smiles_stereo']:
                if field in descriptors and descriptors[field]:
                    logger.info(f"Successfully fetched SMILES for CCD '{ccd_code}' from RCSB API (field: {field})")
                    return descriptors[field]
            
            # If we get here, no SMILES field was found
            logger.warning(f"No SMILES found in any descriptor field for CCD '{ccd_code}'")
            return None
            
        except requests.HTTPError as http_err:
            # Extract status code directly from the HTTPError object
            status = http_err.response.status_code
            if status == 404:
                # CCD code not found - no need to retry
                logger.warning(f"CCD code '{ccd_code}' not found in RCSB database (404)")
                return None
            else:
                # Server error, might be transient
                logger.warning(f"HTTP error fetching data for CCD '{ccd_code}': {http_err}. Attempt {attempt+1}/{max_retries}")
                
        except (requests.ConnectionError, requests.Timeout) as err:
            # Network-related errors
            logger.warning(f"Network error fetching data for CCD '{ccd_code}': {err}. Attempt {attempt+1}/{max_retries}")
            
        except Exception as e:
            # Unexpected errors
            logger.error(f"Unexpected error fetching SMILES for CCD '{ccd_code}': {e}")
            return None
            
        # Only sleep if we're going to retry
        if attempt < max_retries - 1:
            # Simpler exponential backoff
            delay = 2 ** attempt
            time.sleep(delay)
    
    # If we've exhausted all retries
    logger.error(f"Failed to fetch SMILES for CCD '{ccd_code}' after {max_retries} attempts")
    return None

def af3_json_to_chai_fasta(json_path: str | Path, fasta_path: str | Path) -> Tuple[bool, Dict[str, Any]]:
    """
    Convert an AlphaFold-3 input JSON or _data.json file to a Chai-style FASTA.
    Handles proteins, RNA, DNA, and ligands (both SMILES and CCD codes).
    For CCD codes, it fetches SMILES strings from the RCSB PDB API.
    If multiple IDs are provided for a single CCD code entry, it creates a FASTA entry for each ID.

    Args:
        json_path: Path to the input AF3 JSON file.
        fasta_path: Path where the output Chai-style FASTA file will be written.
    
    Returns:
        A tuple containing (success: bool, result_info: Dict) where result_info contains
        statistics like count of entries by type, skipped entries, and potential warnings.
    """
    result_info = {
        "success": False,
        "total_entries": 0,
        "entries_by_type": {
            "protein": 0,
            "rna": 0,
            "dna": 0, 
            "ligand": 0,
            "glycan": 0
        },
        "skipped_entries": [],
        "warnings": []
    }

    try:
        json_path_obj = Path(json_path)
        fasta_path_obj = Path(fasta_path)

        if not json_path_obj.is_file():
            error_msg = f"AF3 JSON file not found: {json_path_obj}"
            logger.error(error_msg)
            result_info["warnings"].append(error_msg)
            return False, result_info

        j = json.loads(json_path_obj.read_text(encoding="utf-8"))
        fastas_data = [] # Store tuples of (header, sequence_or_smiles)

        if "sequences" not in j:
            error_msg = f"'sequences' key not found in AF3 JSON: {json_path_obj}"
            logger.error(error_msg)
            result_info["warnings"].append(error_msg)
            return False, result_info

        for entry in j["sequences"]:
            entity_type = None
            details_block = None

            # Determine the entity type and get its details block
            for potential_type in ("protein", "rna", "dna", "ligand"):
                if potential_type in entry:
                    entity_type = potential_type
                    details_block = entry[entity_type]
                    break
            
            if not entity_type or not details_block:
                warning_msg = f"Unknown entity type or malformed entry in AF3 JSON: {entry}"
                logger.warning(warning_msg)
                result_info["warnings"].append(warning_msg)
                result_info["skipped_entries"].append(f"unknown_type:{entry}")
                continue

            chain_ids = details_block.get("id")
            if not chain_ids:
                warning_msg = f"Entry for '{entity_type}' missing 'id'. Entry: {entry}"
                logger.warning(warning_msg)
                result_info["warnings"].append(warning_msg)
                result_info["skipped_entries"].append(f"{entity_type}:missing_id")
                continue

            if not isinstance(chain_ids, list):
                chain_ids = [chain_ids]

            if entity_type in ["protein", "rna", "dna"]:
                sequence = details_block.get("sequence")
                if not sequence:
                    warning_msg = f"Entry for '{entity_type}' missing 'sequence'. Entry: {entry}"
                    logger.warning(warning_msg)
                    result_info["warnings"].append(warning_msg)
                    result_info["skipped_entries"].append(f"{entity_type}:missing_sequence")
                    continue

                for cid in chain_ids:
                    # Chai FASTA format for polymers: >protein|name=A or >rna|name=A etc.
                    header = f"{entity_type}|name={cid}"
                    fastas_data.append((header, sequence))
                    result_info["entries_by_type"][entity_type] += 1
                    result_info["total_entries"] += 1

            elif entity_type == "ligand":
                smiles_str = details_block.get("smiles")
                ccd_codes = details_block.get("ccdCodes")

                if smiles_str: # Direct SMILES provided
                    if not isinstance(smiles_str, str):
                        warning_msg = f"Ligand entry for ID(s) {chain_ids} has non-string SMILES: {smiles_str}. Skipping."
                        logger.warning(warning_msg)
                        result_info["warnings"].append(warning_msg)
                        result_info["skipped_entries"].append(f"ligand:invalid_smiles")
                        continue

                    for cid in chain_ids:
                        header = f"ligand|name={cid}"
                        fastas_data.append((header, smiles_str))
                        result_info["entries_by_type"]["ligand"] += 1
                        result_info["total_entries"] += 1

                elif ccd_codes:
                    if not isinstance(ccd_codes, list) or not ccd_codes:
                        warning_msg = f"Ligand entry for ID(s) {chain_ids} has invalid 'ccdCodes': {ccd_codes}. Skipping."
                        logger.warning(warning_msg)
                        result_info["warnings"].append(warning_msg)
                        result_info["skipped_entries"].append(f"ligand:invalid_ccd_codes")
                        continue
                    
                    # Case 1: Many IDs, one CCD (e.g., "id": ["G","H","I"], "ccdCodes": ["ATP"])
                    # Create multiple entries with the same SMILES
                    if len(ccd_codes) == 1:
                        ccd_code = ccd_codes[0].upper()
                        # Special handling for glycan rings - create glycan records instead of ligand records with SMILES
                        if ccd_code in GLYCAN_RING_CODES:
                            logger.info(f"Treating CCD '{ccd_code}' as a glycan for chain(s) {chain_ids}")
                            for cid in chain_ids:
                                header = f"glycan|name={cid}"
                                # For glycans, we just use the CCD code as the content
                                fastas_data.append((header, ccd_code))
                                result_info["entries_by_type"]["glycan"] += 1
                                result_info["total_entries"] += 1
                        else:
                            # Normal ligand handling with SMILES
                            smiles = fetch_ccd_smiles(ccd_code)
                            if smiles:
                                for cid in chain_ids:
                                    header = f"ligand|name={cid}"
                                    fastas_data.append((header, smiles))
                                    result_info["entries_by_type"]["ligand"] += 1
                                    result_info["total_entries"] += 1
                                    logger.info(f"Added ligand '{cid}' with SMILES from CCD '{ccd_code}'")
                            else:
                                warning_msg = f"Could not get SMILES for CCD '{ccd_code}'. Skipping ligand(s) {chain_ids}."
                                logger.warning(warning_msg)
                                result_info["warnings"].append(warning_msg)
                                result_info["skipped_entries"].append(f"ligand:missing_smiles:{ccd_code}")
                    
                    # Case 2: One ID, many CCDs (e.g., "id": "J", "ccdCodes": ["ATP","MG"])
                    # Create separate FASTA entries for each CCD with ID_CCD naming
                    elif len(chain_ids) == 1:
                        cid = chain_ids[0]
                        for ccd in ccd_codes:
                            ccd = ccd.upper()
                            if ccd in GLYCAN_RING_CODES:
                                header = f"glycan|name={cid}_{ccd}"
                                fastas_data.append((header, ccd))
                                result_info["entries_by_type"]["glycan"] += 1
                                result_info["total_entries"] += 1
                                logger.info(f"Added glycan '{header}'")
                            else:
                                smiles = fetch_ccd_smiles(ccd)
                                if smiles:
                                    # Use ID_CCD naming convention for composite ligands
                                    header = f"ligand|name={cid}_{ccd}"
                                    fastas_data.append((header, smiles))
                                    result_info["entries_by_type"]["ligand"] += 1
                                    result_info["total_entries"] += 1
                                    logger.info(f"Added composite ligand component '{cid}_{ccd}' with SMILES from CCD '{ccd}'")
                                else:
                                    warning_msg = f"Could not get SMILES for CCD '{ccd}' in composite ligand '{cid}'. Skipping this component."
                                    logger.warning(warning_msg)
                                    result_info["warnings"].append(warning_msg)
                                    result_info["skipped_entries"].append(f"ligand:missing_smiles:{ccd}")
                    
                    # Case 3: Many IDs, many CCDs of same length (e.g., "id": ["A","B"], "ccdCodes": ["ATP","ADP"])
                    # Zip them together for 1:1 mapping
                    elif len(chain_ids) == len(ccd_codes):
                        logger.info(f"Processing 1:1 mapping of chain IDs to CCD codes: {chain_ids} -> {ccd_codes}")
                        for cid, ccd in zip(chain_ids, ccd_codes):
                            ccd = ccd.upper()
                            if ccd in GLYCAN_RING_CODES:
                                header = f"glycan|name={cid}"
                                fastas_data.append((header, ccd))
                                result_info["entries_by_type"]["glycan"] += 1
                                result_info["total_entries"] += 1
                                logger.info(f"Added glycan '{header}' (1:1 mapping)")
                            else:
                                smiles = fetch_ccd_smiles(ccd)
                                if smiles:
                                    header = f"ligand|name={cid}"
                                    fastas_data.append((header, smiles))
                                    result_info["entries_by_type"]["ligand"] += 1
                                    result_info["total_entries"] += 1
                                    logger.info(f"Added ligand '{cid}' with SMILES from CCD '{ccd}' (1:1 mapping)")
                                else:
                                    warning_msg = f"Could not get SMILES for CCD '{ccd}' mapped to ID '{cid}'. Skipping this ligand."
                                    logger.warning(warning_msg)
                                    result_info["warnings"].append(warning_msg)
                                    result_info["skipped_entries"].append(f"ligand:missing_smiles:{ccd}")
                    
                    # Case 4: Ambiguous mapping (e.g., "id": ["A","B","C"], "ccdCodes": ["ATP","ADP"])
                    # Log error for ambiguous mapping
                    else:
                        warning_msg = f"Ambiguous ligand mapping: ids={chain_ids}, ccdCodes={ccd_codes}. Cannot determine correct pairing. Skipping these ligands."
                        logger.error(warning_msg)
                        result_info["warnings"].append(warning_msg)
                        result_info["skipped_entries"].append(f"ligand:ambiguous_mapping")
                else:
                    warning_msg = f"Ligand entry for ID(s) {chain_ids} has no 'smiles' or 'ccdCodes'. Entry: {entry}"
                    logger.warning(warning_msg)
                    result_info["warnings"].append(warning_msg)
                    result_info["skipped_entries"].append(f"ligand:no_definition")
            else:
                warning_msg = f"Unrecognized entity_type '{entity_type}' encountered. Entry: {entry}"
                logger.warning(warning_msg)
                result_info["warnings"].append(warning_msg)
                result_info["skipped_entries"].append(f"unknown_entity:{entity_type}")

        if not fastas_data:
            warning_msg = f"No valid sequence or ligand entries found in AF3 JSON {json_path_obj} to convert to Chai FASTA."
            logger.warning(warning_msg)
            result_info["warnings"].append(warning_msg)
            return False, result_info

        fasta_string = "".join(f">{header}\n{seq_or_smiles}\n" for header, seq_or_smiles in fastas_data)
        
        fasta_path_obj.parent.mkdir(parents=True, exist_ok=True)
        fasta_path_obj.write_text(fasta_string, encoding="utf-8")
        logger.info(f"Wrote {len(fastas_data)} Chai-style FASTA records from {json_path_obj} to {fasta_path_obj}")
        
        result_info["success"] = True
        return True, result_info
        
    except json.JSONDecodeError as e:
        error_msg = f"Error decoding AF3 JSON file {json_path}: {e}"
        logger.error(error_msg)
        result_info["warnings"].append(error_msg)
        return False, result_info
    except IOError as e:
        error_msg = f"Error writing Chai FASTA file to {fasta_path}: {e}"
        logger.error(error_msg)
        result_info["warnings"].append(error_msg)
        return False, result_info
    except Exception as e:
        error_msg = f"An unexpected error occurred during AF3 JSON to Chai FASTA conversion: {e}"
        logger.error(error_msg, exc_info=True)
        result_info["warnings"].append(error_msg)
        return False, result_info 

def job_input_to_chai_fasta(job_input: JobInput, fasta_path: str | Path) -> bool:
    """
    Generates a Chai-style FASTA file directly from a JobInput object.
    This ensures consistent FASTA formatting for Chai-1, regardless of the MSA source.
    """
    fasta_path_obj = Path(fasta_path)
    try:
        with open(fasta_path_obj, "w") as f:
            for seq_info in job_input.sequences:
                # Determine the entity type for the header.
                # Default to 'protein' if it's a polymer, otherwise use the specific type.
                entity_type = seq_info.molecule_type
                if entity_type not in ["ligand", "glycan", "protein", "rna", "dna"]:
                    # Fallback for unknown polymer types if any, default to protein
                    entity_type = "protein"

                header = f">{entity_type}|name={seq_info.chain_id}"
                # The sequence for ligands is their SMILES or CCD code
                f.write(f"{header}\n")
                f.write(f"{seq_info.sequence}\n")
        
        logger.info(f"Successfully generated Chai-style FASTA from JobInput at: {fasta_path_obj}")
        return True
    except Exception as e:
        logger.error(f"Failed to generate Chai-style FASTA from JobInput: {e}", exc_info=True)
        return False 