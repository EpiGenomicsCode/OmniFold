import json
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

def af3_json_to_chai_fasta(json_path: str | Path, fasta_path: str | Path) -> bool:
    """
    Convert an AlphaFold-3 input JSON or _data.json file to a Chai-style FASTA.

    Args:
        json_path: Path to the input AF3 JSON file.
        fasta_path: Path where the output Chai-style FASTA file will be written.
    
    Returns:
        True if conversion was successful and file was written, False otherwise.
    """
    try:
        json_path_obj = Path(json_path)
        fasta_path_obj = Path(fasta_path)

        if not json_path_obj.is_file():
            logger.error(f"AF3 JSON file not found: {json_path_obj}")
            return False

        j = json.loads(json_path_obj.read_text())
        fastas_data = [] # Store tuples of (header, sequence)

        if "sequences" not in j:
            logger.error(f"'sequences' key not found in AF3 JSON: {json_path_obj}")
            return False

        for entry in j["sequences"]:
            found_key = False
            for key in ("protein", "rna", "dna", "ligand"):
                if key in entry:
                    block = entry[key]
                    # Ensure 'id' and 'sequence' keys exist
                    if "id" not in block or "sequence" not in block:
                        logger.warning(f"Entry for '{key}' missing 'id' or 'sequence'. Entry: {entry}")
                        continue # Skip this problematic block

                    ids = block["id"] if isinstance(block["id"], list) else [block["id"]]
                    seq = block["sequence"]
                    for cid in ids:
                        header = f"{key}|name=chain_{cid}"
                        fastas_data.append((header, seq))
                    found_key = True
                    break 
            if not found_key:
                logger.warning(f"Unknown polymer type or malformed entry in AF3 JSON: {entry}")

        if not fastas_data:
            logger.warning(f"No valid sequence entries found in AF3 JSON {json_path_obj} to convert to FASTA.")
            return False # Or write an empty file if that's preferred

        fasta_string = "".join(f">{header}\n{seq}\n" for header, seq in fastas_data)
        
        fasta_path_obj.parent.mkdir(parents=True, exist_ok=True)
        fasta_path_obj.write_text(fasta_string)
        logger.info(f"Wrote {len(fastas_data)} Chai-style FASTA records from {json_path_obj} to {fasta_path_obj}")
        return True
        
    except json.JSONDecodeError as e:
        logger.error(f"Error decoding AF3 JSON file {json_path}: {e}")
        return False
    except IOError as e:
        logger.error(f"Error writing Chai FASTA file to {fasta_path}: {e}")
        return False
    except Exception as e:
        logger.error(f"An unexpected error occurred during AF3 JSON to Chai FASTA conversion: {e}", exc_info=True)
        return False 