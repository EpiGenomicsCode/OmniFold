import os
import logging
import subprocess
from typing import List, Dict

logger = logging.getLogger(__name__)

def detect_available_gpus() -> List[int]:
    """
    Detect available GPUs and return their IDs.
    
    Returns:
        List of available GPU IDs
    """
    cuda_devices = os.environ.get("CUDA_VISIBLE_DEVICES")
    if cuda_devices:
        return [int(d) for d in cuda_devices.split(",")]
    
    try:
        import torch
        return list(range(torch.cuda.device_count()))
    except ImportError:
        try:
            result = subprocess.run(
                ["nvidia-smi", "-L"],
                capture_output=True,
                text=True,
                check=True
            )
            return list(range(len(result.stdout.strip().split("\n"))))
        except Exception as e:
            logger.error(f"Failed to detect GPUs: {e}")
            return []

def assign_gpus_to_models(requested_model_names: List[str], force_sequential: bool = False) -> Dict[str, any]:
    """
    Assign GPUs to models based on availability and sequential flag.
    
    Args:
        requested_model_names: List of model names that need GPU assignment (e.g., ['alphafold3', 'chai1'])
        force_sequential: If True, all models are assigned to the first available GPU.
        
    Returns:
        Dictionary mapping model names to GPU IDs (str, e.g. "0", "1") or a placeholder string.
    """
    all_known_models_ordered = ["alphafold3", "boltz1", "chai1"] # Defines priority or typical order
    
    available_gpus_str = [str(gpu_id) for gpu_id in detect_available_gpus()] # Ensure string IDs

    if not available_gpus_str:
        logger.warning("No GPUs available or nvidia-smi failed. Assigning all models to run sequentially on CPU/placeholder GPU.")
        return {model_name: "gpu_placeholder_0" for model_name in requested_model_names}

    num_gpus = len(available_gpus_str)
    num_models_to_run = len(requested_model_names)

    if force_sequential or num_models_to_run == 0:
        if num_models_to_run == 0:
            return {}
        logger.info(f"Assigning all {num_models_to_run} requested models sequentially to {available_gpus_str[0]}.")
        return {model_name: available_gpus_str[0] for model_name in requested_model_names}
    
    # Filter and order the requested models based on our known model list
    ordered_and_filtered_requested_models = [m for m in all_known_models_ordered if m in requested_model_names]
    
    # Add any requested models that weren't in all_known_models_ordered to the end (less prioritized)
    for model_name in requested_model_names:
        if model_name not in ordered_and_filtered_requested_models:
            ordered_and_filtered_requested_models.append(model_name)
            logger.warning(f"Model '{model_name}' was requested but is not in the predefined ordered list. Appending for GPU assignment.")
    
    if not ordered_and_filtered_requested_models:
        logger.warning("No models to assign GPUs to after filtering.")
        return {}

    assignments = {}
    gpu_idx = 0
    for model_name in ordered_and_filtered_requested_models:
        assignments[model_name] = available_gpus_str[gpu_idx % num_gpus]
        gpu_idx += 1
        
    logger.info(f"GPU assignments: {assignments}")
    return assignments

def set_gpu_visibility(gpu_id: any) -> None:
    """
    Set CUDA_VISIBLE_DEVICES environment variable to restrict GPU visibility.
    
    Args:
        gpu_id: ID of the GPU to make visible
    """
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id) 