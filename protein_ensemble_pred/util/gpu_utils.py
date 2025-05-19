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

def assign_gpus_to_models(num_models: int, force_sequential: bool = False) -> Dict[str, int]:
    """
    Assign GPUs to models based on availability and sequential flag.
    
    Args:
        num_models: Number of models that need GPU assignment
        force_sequential: If True, all models are assigned to the first available GPU.
        
    Returns:
        Dictionary mapping model names to GPU IDs
    """
    available_gpus = detect_available_gpus()
    if not available_gpus:
        logger.error("No GPUs detected. Cannot assign GPUs to models.")
        return {}
    
    model_keys = ["alphafold3", "boltz1"][:num_models] 

    if force_sequential or len(available_gpus) < num_models:
        if len(available_gpus) < num_models and not force_sequential:
            logger.warning(
                f"Only {len(available_gpus)} GPU(s) available for {num_models} models. "
                "Will run models sequentially on the first available GPU."
            )
        elif force_sequential:
            logger.info(f"Forcing sequential execution on GPU {available_gpus[0]} for {num_models} models.")
        
        assignments = {}
        for i in range(num_models):
            if i < len(model_keys):
                 assignments[model_keys[i]] = available_gpus[0]
            else: 
                 assignments[f"model_{i+1}"] = available_gpus[0] 
        return assignments
    
    assignments = {}
    for i in range(num_models):
        if i < len(model_keys):
            assignments[model_keys[i]] = available_gpus[i]
        else:
            assignments[f"model_{i+1}"] = available_gpus[i]
    logger.info(f"Assigned GPUs: {assignments}")
    return assignments

def set_gpu_visibility(gpu_id: int) -> None:
    """
    Set CUDA_VISIBLE_DEVICES environment variable to restrict GPU visibility.
    
    Args:
        gpu_id: ID of the GPU to make visible
    """
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id) 