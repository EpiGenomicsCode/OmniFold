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
    # First check CUDA_VISIBLE_DEVICES
    cuda_devices = os.environ.get("CUDA_VISIBLE_DEVICES")
    if cuda_devices:
        return [int(d) for d in cuda_devices.split(",")]
    
    # If not set, try to detect all GPUs
    try:
        import torch
        return list(range(torch.cuda.device_count()))
    except ImportError:
        # Fallback to nvidia-smi if torch not available
        try:
            result = subprocess.run(
                ["nvidia-smi", "-L"],
                capture_output=True,
                text=True,
                check=True
            )
            # Parse output to get GPU count
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
        # If running sequentially is allowed without GPU, this could return None for GPU IDs
        # For now, sticking to the original behavior of requiring GPUs if models are to be run.
        logger.error("No GPUs detected. Cannot assign GPUs to models.")
        # Return empty or raise, depending on how Orchestrator handles no GPU assignment
        return {}
    
    # Model names list - assuming fixed for now, but could be made dynamic
    # The number of models in this list should ideally match num_models if logic is complex.
    # For simple 2-model case, it's okay.
    model_keys = ["alphafold3", "boltz1"][:num_models] # Ensure we only assign for num_models

    if force_sequential or len(available_gpus) < num_models:
        if len(available_gpus) < num_models and not force_sequential:
            logger.warning(
                f"Only {len(available_gpus)} GPU(s) available for {num_models} models. "
                "Will run models sequentially on the first available GPU."
            )
        elif force_sequential:
            logger.info(f"Forcing sequential execution on GPU {available_gpus[0]} for {num_models} models.")
        
        # Assign the first available GPU to all requested models
        assignments = {}
        for i in range(num_models):
            if i < len(model_keys):
                 assignments[model_keys[i]] = available_gpus[0]
            else: # Should not happen if model_keys is sliced by num_models
                 assignments[f"model_{i+1}"] = available_gpus[0] 
        return assignments
    
    # Assign different GPUs to each model if enough are available and not forced sequential
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