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

def assign_gpus_to_models(num_models: int) -> Dict[str, int]:
    """
    Assign GPUs to models based on availability.
    
    Args:
        num_models: Number of models that need GPU assignment
        
    Returns:
        Dictionary mapping model names to GPU IDs
        
    Raises:
        RuntimeError: If no GPUs are detected
    """
    available_gpus = detect_available_gpus()
    if not available_gpus:
        raise RuntimeError("No GPUs detected. At least one GPU is required.")
    
    if len(available_gpus) < num_models:
        logger.warning(
            f"Only {len(available_gpus)} GPU(s) available for {num_models} models. "
            "Will run models sequentially."
        )
        # Return same GPU for all models - they'll run sequentially
        return {model: available_gpus[0] for model in ["alphafold3", "boltz1"]}
    
    # Assign different GPUs to each model
    return {
        "alphafold3": available_gpus[0],
        "boltz1": available_gpus[1]
    }

def set_gpu_visibility(gpu_id: int) -> None:
    """
    Set CUDA_VISIBLE_DEVICES environment variable to restrict GPU visibility.
    
    Args:
        gpu_id: ID of the GPU to make visible
    """
    os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id) 