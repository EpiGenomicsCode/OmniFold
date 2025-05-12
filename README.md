# Protein Ensemble Prediction CLI

## Overview

This command-line application simplifies running ensemble protein structure predictions using both **AlphaFold3** and **Boltz-1** on High-Performance Computing (HPC) clusters. 

**The key advantage:** You provide your target sequence(s) in **one** of the supported input formats (FASTA, AlphaFold3 JSON, or Boltz-1 YAML), specify how Multiple Sequence Alignments (MSAs) should be obtained (or let the tool generate them automatically), and the application handles the rest:

*   **Internal Conversion:** Automatically converts your input into the specific formats required by both AlphaFold3 and Boltz-1.
*   **Unified MSA Handling:** Manages MSA generation or reuse consistently for both models based on your choice (using the AF3 pipeline or ColabFold method).
*   **Parallel Execution:** Orchestrates predictions with AlphaFold3 and Boltz-1, potentially in parallel on different GPUs.
*   **Containerized Runs:** Executes models reliably within their Singularity containers.
*   **Organized Output:** Saves the native outputs from each model into a structured output directory.

This eliminates the need for manual format conversions and separate pipeline runs for each model, streamlining your ensemble prediction workflow.

## Requirements

### Software

*   **Singularity (or Apptainer):** Must be installed on the HPC system to run the containerized models.
*   **Python 3.9+:** Required for the CLI application itself.
*   **Dependencies:**  A full list can be found in `requirements.txt`.

### Input Formats

This tool accepts inputs in FASTA, AlphaFold3 JSON, and Boltz-1 YAML formats. 

*   For detailed **FASTA** formatting guidelines, see [docs/fasta.md](docs/fasta.md).
*   For the official **AlphaFold3 JSON** input specification, please refer to the [AlphaFold3 Input Documentation](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md).
*   For the official **Boltz-1 YAML** input specification, please refer to the [Boltz Prediction Documentation](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md).

### Model and Data Paths

The following paths must be provided as command-line arguments when running the application:

*   **AlphaFold3 Singularity Image (`--alphafold3_sif_path`):** Absolute path to the AlphaFold3 Singularity image file (`.sif`). 
*   **Boltz-1 Singularity Image (`--boltz1_sif_path`):** Absolute path to the Boltz-1 Singularity image file (`.sif`).
*   **AlphaFold3 Model Weights (`--alphafold3_model_weights_dir`):** Absolute path to the directory containing the downloaded AlphaFold3 model parameters/weights.

### Optional (but often necessary) Paths

*   **AlphaFold3 Databases (`--alphafold3_database_dir`):** Absolute path to the root directory containing AlphaFold3 databases (e.g., UniRef, MGnify, PDB, etc.). This is required if the application needs to run the AlphaFold3 MSA generation pipeline.
*   **ColabFold MSA Server URL (`--colabfold_msa_server_url`):** URL for a ColabFold MMseqs2 API server. If using the `colabfold` MSA method, and you wish to use a specific (e.g., local) server, provide its URL. If not provided when `msa_method` is `colabfold`, Boltz will use its internal default server URL.

## Installation

1.  **Clone the repository:**
    ```bash
    git clone <repository_url>
    cd protein-ensemble-pred
    ```
2.  **Install Python dependencies:**
    ```bash
    pip install -r requirements.txt
    ```
3.  **Ensure Singularity Images and Data are Accessible:**
    *   Download or build the AlphaFold3 and Boltz-1 Singularity images. **We recommend using the pre-built images available from [https://github.com/EpiGenomicsCode/AF3-Container](https://github.com/EpiGenomicsCode/AF3-Container).** Place the `.sif` files in accessible locations on your HPC system.
    *   **Download the AlphaFold3 Model Weights:** Access to the official AlphaFold3 model parameters requires registration for non-commercial use via the [AlphaFold 3 Model Parameters Request Form](https://docs.google.com/forms/d/e/1FAIpQLSfWZAgo1aYk0O4MuAXZj8xRQ8DafeFJnldNOnh_13qAx2ceZw/viewform). Ensure you meet the terms and download the weights to an accessible directory path.
    *   **Download AlphaFold3 Databases (if needed):** If you plan to use the `alphafold3` MSA generation method, you must download the required databases. Use the official [fetch_databases.sh script](https://github.com/google-deepmind/alphafold3/blob/main/fetch_databases.sh) provided by Google DeepMind. Ensure the databases are stored in an accessible directory path.

## Basic Usage

To run a prediction:

```bash
python -m protein_ensemble_pred.cli \
    --input_file /path/to/your/input.fasta \
    --output_dir /path/to/your/output_directory \
    --alphafold3_sif_path /path/to/alphafold3.sif \
    --boltz1_sif_path /path/to/boltz1.sif \
    --alphafold3_model_weights_dir /path/to/af3_weights \
    --alphafold3_database_dir /path/to/af3_databases # (If AF3 MSA generation is needed)
    # Add other optional parameters as needed (e.g., --msa_method, model-specific params)
```

### Example Command:

```bash
python -m protein_ensemble_pred.cli \
    --input_file examples/T1084.fasta \
    --output_dir results/T1084_output \
    --alphafold3_sif_path /apps/containers/alphafold3.sif \
    --boltz1_sif_path /apps/containers/boltz1.sif \
    --alphafold3_model_weights_dir /data/alphafold3_params \
    --alphafold3_database_dir /data/alphafold_databases \
    --msa_method alphafold3 \
    --log_level INFO
```

Refer to the command-line help for a full list of options and their descriptions:

```bash
python -m protein_ensemble_pred.cli --help
```

For detailed information on FASTA input formatting, see [docs/fasta.md](docs/fasta.md).

## Output Structure

The application will create the specified output directory. Inside this directory, you will typically find:

*   Subdirectories for AlphaFold3 and Boltz-1 containing their respective native output files (PDB/CIF structures, confidence scores, etc.).
*   Configuration files generated for each model.
*   Log files (`ensemble_prediction.log`, `alphafold3_run.log`, `boltz_run.log`).
*   If MSAs were generated, intermediate MSA files may also be present in a subdirectory (e.g., `msa_intermediate_files`).

## How it Works

1.  **Unified Input Handling:** Parses your single input file (FASTA, AF3 JSON, or Boltz YAML) and standardizes the sequence and chain information internally.
2.  **MSA Management:** Determines if MSAs are needed based on your input and `--msa_method` flag. 
    *   If `msa_method` is `alphafold3` (default), it runs the AlphaFold3 data pipeline using its Singularity container to generate MSAs suitable for both models.
    *   If `msa_method` is `colabfold`, it invokes the Boltz Singularity container with flags to use its MSA server functionality (which typically calls a ColabFold API) to retrieve MSAs, again making them available for both models.
    *   Existing MSAs from the input file can also be used, bypassing generation.
3.  **Configuration Generation:** Creates the specific input files (AF3 JSON, Boltz YAML) required by each model, incorporating the standardized sequence data and consistent MSA information.
4.  **Orchestration & Execution:** 
    *   Detects available GPUs.
    *   Assigns GPUs to models (one per model if available) for parallel execution, or runs sequentially on a single GPU.
    *   Constructs and executes `singularity run/exec` commands for both AlphaFold3 and Boltz-1, binding necessary directories (input configs, output, model weights, databases).
5.  **Output Collection:** Gathers results and logs from both model runs into the specified output directory.
