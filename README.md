# OmniFold

## Overview

This command-line application simplifies running ensemble protein structure predictions using **AlphaFold3**, **Boltz-2**, and **Chai-1** on High-Performance Computing (HPC) clusters. 

**The key advantage:** You provide your target sequence(s) in **one** of the supported input formats (FASTA, AlphaFold3 JSON, or Boltz-2 YAML), specify how Multiple Sequence Alignments (MSAs) should be obtained (or let the tool generate them automatically), and the application handles the rest:

*   **Internal Conversion:** Automatically converts your input into the specific formats required by each model.
*   **Unified MSA Handling:** Manages MSA generation or reuse consistently.
    *   If the AlphaFold3 MSA pipeline (Hmmer) is used, the unpaired and paired A3M files are extracted to form the csv file for Boltz to consume.
    *   Additionally, if Chai-1 is being run, the AlphaFold3 internal a3m states are extracted to suit Chai-1's PQT format.
    *   When using the `colabfold` MSA method, the API is queried only once. The resulting unpaired and paired MSAs are passed to the AlphaFold3 pipeline, while Boltz and Chai-1 use their native ColabFold integrations with the same cached MSA files.
*   **Parallel Execution:** Orchestrates predictions with AlphaFold3, Boltz-2, and Chai-1, potentially in parallel on different GPUs.
*   **Containerized Runs:** Executes models reliably within their Singularity containers.
*   **Organized Output:** Saves the native outputs from each model into a structured output directory.
*   **Automated Reporting:** Automatically generates a comprehensive, interactive HTML report (`OmniFold_Report.zip`) comparing all model outputs, including metrics like `ipSAE` and `pDockQ`.

This eliminates the need for manual format conversions and separate pipeline runs for each model, streamlining your ensemble prediction workflow.

## Requirements

### Software

*   **Singularity (or Apptainer):** Must be installed on the HPC system to run the containerized models.
*   **Python 3.9+:** Required for the CLI application itself.
*   **Dependencies:**  A full list can be found in `requirements.txt`.

### Input Formats

This tool accepts inputs in FASTA, AlphaFold3 JSON, and Boltz-2 YAML formats. 

*   For detailed **FASTA** formatting guidelines, see [docs/fasta.md](docs/fasta.md).
*   For the official **AlphaFold3 JSON** input specification, please refer to the [AlphaFold3 Input Documentation](https://github.com/google-deepmind/alphafold3/blob/main/docs/input.md).
*   For the official **Boltz-2 YAML** input specification, please refer to the [Boltz Prediction Documentation](https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md).
*   For **Chai-1 FASTA** input format, the tool generates it automatically if using AlphaFold3 MSAs. The headers follow the pattern `>protein|name=chain_X`, `>rna|name=chain_Y`, etc. A standard fasta format sticking to the [documentation](docs/fasta.md) is recommended but not required. 

### Model and Data Paths

The following paths must be provided as command-line arguments if the respective model is intended to be run:

*   **AlphaFold3 Singularity Image (`--alphafold3_sif_path`):** Absolute path to the AlphaFold3 Singularity image file (`.sif`). If not provided, AlphaFold3 will be skipped.
*   **Boltz-2 Singularity Image (`--boltz1_sif_path`):** Absolute path to the Boltz Singularity image file (`.sif`). If not provided, Boltz will be skipped.
*   **Chai-1 Singularity Image (`--chai1_sif_path`):** Absolute path to the Chai-1 Singularity image file (`.sif`). If not provided, Chai-1 will be skipped.
*   **AlphaFold3 Model Weights (`--alphafold3_model_weights_dir`):** Absolute path to the directory containing the downloaded AlphaFold3 model parameters/weights. (Required if running AlphaFold3).

### Optional (but often necessary) Paths

*   **AlphaFold3 Databases (`--alphafold3_database_dir`):** Absolute path to the root directory containing AlphaFold3 databases (e.g., UniRef, MGnify, PDB, etc.). This is required if the application needs to run the AlphaFold3 MSA generation pipeline.
*   **ColabFold MSA Server URL (`--colabfold_msa_server_url`):** URL for a ColabFold MMseqs2 API server. If using the `colabfold` MSA method, and you wish to use a specific (e.g., local) server, provide its URL. If not provided when `msa_method` is `colabfold`, Boltz will use its internal default server URL.

## Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/EpiGenomicsCode/OmniFold.git
    cd OmniFold
    ```
2.  **Install Python Dependencies:** It is highly recommended to use a virtual environment.
    ```bash
    python3 -m venv .venv
    source .venv/bin/activate
    pip install .
    ```
3.  **Install PAE Viewer Dependencies:** To enable the automatic generation of the interactive PAE viewers, you must also install the Node.js dependencies. This is a one-time setup:
    ```bash
    cd omnifold/html_report/pae-viewer
    npm install
    cd ../../..  # Return to the project root directory
    ```
4.  **Ensure Singularity Images and Data are Accessible:**
    *   Download or build the AlphaFold3, Boltz-2, and Chai-1 Singularity images. **We recommend using the pre-built images available from [Protein Structure Containers](https://github.com/EpiGenomicsCode/ProteinStruct-Containers).** Place the `.sif` files in accessible locations on your HPC system.
    *   **Download the AlphaFold3 Model Weights:** Access to the official AlphaFold3 model parameters requires registration for non-commercial use via the [AlphaFold 3 Model Parameters Request Form](https://docs.google.com/forms/d/e/1FAIpQLSfWZAgo1aYk0O4MuAXZj8xRQ8DafeFJnldNOnh_13qAx2ceZw/viewform). Ensure you meet the terms and download the weights to an accessible directory path.
    *   **Download AlphaFold3 Databases (if needed):** If you plan to use the `alphafold3` MSA generation method, you must download the required databases. Use the official [fetch_databases.sh script](https://github.com/google-deepmind/alphafold3/blob/main/fetch_databases.sh) provided by Google DeepMind. Ensure the databases are stored in an accessible directory path.

## Basic Usage

To run a prediction, use the `omnifold` command followed by your desired arguments:

```bash
omnifold \
    --input_file /path/to/your/input.fasta \
    --output_dir /path/to/your/output_directory \
    --alphafold3_sif_path /path/to/alphafold3.sif \
    --boltz1_sif_path /path/to/boltz1.sif \
    --chai1_sif_path /path/to/chai1.sif \
    --alphafold3_model_weights_dir /path/to/af3_weights \
    --alphafold3_database_dir /path/to/af3_databases \
    --msa_method alphafold3 \
    --log_level INFO
    # Chai-1 specific parameters (examples, actual flags might differ):
    # --chai1_num_recycling_steps 6 \
    # --chai1_use_msa_server \
```

### Example Command:

The following command runs a prediction for `tests/files/multi_chain_1.fasta`, attempting to run all three models (AlphaFold3, Boltz-1, Chai-1), saving results to `results/multi_chain_1_output`, and specifies paths for Singularity images, model weights, and databases. It uses the `alphafold3` method for MSA generation and sets the log level to `INFO`.

```bash
omnifold --input_file tests/files/multi_chain_1.fasta \
    --output_dir results/multi_chain_1_output \
    --alphafold3_sif_path /../containers/alphafold3_x86.sif \
    --boltz1_sif_path /../containers/boltz_x86.sif \
    --chai1_sif_path /../containers/chai_lab_x86.sif \
    --alphafold3_model_weights_dir /../alphafold3_weights \
    --alphafold3_database_dir /../databases \
    --msa_method alphafold3 \
    --log_level INFO
    # Add any desired Chai-1 specific arguments here
```

Refer to the command-line help for a full list of options and their descriptions:

```bash
omnifold --help
```

For detailed information on FASTA input formatting, see [docs/fasta.md](docs/fasta.md).

## Output Structure

The application will create the specified output directory. Inside this directory, you will typically find:

*   Subdirectories for AlphaFold3, Boltz, and Chai-1 containing their respective native output files (PDB/CIF structures, confidence scores, etc.).
*   Configuration files generated for each model.
*   Log files (`ensemble_prediction.log`, `alphafold3_run.log`, `boltz_run.log`, `chai1_run.log`).
*   If MSAs were generated, intermediate MSA files may also be present in a subdirectory (e.g., `msa_intermediate_files`, `msas_forChai`).

## HTML Report Generation

At the end of a successful pipeline run, the tool automatically generates a comprehensive `OmniFold_Report.zip` file in your output directory. This ZIP file is shareable and contains:
1.  `final_report.html`: A detailed comparison of all model predictions.
2.  `pae_viewers/`: A directory containing standalone HTML files, one for each model's prediction. Each file opens an interactive viewer that couples the 3D structure with its corresponding PAE (Predicted Aligned Error) matrix. Users can select regions on the PAE plot to highlight them on the 3D model and vice-versa.

The report includes:
-   A summary table of key metrics (pLDDT, pTM, ipTM) for the best model from each method.
-   An interface confidence table with `ipSAE` and `pDockQ` scores for all chain pairs. The `ipSAE` score is calculated based on the method described by Dunbrack et al.
-   An interactive, overlaid pLDDT plot to compare per-residue confidence across models.

### Chai-1 PAE & ipSAE Support
To enable full interface analysis, this tool uses modified Chai-1 scripts (located in `omnifold/chai1_modifications`) that are bound into the container at runtime. This ensures that Chai-1 outputs the PAE matrix required for `ipSAE` and `pDockQ` calculations.

## How it Works

1.  **Unified Input Handling:** Parses your single input file (FASTA, AF3 JSON, or Boltz YAML) and standardizes the sequence and chain information internally.
2.  **MSA Management:** Determines if MSAs are needed based on your input and `--msa_method` flag.
    *   If `msa_method` is `alphafold3` (default), it runs the AlphaFold3 data pipeline using its Singularity container. This tool utilizes modified versions of AlphaFold3's internal `pipeline.py` and `run_alphafold.py` scripts (bound from `omnifold/singularity_af3/...` into the container at runtime) to ensure comprehensive A3M file generation (e.g., UniRef90, MGnify, etc.) for each chain within the standard AlphaFold3 output structure. Caching is restructured to allow for this change.
        *   The resulting AlphaFold3 `_data.json` is parsed to extract per-protein A3M files from the generated MSAs, which are then made available for Boltz.
        *   If Chai-1 is to be run, the A3M files from the AlphaFold3 output (`msas/chain_X/*.a3m`) are converted into Chai-1's PQT format (`msas_forChai/*.pqt`).
    *   If `msa_method` is `colabfold`, the tool queries the ColabFold API once to retrieve MSAs. The resulting A3M files are supplied to the AlphaFold3 pipeline, while Boltz and Chai-1 consume the same cached MSAs using their native integrations.
    *   Existing MSAs from the input file can also be used, bypassing generation.
3.  **Configuration Generation:** Creates the specific input files (AF3 JSON, Boltz YAML, Chai-1 FASTA) required by each model, incorporating the standardized sequence data and consistent MSA information.
4.  **Orchestration & Execution:**
    *   Determines which models (AlphaFold3, Boltz, Chai-1) to run based on whether their respective Singularity image paths (`--alphafold3_sif_path`, `--boltz1_sif_path`, `--chai1_sif_path`) have been provided by the user. Models without a SIF path are skipped.
    *   Detects available GPUs.
    *   Assigns GPUs to the selected models (one per model if available and different models are run, or runs sequentially on a single GPU).
    *   Constructs and executes `singularity run/exec` commands for the selected models, binding necessary directories (input configs, output, model weights, databases, MSAs).
5.  **Output Collection:** Gathers results and logs from all executed model runs into the specified output directory.
6.  **Report Generation:** Finally, the tool generates a comprehensive `OmniFold_Report.zip` file containing a detailed HTML report that compares key metrics (`pLDDT`, `ipSAE`, `pDockQ`) and includes interactive, shareable PAE viewers that link the 3D structure to the PAE matrix for in-depth analysis.

## Acknowledgements
- AlphaFold by DeepMind Technologies Limited
- Boltz-2 by Passaro, Saro, et al. "Boltz-2: Towards Accurate and Efficient Binding Affinity Prediction" bioRxiv (2025)
- Chai 1 by Chai Discovery, Inc.
- The interactive PAE viewer included in the HTML report is adapted from the original [PAE Viewer](https://gitlab.gwdg.de/general-microbiology/pae-viewer) developed by the Department of General Microbiology (Institute of Microbiology and Genetics, Georg August University of Göttingen) under the direction of Jörg Stülke.
- The research project is generously funded by Cornell University BRC Epigenomics Core Facility (RRID:SCR_021287), Penn State Institute for Computational and Data Sciences (RRID:SCR_025154) , Penn State University Center for Applications of Artificial Intelligence and Machine Learning to Industry Core Facility (AIMI) (RRID:SCR_022867) and supported by a gift to AIMI research from Dell Technologies.
- Computational support was provided by NSF ACCESS to William KM Lai and Gretta Kellogg through BIO230041

## References
- Dunbrack RL Jr. Rēs ipSAE loquunt: What's wrong with AlphaFold's ipTM score and how to fix it. bioRxiv [Preprint]. 2025 Feb 14:2025.02.10.637595. doi: 10.1101/2025.02.10.637595. PMID: 39990437; PMCID: PMC11844409. [Available at: https://pmc.ncbi.nlm.nih.gov/articles/PMC11844409/](https://pmc.ncbi.nlm.nih.gov/articles/PMC11844409/)
