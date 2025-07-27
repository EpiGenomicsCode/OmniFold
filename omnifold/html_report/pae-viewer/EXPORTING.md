# Creating a Standalone HTML Report

This document explains how to generate a single, self-contained HTML file that bundles all necessary CSS, JavaScript, images, and data for the PAE Viewer. This file can be shared and opened on any modern browser without requiring a running web server.

## Setup

To use the export script, you need:
1.  **Python 3**: Version 3.9 or higher.
2.  **Node.js and npm**: Required to install the JavaScript bundler. You can download them from [nodejs.org](https://nodejs.org/).

Once Python and Node.js are installed, open a terminal in the project's root directory and run the following command to install the required `esbuild` dependency:

```bash
npm install
```

## Usage

After the setup is complete, you can generate a standalone report by running the `export_report.py` script. You must provide paths to the structure file, scores file, and chain labels.

For example:
```bash
python3 resources/scripts/export_report.py \\
  --structure resources/sample/GatA-GatB/fold_gata_gatb_model_0.cif \\
  --scores resources/sample/GatA-GatB/fold_gata_gatb_full_data_0.json \\
  --crosslinks resources/sample/GatA-GatB/GatA-GatB.csv \\
  --labels "GatA;GatB" \\
  --output GatA-GatB --open
```

Alternatively, you can run the entire command on a single line to avoid any shell interpretation issues:
```bash
python3 resources/scripts/export_report.py --structure resources/sample/GatA-GatB/fold_gata_gatb_model_0.cif --scores resources/sample/GatA-GatB/fold_gata_gatb_full_data_0.json --crosslinks resources/sample/GatA-GatB/GatA-GatB.csv --labels "GatA;GatB" --output GatA-GatB --open
```

- The script will create a file named `GatA-GatB_report_[timestamp].html` in the root directory.
- The `--open` flag automatically opens the generated report in your default web browser.

## Processing "Chai" Model Outputs

The script also includes a special mode to directly process the output directories from the "Chai" model.

### Chai Usage

When using this mode, you only need to provide the path to the output directory and the model index you wish to visualize. The script will automatically find the correct structure, PAE, and scores files.

```bash
python3 resources/scripts/export_report.py \\
  --chai-input-dir /path/to/your/chai1/ \\
  --model-index 0 \\
  --open
```

- `--chai-input-dir`: The path to the directory containing the model outputs (e.g., `scores.model_idx_0.npz`, `pae.model_idx_0.npz`, etc.).
- `--model-index`: The integer index of the model you want to process (e.g., 0, 1, 2, ...). Defaults to 0.
- When using `--chai-input-dir`, you do not need to provide `--structure`, `--scores`, `--labels`, or `--crosslinks`.

## Processing "Boltz" Model Outputs

The script also supports processing output directories from "Boltz" models.

### Boltz Usage

This mode requires the path to the `predictions` directory and the specific model name prefix used in the output files.

```bash
python3 resources/scripts/export_report.py \\
  --boltz-input-dir /path/to/your/predictions/ \\
  --boltz-model-name "Hemoglobin_tetramer_boltz_inference_generated_with_msas" \\
  --model-index 0 \\
  --open
```

- `--boltz-input-dir`: The path to the directory containing the model-specific sub-directory (e.g., `.../predictions/`).
- `--boltz-model-name`: The unique name of the model run, which is part of every filename.
- `--model-index`: The integer index of the model you want to process (e.g., 0, 1, 2, ...). Defaults to 0.
- When using `--boltz-input-dir`, you do not need to provide the standard data arguments. 