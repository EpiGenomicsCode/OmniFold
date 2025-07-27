#!/bin/bash
# Exit immediately if a command exits with a non-zero status.
set -e

# --- Configuration ---
ENV_NAME="omnifold_env"
ENV_YML="environment.yml"

# --- Helper Functions ---
print_info() {
    echo -e "\\033[1;34m[INFO]\\033[0m $1"
}

print_success() {
    echo -e "\\033[1;32m[SUCCESS]\\033[0m $1"
}

print_warning() {
    echo -e "\\033[1;33m[WARNING]\\033[0m $1"
}

command_exists() {
    command -v "$1" &> /dev/null
}

# --- Main Logic ---
print_info "Starting OmniFold Setup..."

# 1. Check for Conda
if ! command_exists conda; then
    print_warning "Conda is not installed or not in your PATH."
    echo "Please install Miniconda or Anaconda first by following the instructions at:"
    echo "https://docs.conda.io/projects/conda/en/latest/user-guide/install/"
    exit 1
fi

# 2. Create or Update Conda Environment
if conda env list | grep -q "^${ENV_NAME}\\s"; then
    print_info "Conda environment '${ENV_NAME}' already exists. Updating..."
    conda env update --file "${ENV_YML}" --name "${ENV_NAME}" --prune
else
    print_info "Creating new Conda environment '${ENV_NAME}' from ${ENV_YML}..."
    conda env create --file "${ENV_YML}"
fi

# 3. Activate the environment to run subsequent commands
print_info "Activating the Conda environment..."
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "${ENV_NAME}"

# 4. Install PAE Viewer dependencies
PAE_VIEWER_DIR="omnifold/html_report/pae-viewer"
if [ -d "$PAE_VIEWER_DIR" ]; then
    print_info "Installing Node.js dependencies for the PAE Viewer..."
    cd "$PAE_VIEWER_DIR"
    npm install
    cd - > /dev/null # Return to previous directory quietly
else
    print_warning "PAE Viewer directory not found at ${PAE_VIEWER_DIR}. Skipping npm install."
fi

# 5. Install the OmniFold package itself using pip
print_info "Installing the OmniFold package in editable mode..."
pip install -e .

print_success "OmniFold setup is complete!"
print_info "To use OmniFold, first activate the environment with:"
print_info "conda activate ${ENV_NAME}" 