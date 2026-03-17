#!/bin/bash

# SETUP.sh - Environment and Permissions Setup for PolySNiP-nf

echo "=========================================="
echo "   PolySNiP-nf Setup Script"
echo "=========================================="

# 1. Ensure Python scripts in bin/ are executable
echo "[1/3] Setting permissions for scripts in bin/..."
chmod +x bin/*.py
if [ $? -eq 0 ]; then
    echo "Success: bin/ scripts are now executable."
else
    echo "Error: Failed to set permissions."
    exit 1
fi

# 2. Create Conda environment
echo "[2/3] Creating conda environment 'polysnip-nf' from environment.yml..."
if command -v conda &> /dev/null; then
    conda env create -f environment.yml
    if [ $? -eq 0 ]; then
        echo "Success: Conda environment 'polysnip-nf' created."
    else
        echo "Note: If environment already exists, you can update it with 'conda env update -f environment.yml'."
    fi
else
    echo "Warning: 'conda' command not found. Please ensure Conda is installed and in your PATH."
    echo "Once installed, run: conda env create -f environment.yml"
fi

# 3. Create expected directories
echo "[3/3] Preparing data and refs directories..."
mkdir -p data refs

echo "=========================================="
echo "   Setup Complete!"
echo "=========================================="
echo "To activate the environment manually:"
echo "  conda activate polysnip-nf"
echo ""
echo "To run the pipeline via Slurm:"
echo "  sbatch submit.sh"
echo "=========================================="
