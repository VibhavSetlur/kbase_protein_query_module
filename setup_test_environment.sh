#!/bin/bash
# Setup script for KBase Protein Query Module test environment
# This script sets up the conda environment with proper dependencies

set -e  # Exit on any error

echo "Setting up KBase Protein Query Module test environment..."

# Source conda
source /opt/anaconda3/etc/profile.d/conda.sh

# Remove existing environment if it exists
echo "Removing existing environment if it exists..."
conda env remove -n kbase_protein_query_module_env -y 2>/dev/null || true

# Create new environment with Python 3.8
echo "Creating conda environment with Python 3.8..."
conda create -y -n kbase_protein_query_module_env python=3.8

# Install essential packages first
echo "Installing essential packages..."
conda run -n kbase_protein_query_module_env pip install --user -r requirements-essential.txt

# Install torch CPU version (to avoid CUDA dependencies)
echo "Installing PyTorch CPU version..."
conda run -n kbase_protein_query_module_env pip install --user torch --index-url https://download.pytorch.org/whl/cpu

# Install additional ML packages
echo "Installing ML packages..."
conda run -n kbase_protein_query_module_env pip install --user transformers networkx scikit-learn --no-deps

# Verify installation
echo "Verifying installation..."
conda run -n kbase_protein_query_module_env python -c "
import torch
import numpy as np
import pandas as pd
import networkx as nx
import pytest
print('✓ All essential packages installed successfully')
print(f'✓ Python version: {torch.version.__version__ if hasattr(torch, \"version\") else \"3.8\"}')
print(f'✓ PyTorch version: {torch.__version__}')
print(f'✓ CUDA available: {torch.cuda.is_available()}')
"

echo ""
echo "✅ Test environment setup complete!"
echo ""
echo "To activate the environment:"
echo "  conda activate kbase_protein_query_module_env"
echo ""
echo "To run tests:"
echo "  conda run -n kbase_protein_query_module_env python -m pytest test/ -v"
echo ""
echo "Note: Some packages like faiss-cpu, plotly, and matplotlib are commented out"
echo "in requirements.txt due to disk space constraints on this server."

