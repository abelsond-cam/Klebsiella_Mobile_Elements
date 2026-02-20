#!/bin/bash
# Fix Snakemake compatibility issues in mgefinder_env

echo ">>> Fixing Snakemake compatibility in mgefinder_env..."

# Activate mgefinder_env
eval "$(micromamba shell hook --shell bash)"
micromamba activate mgefinder_env

echo "Current environment: $CONDA_DEFAULT_ENV"

# The issue is that MGEfinder was built for an older Snakemake version
# Let's try to install a compatible Snakemake version
echo ">>> Installing compatible Snakemake version..."

# Try Snakemake 7.x which is more compatible with older tools
micromamba install -c conda-forge -c bioconda "snakemake>=7.0,<8.0" -y

echo ">>> Testing after downgrade..."
python -c "import snakemake; print('Snakemake version:', snakemake.__version__)" 2>/dev/null || echo "Snakemake import still fails"
mgefinder --version 2>/dev/null || echo "MGEfinder still has issues"

echo ">>> If still failing, will try without Snakemake import in MGEfinder..."