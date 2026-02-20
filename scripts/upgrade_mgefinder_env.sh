#!/bin/bash
# Upgrade mgefinder_env to include all needed packages for the pipeline

echo ">>> Upgrading mgefinder_env with pipeline dependencies..."

# Activate mgefinder_env
eval "$(micromamba shell hook --shell bash)"
micromamba activate mgefinder_env

echo "Current environment: $CONDA_DEFAULT_ENV"
echo "Python version: $(python --version)"

echo ">>> Installing snakemake and dependencies..."
micromamba install -c conda-forge -c bioconda snakemake pulp pyyaml pandas -y

echo ">>> Installing fastq-dl for downloads..."
micromamba install -c conda-forge -c bioconda fastq-dl -y

echo ">>> Checking installations..."
echo "Snakemake: $(snakemake --version 2>/dev/null || echo 'FAILED')"
echo "MGEfinder: $(mgefinder --version 2>/dev/null || echo 'FAILED')"
echo "fastq-dl: $(fastq-dl --version 2>/dev/null || echo 'FAILED')"

echo ">>> Testing pulp..."
python -c "import pulp; print(f'pulp {pulp.__version__} - listSolvers available: {hasattr(pulp, \"listSolvers\")}')" 2>/dev/null || echo "pulp test failed"

echo ">>> Upgrade complete!"