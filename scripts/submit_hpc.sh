#!/bin/bash
#SBATCH --job-name=mgefinder
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --output=mgefinder_%j.log
#SBATCH --error=mgefinder_%j.err

# HPC job script for MGEfinder pipeline
# Adjust SBATCH parameters above for your system

echo "=== MGEfinder Pipeline HPC Job ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Started: $(date)"
echo "Node: $(hostname)"

# Change to project directory
cd $SLURM_SUBMIT_DIR

# Activate environment
eval "$(micromamba shell hook --shell bash)"
micromamba activate mgefinder_env

# Run pipeline
echo ">>> Running pipeline..."
python3 src/run_pipeline.py --config config/config.yaml --verbose --jobs $SLURM_CPUS_PER_TASK

echo ">>> Job completed: $(date)"