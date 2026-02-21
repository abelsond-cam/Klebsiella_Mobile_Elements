#!/bin/bash
#SBATCH --job-name=mgefinder
#SBATCH --time=04:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --partition=icelake
#SBATCH --account=FLOTO-SL2-CPU
#SBATCH --output=mgefinder_trial_%j.log
#SBATCH --error=mgefinder_trial_%j.err

# HPC job script for MGEfinder pipeline (runs rows 0, 1, 2 with --skip-download).
# Submit: sbatch scripts/submit_hpc.sh
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

# Run pipeline for all rows
for row in 0 1 2; do
    echo ">>> Running pipeline for row $row..."
    python3 src/run_pipeline.py --config config/config.yaml --verbose --jobs $SLURM_CPUS_PER_TASK --row $row --skip-download
    if [ $? -ne 0 ]; then
        echo "ERROR: Pipeline failed for row $row"
        exit 1
    fi
done

echo ">>> Job completed: $(date)"