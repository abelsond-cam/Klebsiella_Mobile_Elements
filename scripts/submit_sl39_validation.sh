#!/bin/bash
#SBATCH --job-name=mge_sl39_val
#SBATCH --time=01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=icelake-himem
#SBATCH --account=FLOTO-PROJECT-K-SL2-CPU
#SBATCH --output=mge_sl39_val_%j.log
#SBATCH --error=mge_sl39_val_%j.err

# SL39 one-sample validation run (metadata-driven selection).
# Sample/reference come from config/config.yaml: sublineage=SL39, validate_n=1,
# reference_sample_name=GCF_000814305.1_ASM81430v1_genomic.
# FASTQ for the sample is downloaded beforehand, so --skip-download.
# run_pipeline.py activates the correct micromamba env per phase itself.
# Submit: sbatch scripts/submit_sl39_validation.sh

set -euo pipefail
echo "=== MGEfinder SL39 validation ==="
echo "Job ID: ${SLURM_JOB_ID:-NA}  Node: $(hostname)  Started: $(date)"

cd "$SLURM_SUBMIT_DIR"
eval "$(micromamba shell hook --shell bash)"
micromamba activate mgefinder_env

python3 src/run_pipeline.py \
    --config config/config.yaml \
    --verbose \
    --jobs "$SLURM_CPUS_PER_TASK" \
    --skip-download

echo ">>> Job completed: $(date)"
