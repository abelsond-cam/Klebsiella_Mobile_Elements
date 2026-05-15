#!/bin/bash
#SBATCH --job-name=mge_sl39_smoke
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --partition=icelake-himem
#SBATCH --account=FLOTO-PROJECT-K-SL2-CPU
#SBATCH --output=mge_sl39_smoke_%j.log
#SBATCH --error=mge_sl39_smoke_%j.err

# SL39 stream_fastq SMOKE TEST (2 samples) before committing to the full 36h run.
# Validates the two unproven things in one go:
#   1. compute-node -> ENA outbound internet (fastq-dl inside the SLURM job)
#   2. the download_fastq -> bwa -> temp()-delete streaming mechanism end-to-end
# Writes to a THROWAWAY output dir so it never touches the real SL39 results.
# Submit: sbatch scripts/submit_sl39_smoketest.sh

set -euo pipefail
echo "=== MGEfinder SL39 stream_fastq smoke test (n=2) ==="
echo "Job ID: ${SLURM_JOB_ID:-NA}  Node: $(hostname)  Started: $(date)"

cd "$SLURM_SUBMIT_DIR"
eval "$(micromamba shell hook --shell bash)"
micromamba activate mgefinder_env

python3 src/run_pipeline.py \
    --config config/config.yaml \
    --output-dir /home/dca36/rds/rds-floto-bacterial-4k08a2yyQLw/david/processed/mgefinder/SL39_smoketest \
    --stream-fastq \
    --test-n 2 \
    --verbose \
    --jobs "$SLURM_CPUS_PER_TASK"

echo ">>> Smoke test completed: $(date)"
