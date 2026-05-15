#!/bin/bash
#SBATCH --job-name=mge_sl39_full
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=76
#SBATCH --partition=icelake-himem
#SBATCH --account=FLOTO-PROJECT-K-SL2-CPU
#SBATCH --output=mge_sl39_full_%j.log
#SBATCH --error=mge_sl39_full_%j.err

# SL39 FULL run (~900 samples), metadata-driven, stream_fastq.
# Only submit after submit_sl39_smoketest.sh has passed (proves compute-node
# internet + the streaming mechanism).
#   --test-n -1  -> --limit -1 -> all SL39 samples (no cap)
#   --stream-fastq -> per-sample download, temp()-deleted right after bwa,
#                     so peak FASTQ disk ~= --jobs samples in flight (not ~750 GB)
#   --output-dir -> isolated SL39 tree (won't touch the old trial working dir)
# bwa is single-threaded, so parallelism = --jobs = concurrent samples.
# Submit: sbatch scripts/submit_sl39_full.sh

set -euo pipefail
echo "=== MGEfinder SL39 FULL run (all samples) ==="
echo "Job ID: ${SLURM_JOB_ID:-NA}  Node: $(hostname)  Started: $(date)"
echo "CPUs/jobs: ${SLURM_CPUS_PER_TASK}  Walltime: 36h"

cd "$SLURM_SUBMIT_DIR"
eval "$(micromamba shell hook --shell bash)"
micromamba activate mgefinder_env

python3 src/run_pipeline.py \
    --config config/config.yaml \
    --output-dir /home/dca36/rds/rds-floto-bacterial-4k08a2yyQLw/david/processed/mgefinder/SL39 \
    --stream-fastq \
    --test-n -1 \
    --verbose \
    --jobs "$SLURM_CPUS_PER_TASK"

echo ">>> Full SL39 run completed: $(date)"
