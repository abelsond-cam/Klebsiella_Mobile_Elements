#!/bin/bash
#SBATCH --job-name=mge_pool
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=16
#SBATCH --partition=icelake-himem
#SBATCH --account=FLOTO-PROJECT-K-SL2-CPU
#SBATCH --output=mge_%x_%j.log
#SBATCH --error=mge_%x_%j.err

# Shared per-SL read-pool PRODUCER. $1 = Sublineage (e.g. SL39).
# Downloads the SL's usable reads ONCE (chunked, parallel, retried) into a pool
# shared by every reference-level consumer. fastq-dl env activated ONCE here
# (no per-sample micromamba -> no mamba-lock storm). Compute node has ENA
# access. Idempotent: requeue resumes via ready/ + dropped_samples.tsv.
# Submit: sbatch --job-name=mge_pool_SL39 scripts/submit_readpool_producer.sh SL39

set -euo pipefail
SL="${1:?Usage: sbatch submit_readpool_producer.sh <Sublineage>}"
POOL="/home/dca36/rds/rds-floto-bacterial-4k08a2yyQLw/david/processed/mgefinder/_readpool/${SL}"

echo "=== read-pool producer: ${SL} ==="
echo "Job ${SLURM_JOB_ID:-NA}  Node $(hostname)  Pool ${POOL}  Started $(date)"

cd "$SLURM_SUBMIT_DIR"
eval "$(micromamba shell hook --shell bash)"
micromamba activate fastq-dl

python3 src/run_readpool_producer.py \
    --config config/config.yaml \
    --sublineage "$SL" \
    --pool-dir "$POOL" \
    --chunk-size 100 \
    --workers "${SLURM_CPUS_PER_TASK:-16}" \
    --max-retries 4 \
    --backoff-base 20

echo ">>> producer ${SL} finished: $(date)"
