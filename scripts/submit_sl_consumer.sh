#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=76
#SBATCH --partition=icelake-himem
#SBATCH --account=FLOTO-PROJECT-K-SL2-CPU
#SBATCH --output=mge_%x_%j.log
#SBATCH --error=mge_%x_%j.err

# Per-chunk CONSUMER. $1=Sublineage $2=ref level {f|d} $3=chunk index k.
# NON-STREAM, reads from the shared pool (producer-owned). Runs ONLY this
# chunk's explicit per-sample 03.inferseq targets (no cohort aggregation).
# Waits (bounded) for the producer's chunk_<k>.READY before starting Snakemake.
# Submit (chained): sbatch --job-name=mge_SL39_f_c0 --dependency=after:<P0> \
#   scripts/submit_sl_consumer.sh SL39 f 0

set -euo pipefail
SL="${1:?need Sublineage}"; LVL="${2:?need f|d}"; K="${3:?need chunk index}"
case "$LVL" in f|d) ;; *) echo "level must be f|d" >&2; exit 1;; esac
POOL="/home/dca36/rds/rds-floto-bacterial-4k08a2yyQLw/david/processed/mgefinder/_readpool/${SL}"
WD="/home/dca36/rds/rds-floto-bacterial-4k08a2yyQLw/david/processed/mgefinder/${SL}/ref_${LVL}"
READY="${POOL}/chunks/chunk_${K}.READY"

echo "=== consumer ${SL} ref_${LVL} chunk ${K} ==="
echo "Job ${SLURM_JOB_ID:-NA}  Node $(hostname)  Started $(date)"

# Bounded wait for the producer to have downloaded this chunk (it runs ahead).
deadline=$(( $(date +%s) + 30*3600 ))
while [ ! -f "$READY" ]; do
    if [ "$(date +%s)" -ge "$deadline" ]; then
        echo "ERROR: chunk_${K}.READY not present after 30h wait: $READY" >&2
        exit 1
    fi
    echo "waiting for $READY ... $(date)"; sleep 120
done
echo ">>> chunk ${K} READY; starting compute $(date)"

cd "$SLURM_SUBMIT_DIR"
eval "$(micromamba shell hook --shell bash)"
micromamba activate mgefinder_env

python3 src/run_pipeline.py \
    --config config/config.yaml \
    --sublineage "$SL" \
    --reference-level "$LVL" \
    --output-dir "$WD" \
    --read-pool "$POOL" \
    --consumer-tag "ref_${LVL}" \
    --chunk "$K" \
    --drop-file "${POOL}/dropped_samples.tsv" \
    --test-n -1 \
    --verbose \
    --jobs "$SLURM_CPUS_PER_TASK"

echo ">>> consumer ${SL} ref_${LVL} chunk ${K} done: $(date)"
