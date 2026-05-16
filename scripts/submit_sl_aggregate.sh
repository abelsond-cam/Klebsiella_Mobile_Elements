#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=76
#SBATCH --partition=icelake-himem
#SBATCH --account=FLOTO-PROJECT-K-SL2-CPU
#SBATCH --output=mge_%x_%j.log
#SBATCH --error=mge_%x_%j.err

# Cohort AGGREGATION (one per reference level). $1=Sublineage $2=ref {f|d}.
# All per-sample TSVs already persisted by the chunk consumers. Regenerates the
# dataset to the SUCCESSFUL set (drop-file) so SAMPLES excludes permanently
# undownloadable accessions and expand(sample=SAMPLES) cannot deadlock, then
# runs make_database -> clusterseq -> genotype -> summarize -> makefasta.
# Submit: sbatch --job-name=mge_SL39_f_agg --dependency=afterok:<lastchunk>:<P0> \
#   scripts/submit_sl_aggregate.sh SL39 f

set -euo pipefail
SL="${1:?need Sublineage}"; LVL="${2:?need f|d}"
case "$LVL" in f|d) ;; *) echo "level must be f|d" >&2; exit 1;; esac
POOL="/home/dca36/rds/rds-floto-bacterial-4k08a2yyQLw/david/processed/mgefinder/_readpool/${SL}"
WD="/home/dca36/rds/rds-floto-bacterial-4k08a2yyQLw/david/processed/mgefinder/${SL}/ref_${LVL}"

echo "=== aggregate ${SL} ref_${LVL} ==="
echo "Job ${SLURM_JOB_ID:-NA}  Node $(hostname)  Started $(date)"

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
    --aggregate-only \
    --drop-file "${POOL}/dropped_samples.tsv" \
    --test-n -1 \
    --verbose \
    --jobs "$SLURM_CPUS_PER_TASK"

echo ">>> aggregate ${SL} ref_${LVL} done: $(date)"
echo "    Results: ${WD}/results/"
