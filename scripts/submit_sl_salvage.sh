#!/bin/bash
#SBATCH --time=07:30:00
#SBATCH --cpus-per-task=76
#SBATCH --partition=icelake-himem
#SBATCH --account=FLOTO-PROJECT-K-SL2-CPU
#SBATCH --output=mge_%x_%j.log
#SBATCH --error=mge_%x_%j.err

# SHORT salvage run before HPC maintenance. $1=Sublineage $2=ref {f|d} $3=N.
# Uses the validated hardened STREAMING path (one self-consistent dataset
# drives both download_fastq targets and SAMPLES, so no producer/consumer
# mismatch) capped at N samples so it completes well inside the wall and
# yields a real N-sample cohort (clusterseq/genotype/makefasta) to analyse
# during the downtime. Isolated output dir (won't collide with the 36h
# decoupled jobs). Submit:
#   sbatch --job-name=mge_SL39_f_sv scripts/submit_sl_salvage.sh SL39 f 60

set -euo pipefail
SL="${1:?need Sublineage}"; LVL="${2:?need f|d}"; N="${3:?need sample cap N}"
case "$LVL" in f|d) ;; *) echo "level must be f|d" >&2; exit 1;; esac
WD="/home/dca36/rds/rds-floto-bacterial-4k08a2yyQLw/david/processed/mgefinder/${SL}/ref_${LVL}_salvage${N}"

echo "=== SALVAGE ${SL} ref_${LVL} N=${N} ==="
echo "Job ${SLURM_JOB_ID:-NA}  Node $(hostname)  Started $(date)  Wall 7.5h"

cd "$SLURM_SUBMIT_DIR"
eval "$(micromamba shell hook --shell bash)"
micromamba activate mgefinder_env

python3 src/run_pipeline.py \
    --config config/config.yaml \
    --sublineage "$SL" \
    --reference-level "$LVL" \
    --output-dir "$WD" \
    --stream-fastq \
    --test-n "$N" \
    --verbose \
    --jobs "$SLURM_CPUS_PER_TASK"

echo ">>> SALVAGE ${SL} ref_${LVL} done: $(date)"
echo "    Results: ${WD}/results/"
