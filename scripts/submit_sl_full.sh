#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=76
#SBATCH --partition=icelake-himem
#SBATCH --account=FLOTO-PROJECT-K-SL2-CPU
#SBATCH --output=mge_%x_%j.log
#SBATCH --error=mge_%x_%j.err

# Generic full-SL MGEfinder run, parameterised by Sublineage + reference level.
#   $1 = Sublineage (e.g. SL258). Required.
#   $2 = reference level: f = ref_f (MGH 78578 basal, DEFAULT) | d = ref_d
#        (per-SL level-d). Passed as --reference-level so the reference is
#        DETERMINISTIC for this job regardless of later config/code changes.
# Reference is resolved per-SL from best_reference_per_sample.csv (joined by
# sample_accession, QC'd unique). stream_fastq, all samples (--test-n -1).
# Output dir is isolated per SL AND per level so f/d runs never collide:
#   .../mgefinder/<SL>/ref_<level>
# Submit both waves as pending, e.g.:
#   sbatch --job-name=mge_SL39_f scripts/submit_sl_full.sh SL39 f
#   sbatch --job-name=mge_SL39_d scripts/submit_sl_full.sh SL39 d
# Large SLs (>~1500 samples) may hit the 36h wall mid-run; per-sample TSVs in
# wd are NOT temp() (only FASTQ is), so a resubmit resumes via Snakemake.

set -euo pipefail
SL="${1:?Usage: sbatch --job-name=mge_<SL>_<lvl> submit_sl_full.sh <Sublineage> <f|d>}"
LVL="${2:-f}"
case "$LVL" in f|d) ;; *) echo "reference level must be f or d, got '$LVL'" >&2; exit 1;; esac
OUTDIR="/home/dca36/rds/rds-floto-bacterial-4k08a2yyQLw/david/processed/mgefinder/${SL}/ref_${LVL}"

echo "=== MGEfinder FULL run: ${SL} (reference-level ${LVL}) ==="
echo "Job ID: ${SLURM_JOB_ID:-NA}  Node: $(hostname)  Started: $(date)"
echo "Output dir: ${OUTDIR}  CPUs/jobs: ${SLURM_CPUS_PER_TASK}  Walltime: 36h"

cd "$SLURM_SUBMIT_DIR"
eval "$(micromamba shell hook --shell bash)"
micromamba activate mgefinder_env

python3 src/run_pipeline.py \
    --config config/config.yaml \
    --sublineage "$SL" \
    --reference-level "$LVL" \
    --output-dir "$OUTDIR" \
    --stream-fastq \
    --test-n -1 \
    --verbose \
    --jobs "$SLURM_CPUS_PER_TASK"

echo ">>> ${SL} (ref-level ${LVL}) run completed: $(date)"
