#!/bin/bash
#SBATCH --time=36:00:00
#SBATCH --cpus-per-task=76
#SBATCH --partition=icelake-himem
#SBATCH --account=FLOTO-PROJECT-K-SL2-CPU
#SBATCH --output=mge_%x_%j.log
#SBATCH --error=mge_%x_%j.err

# Generic full-SL MGEfinder run, parameterised by Sublineage.
#   $1 = Sublineage (e.g. SL258). Required.
# Same external reference as SL39 (config reference_sample_name), stream_fastq,
# all samples in the SL (--test-n -1), isolated output dir .../mgefinder/<SL>.
# Set the job name at submit so logs/dirs are per-SL:
#   sbatch --job-name=mge_SL258 scripts/submit_sl_full.sh SL258
# Large SLs (>~1500 samples) may hit the 36h wall mid-run; per-sample TSVs in
# wd are NOT temp() (only FASTQ is), so a resubmit resumes via Snakemake.

set -euo pipefail
SL="${1:?Usage: sbatch --job-name=mge_<SL> submit_sl_full.sh <Sublineage>}"
OUTDIR="/home/dca36/rds/rds-floto-bacterial-4k08a2yyQLw/david/processed/mgefinder/${SL}"

echo "=== MGEfinder FULL run: ${SL} ==="
echo "Job ID: ${SLURM_JOB_ID:-NA}  Node: $(hostname)  Started: $(date)"
echo "Output dir: ${OUTDIR}  CPUs/jobs: ${SLURM_CPUS_PER_TASK}  Walltime: 36h"

cd "$SLURM_SUBMIT_DIR"
eval "$(micromamba shell hook --shell bash)"
micromamba activate mgefinder_env

python3 src/run_pipeline.py \
    --config config/config.yaml \
    --sublineage "$SL" \
    --output-dir "$OUTDIR" \
    --stream-fastq \
    --test-n -1 \
    --verbose \
    --jobs "$SLURM_CPUS_PER_TASK"

echo ">>> ${SL} run completed: $(date)"
