#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=76
#SBATCH --partition=icelake-himem
#SBATCH --account=FLOTO-PROJECT-K-SL2-CPU
#SBATCH --output=mge_%x_%j.log
#SBATCH --error=mge_%x_%j.err

# SELF-HEALING salvage. $1=Sublineage $2=ref {f|d} $3=N.
# One job, three phases, no manual drop-file ever:
#   1. Streaming per-sample (best-effort, --keep-going): every downloadable
#      sample gets find/pair/inferseq written; dead accessions just don't.
#      The cohort aggregation in this phase is EXPECTED to be unmet (Snakemake
#      exits non-zero) — that is fine, we ignore it.
#   2. Reconcile FROM DISK: the successful set = samples that produced
#      03.inferseq_overlap. Everything requested-but-missing -> auto drop-file.
#   3. Aggregate-only over the successful set -> makefasta. No reads needed.
# Result: a finished cohort over whatever actually downloaded, automatically.
# Submit: sbatch --job-name=mge_SL23_f_fix scripts/submit_sl_salvage_fixed.sh SL23 f 45

set -uo pipefail
SL="${1:?need Sublineage}"; LVL="${2:?need f|d}"; N="${3:?need sample cap N}"
case "$LVL" in f|d) ;; *) echo "level must be f|d" >&2; exit 1;; esac
BASE=/home/dca36/rds/rds-floto-bacterial-4k08a2yyQLw/david/processed/mgefinder
WD="$BASE/${SL}/ref_${LVL}_sv${N}fix"
DUMMY="/tmp/dummy_pool_${SL}_${LVL}"

echo "=== SELF-HEALING salvage ${SL} ref_${LVL} N=${N} ==="
echo "Job ${SLURM_JOB_ID:-NA}  Node $(hostname)  Started $(date)  WD=$WD"

cd "$SLURM_SUBMIT_DIR"
eval "$(micromamba shell hook --shell bash)"
micromamba activate mgefinder_env

echo ">>> PHASE 1: streaming per-sample (best-effort; aggregation may be unmet)"
python3 src/run_pipeline.py \
    --config config/config.yaml --sublineage "$SL" --reference-level "$LVL" \
    --output-dir "$WD" --stream-fastq --test-n "$N" --verbose \
    --jobs "$SLURM_CPUS_PER_TASK" || echo ">>> PHASE 1 exited non-zero (expected if any dead accessions) — continuing"

echo ">>> PHASE 2: reconcile successful set from disk"
G=$(basename "$(cat "$WD/mgefinder_reference_path.txt" 2>/dev/null)" 2>/dev/null \
      | sed -E 's/\.(fna|fa)(\.gz)?$//')
DS="$WD/mgefinder_dataset.txt"
if [ -z "$G" ] || [ ! -f "$DS" ]; then
    echo "ERROR: no dataset/reference in $WD — phase 1 produced nothing; aborting" >&2
    exit 1
fi
ls "$WD"/mgefinder/"$G"/*/03.inferseq_overlap.*."$G".tsv 2>/dev/null \
  | sed -E "s#.*/mgefinder/$G/([^/]+)/.*#\1#" | sort -u > "$WD/.done_names"
NDONE=$(wc -l < "$WD/.done_names")
echo ">>> per-sample complete: ${NDONE} sample(s)"
if [ "$NDONE" -eq 0 ]; then
    echo "ERROR: zero samples completed per-sample; nothing to aggregate" >&2
    exit 1
fi
DROP="$WD/dropped_samples.tsv"
{ echo -e "sample_accession\tattempts\tlast_error";
  awk -F'\t' 'NR>1{print $2"\t"$3}' "$DS" \
    | awk -F'\t' 'NR==FNR{d[$1]=1;next} !($2 in d){print $1"\t5\tno per-sample output (auto-reconcile)"}' \
        "$WD/.done_names" -; } > "$DROP"
NDROP=$(($(wc -l < "$DROP")-1))
echo ">>> auto drop-file: ${NDROP} requested-but-missing accession(s) -> $DROP"

echo ">>> PHASE 3: aggregate-only over the ${NDONE} successful samples"
python3 src/run_pipeline.py \
    --config config/config.yaml --sublineage "$SL" --reference-level "$LVL" \
    --output-dir "$WD" --read-pool "$DUMMY" --consumer-tag "ref_${LVL}" \
    --aggregate-only --drop-file "$DROP" --test-n "$N" --verbose \
    --jobs "$SLURM_CPUS_PER_TASK"

echo ">>> DONE ${SL} ref_${LVL}: $(ls "$WD"/results/*/04.makefasta.*.all_seqs.fna 2>/dev/null || echo 'NO makefasta') ($(date))"
echo "    Cohort = ${NDONE} samples (auto-reconciled). Results: ${WD}/results/"
