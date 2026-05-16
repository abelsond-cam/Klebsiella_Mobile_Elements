#!/bin/bash
#SBATCH --job-name=mge_poolgc
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=1
#SBATCH --partition=icelake
#SBATCH --account=FLOTO-PROJECT-K-SL2-CPU
#SBATCH --output=mge_%x_%j.log
#SBATCH --error=mge_%x_%j.err

# Pool GC. $1=Sublineage [$2=--final]. Deletes reads/<acc>_{1,2}.fastq.gz only
# when BOTH consumers (ref_f AND ref_d) have a consumed/<tag>/<acc>.done — i.e.
# no consumer still needs that sample's reads. Safe to run repeatedly. --final
# (after both aggregations) also removes the empty pool scaffolding.
# Submit: sbatch --dependency=afterok:<Af>:<Ad> scripts/submit_readpool_gc.sh SL39 --final

set -euo pipefail
SL="${1:?need Sublineage}"; FINAL="${2:-}"
POOL="/home/dca36/rds/rds-floto-bacterial-4k08a2yyQLw/david/processed/mgefinder/_readpool/${SL}"
[ -d "$POOL/reads" ] || { echo "no pool reads dir: $POOL/reads"; exit 0; }

echo "=== pool GC ${SL} ${FINAL} === $(date)"
freed=0
shopt -s nullglob
for r1 in "$POOL"/reads/*_1.fastq.gz; do
    base="$(basename "$r1")"; acc="${base%_1.fastq.gz}"
    if [ -f "$POOL/consumed/ref_f/${acc}.done" ] && [ -f "$POOL/consumed/ref_d/${acc}.done" ]; then
        rm -f "$POOL/reads/${acc}_1.fastq.gz" "$POOL/reads/${acc}_2.fastq.gz"
        freed=$((freed+1))
    fi
done
echo ">>> GC freed reads for ${freed} fully-consumed sample(s)"

if [ "$FINAL" = "--final" ]; then
    rm -rf "$POOL/reads"
    echo ">>> --final: removed $POOL/reads (sentinels + dropped_samples.tsv kept for audit)"
fi
echo ">>> pool GC ${SL} done: $(date)"
