#!/bin/bash
# Orchestrate the decoupled chunked shared-pool run for one Sublineage:
#   producer (ahead) -> per-chunk consumers for ref_f AND ref_d (share the pool)
#   -> per-level cohort aggregation -> final pool GC.
# Usage:  scripts/orchestrate_sl.sh SL39 [f d]
# Run from the repo root on the LOGIN node (only submits jobs; no heavy I/O).

set -euo pipefail
SL="${1:?Usage: orchestrate_sl.sh <Sublineage> [levels...]}"; shift || true
LEVELS=("${@:-f d}"); [ $# -eq 0 ] && LEVELS=(f d)
POOL="/home/dca36/rds/rds-floto-bacterial-4k08a2yyQLw/david/processed/mgefinder/_readpool/${SL}"
cd "$(dirname "$0")/.."

# 1. Materialise chunk files + get chunk count K (light; no downloads).
INFO=$(micromamba run -n snakemake python3 src/run_readpool_producer.py \
        --config config/config.yaml --sublineage "$SL" \
        --pool-dir "$POOL" --chunk-size 100 --print-chunks 2>&1 | tee /dev/stderr)
K=$(echo "$INFO" | grep -oE 'in [0-9]+ chunk' | grep -oE '[0-9]+')
[ -n "$K" ] || { echo "could not determine chunk count" >&2; exit 1; }
echo ">>> ${SL}: ${K} chunk(s); levels: ${LEVELS[*]}"

# 2. Producer (runs ahead of all consumers).
P0=$(sbatch --parsable --job-name="mge_pool_${SL}" \
      scripts/submit_readpool_producer.sh "$SL")
echo "producer P0=$P0"

# 3. Per-level serial chunk chains (each level independent; share the pool).
declare -A LAST
for LVL in "${LEVELS[@]}"; do
    prev=""
    for ((k=0; k<K; k++)); do
        if [ -z "$prev" ]; then dep="after:$P0"; else dep="afterok:$prev"; fi
        jid=$(sbatch --parsable --job-name="mge_${SL}_${LVL}_c${k}" \
              --dependency="$dep" \
              scripts/submit_sl_consumer.sh "$SL" "$LVL" "$k")
        echo "  consumer ${SL} ${LVL} chunk ${k} -> $jid (dep $dep)"
        prev="$jid"
    done
    LAST[$LVL]="$prev"
done

# 4. Per-level cohort aggregation (after that level's last chunk + producer).
declare -A AGG
for LVL in "${LEVELS[@]}"; do
    jid=$(sbatch --parsable --job-name="mge_${SL}_${LVL}_agg" \
          --dependency="afterok:${LAST[$LVL]}:${P0}" \
          scripts/submit_sl_aggregate.sh "$SL" "$LVL")
    AGG[$LVL]="$jid"
    echo "  aggregate ${SL} ${LVL} -> $jid"
done

# 5. Final pool GC after all aggregations.
gcdep="afterok"; for LVL in "${LEVELS[@]}"; do gcdep="${gcdep}:${AGG[$LVL]}"; done
GZ=$(sbatch --parsable --job-name="mge_poolgc_${SL}" --dependency="$gcdep" \
      scripts/submit_readpool_gc.sh "$SL" --final)
echo "  final GC -> $GZ"
echo ">>> ${SL} job graph submitted."
