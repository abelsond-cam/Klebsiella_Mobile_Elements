#!/bin/bash
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=76
#SBATCH --partition=icelake-himem
#SBATCH --account=FLOTO-PROJECT-K-SL2-CPU
#SBATCH --output=mge_%x_%j.log
#SBATCH --error=mge_%x_%j.err

# Thin delegator -> the SELF-HEALING fixed salvage (auto-reconciles the cohort
# to whatever actually downloaded; no manual drop-file recovery needed).
# Kept so existing references / habits still work.
#   $1=Sublineage $2=ref {f|d} $3=N
exec bash scripts/submit_sl_salvage_fixed.sh "$@"
