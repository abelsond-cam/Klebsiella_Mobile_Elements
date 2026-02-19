#!/bin/bash
#SBATCH --job-name=fastq_dl
#SBATCH --time=4:00:00
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --account=FLOTO-SL2-CPU

eval "$(micromamba shell hook --shell bash)"
micromamba activate fastq-dl

# Run from Klebsiella_Mobile_Elements project root so src/run_fastq_download.py and imports resolve
PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$PROJECT_ROOT"
python src/run_fastq_download.py "$@"
