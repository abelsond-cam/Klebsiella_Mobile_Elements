# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A Snakemake pipeline that runs the **MGEfinder** end-to-end workflow to detect mobile genetic
elements (insertion sequences) in *Klebsiella* genomes. For each "comparison set" (one reference
genome plus a list of sample accessions) it aligns sample short reads to the reference, then runs
MGEfinder's find → pair → inferseq → makedatabase → clusterseq → genotype → summarize → makefasta.

Note: the parent `~/CLAUDE.md` describes the wider CSD3/HPC environment. **Unlike other
`~/workspace/` projects, this one is NOT run via `uv`** — it uses micromamba environments and
Snakemake, driven by a Makefile.

## Environments

The pipeline deliberately splits work across three micromamba environments (no `--use-conda`;
the driver activates the right env per phase):

| Env | Used for |
|-----|----------|
| `snakemake` | Data prep: config loading, dataset generation, validation |
| `mgefinder_env` | Snakemake itself **and** all rule jobs (bwa, bowtie2, mgefinder). You maintain this env manually. |
| `fastq-dl` | FASTQ downloads only |

`mgefinder_env` is fragile. MGEfinder 1.0.6 is Python-2-era and breaks against modern
numpy/scipy/pandas/Biopython. Before assuming a pipeline bug, check
[docs/mgefinder_environment.md](docs/mgefinder_environment.md) — it documents the required
source patches to the installed `mgefinder` package (`dependencies.py`, `formatbam.py`,
`pair.py`, `inferseqassembly.py`) and version pins (Biopython <1.78, pandas <2). `bwa` and
`bowtie2` can vanish from the env after pip installs; re-check with
`micromamba run -n mgefinder_env bwa` / `bowtie2-build --version`.

## Common commands

```bash
make help                       # list all targets
make test-dry-skip-dl test-n=1  # safe login-node dry run, 1 sample, no download (start here)
make test test-n=5              # real run, 5 samples, skips FASTQ download
make dry-run                    # full dry run (all samples in set)
make data                       # full real run
make submit-hpc                 # sbatch scripts/submit_hpc.sh
make lint                       # flake8 src (mgefinder_env)
make test_environment           # scripts/diagnose_env.py — checks mgefinder_env health
make recreate-symlinks          # rebuild assembly symlinks after a path change
```

`du` and any heavy I/O must go through SLURM, not the login node (see parent `~/CLAUDE.md`).

## Architecture

**Entrypoint chain:** `make` → `src/run_pipeline.py` → Snakemake → `Snakefile` →
`workflow/Snakefile` → `workflow/mgefinder.end2end.snakefile` (all real rules live in the last
file; the top two just `include:` downward per Snakemake convention).

**`src/run_pipeline.py` runs three phases**, switching micromamba env per phase:
1. FASTQ downloads (`run_fastq_download.py`, `fastq-dl` env) — skippable with `--skip-download`.
2. Dataset generation (`generate_mgefinder_dataset.py`, or `_limited.py` when `--test-n` is set,
   `snakemake` env) — writes `<wd>/mgefinder_dataset.txt`, a 7-column TSV:
   `data_dir, sample_id, sample_name, gff, contigs, fastq1, fastq2`.
3. Snakemake (`mgefinder_env`). The driver merges `config.yaml` + the resolved `genomes` list
   into `config/.mge_merged_config.yaml` (gitignored, regenerated each run) and passes that
   single file plus `--directory <wd>`.

**Config model** ([config/config.yaml](config/config.yaml) is the single source of truth):
- `data_dir` = where the pipeline **reads** (assemblies, fastq, metadata,
  `reference_comparison_sets.tsv`).
- `wd` = where the pipeline **writes** everything (`00.genome`, `00.assembly`, `00.bam`, `bwa`,
  `mgefinder`, `database`, `results` are created under `wd`).
- `genomes` is NOT in `config.yaml` — it is injected at runtime from
  `reference_comparison_sets.tsv` (`reference_sample_name` column). Running Snakemake directly
  without going through `run_pipeline.py`/`make` will fail with a "'genomes' not found" error.
- `reference_comparison_sets.tsv` columns used: `reference_sample_name` (the reference genome)
  and `mge_comparison_set` (comma-separated sample accessions). `n_set_to_run` / `--row`
  selects which row; one Snakemake run per row.
- MGEfinder parameters under `find:`, `pair:`, `inferseq_*:`, etc. mirror MGEfinder's upstream
  `denovo.original.config.yml` — keep them in sync if updating.

**Critical Snakemake gotcha (absolute vs relative paths):** with `--directory` set, a rule's
output path and the consuming rule's input path can resolve to the same file on disk but be
different *strings* in the DAG, so Snakemake won't link them and raises MissingInputException.
Always build paths with the `join(<DIR>, ...)` helpers (e.g. `BAM_DIR`, `GENOME_DIR`) on both
sides — never mix a relative literal with a `join()`ed absolute path. See
[docs/BAM_and_find_rule.md](docs/BAM_and_find_rule.md).

The Snakefile prints `MGEfinder paths: data_dir (read) = ... | wd (write) = ...` at parse time
— check this first when debugging path issues. `# --- SEGMENT: ... ---` comments delimit
pipeline stages for isolated testing.
