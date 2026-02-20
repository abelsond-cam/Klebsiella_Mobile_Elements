#!/usr/bin/env python3
"""Pipeline runner: Snakemake from snakemake env; MGEfinder tools via conda in Snakefile."""

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd
import yaml


def load_config(config_path: Path) -> dict:
    """Load and validate config."""
    with open(config_path) as f:
        config = yaml.safe_load(f)
    
    required = ["wd", "data_dir", "reference_comparison_sets"]
    for key in required:
        if key not in config:
            raise SystemExit(f"Error: missing required config key: {key}")
    
    return config


def get_reference_info(config: dict, row: int = 0, test_n: int = None) -> tuple:
    """Get reference name and comparison set IDs from TSV."""
    data_dir = Path(config.get("data_dir", config["wd"])).resolve()
    refcomp_path = data_dir / config["reference_comparison_sets"]
    
    if not refcomp_path.exists():
        raise SystemExit(f"Error: reference_comparison_sets not found: {refcomp_path}")
    
    df = pd.read_csv(refcomp_path, sep="\t")
    if row >= len(df):
        raise SystemExit(f"Error: row {row} out of range")
    
    row_data = df.iloc[row]
    ref_name = row_data["reference_sample_name"]
    comparison_ids = [s.strip() for s in str(row_data["mge_comparison_set"]).split(",") if s.strip()]
    
    # Limit samples for testing
    if test_n is not None and test_n > 0:
        comparison_ids = comparison_ids[:test_n]
        print(f">>> TEST MODE: Limited to {len(comparison_ids)} samples (--test-n {test_n})")
    
    return ref_name, comparison_ids


def run_with_env(cmd: list, description: str, env_name: str) -> None:
    """Run command with specified micromamba environment activation."""
    print(f">>> {description}")
    print(f"    Environment: {env_name}")
    print(f"    Command: {' '.join(cmd)}")
    
    # Build full command with environment activation
    full_cmd = [
        "bash", "-c", 
        f'eval "$(micromamba shell hook --shell bash)" && micromamba activate {env_name} && {" ".join(cmd)}'
    ]
    
    try:
        subprocess.run(full_cmd, check=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running {description}: {e}", file=sys.stderr)
        sys.exit(1)


def check_first_inputs(wd: Path, merged: dict, ref_name: str) -> None:
    """Verify inputs for the first (sample, genome) exist; exit with clear message if any are missing."""
    data_dir = Path(merged["data_dir"])
    meta = wd / merged.get("mgefinder_dataset", "mgefinder_dataset.txt")
    if not meta.exists():
        print(f"Error: dataset not found: {meta}", file=sys.stderr)
        sys.exit(1)
    with open(meta) as f:
        next(f)
        first = f.readline()
    if not first:
        print("Error: dataset is empty (no samples).", file=sys.stderr)
        sys.exit(1)
    parts = first.strip().split("\t")
    data_dir_str, sample_id, _sample_name, _gff, contigs = parts[0], parts[1], parts[2], parts[3], parts[4]
    has_fastq_cols = len(parts) >= 7
    if has_fastq_cols:
        fastq1_path, fastq2_path = Path(parts[5]), Path(parts[6])
    sample_data_dir = Path(data_dir_str)

    assemblies_dir = data_dir / merged.get("assemblies_dir", "raw/assemblies")
    ref_path = None
    for ext in (".fna", ".fa", ".fna.gz", ".fa.gz"):
        cand = assemblies_dir / (ref_name + ext)
        if cand.exists():
            ref_path = cand
            break
    if ref_path is None:
        print(f"Error: reference assembly not found (looked for {ref_name}{{.fna,.fa,.fna.gz,.fa.gz}} under {assemblies_dir})", file=sys.stderr)
        sys.exit(1)

    if not Path(contigs).exists():
        print(f"Error: sample assembly not found: {contigs}", file=sys.stderr)
        sys.exit(1)

    if not has_fastq_cols:
        print("Error: dataset has no fastq1/fastq2 columns; re-run dataset generation to get 7-column format.", file=sys.stderr)
        sys.exit(1)
    f1, f2 = fastq1_path, fastq2_path
    if not f1.exists() or not f2.exists():
        print(f"Error: FASTQ pair not found for sample {sample_id}:", file=sys.stderr)
        if not f1.exists():
            print(f"  missing: {f1}", file=sys.stderr)
        if not f2.exists():
            print(f"  missing: {f2}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Run MGEfinder pipeline")
    parser.add_argument("--config", type=Path, default=Path("config/config.yaml"))
    parser.add_argument("--row", type=int, default=0, help="Row from reference_comparison_sets")
    parser.add_argument("--fetch-only", action="store_true", help="Only fetch data, don't run pipeline")
    parser.add_argument("--pipeline-only", action="store_true", help="Only run pipeline, skip fetch")
    parser.add_argument("--verbose", action="store_true", help="Use verbose Snakemake output")
    parser.add_argument("--jobs", "-j", type=int, default=1, help="Snakemake jobs")
    parser.add_argument("--test-n", type=int, help="TEST MODE: Limit to first N comparison samples")
    parser.add_argument("--dry-run", action="store_true", help="Snakemake dry run (-n)")
    parser.add_argument("--skip-download", action="store_true", help="Skip FASTQ download step")
    args = parser.parse_args()
    
    # Environment names (micromamba; no --use-conda so Snakemake never invokes conda)
    SNAKE_ENV = "snakemake"        # Data prep (validate, dataset generation)
    MGEFINDER_ENV = "mgefinder_env"  # Snakemake + all rule jobs (bwa, bowtie2, mgefinder); you manage this env yourself
    FASTQ_ENV = "fastq-dl"         # FASTQ downloads only
    
    config = load_config(args.config)
    wd = Path(config["wd"]).resolve()
    ref_name, comparison_ids = get_reference_info(config, args.row, args.test_n)
    
    print(f">>> Pipeline setup: wd={wd}, reference={ref_name}")
    print(f">>> Comparison samples: {len(comparison_ids)} ({', '.join(comparison_ids[:3])}{'...' if len(comparison_ids) > 3 else ''})")
    print(f">>> Data prep: {SNAKE_ENV}; Snakemake (no conda): {MGEFINDER_ENV}")
    
    # PHASE 1: Data preparation (snakemake environment)
    if not args.pipeline_only:
        print("\n" + "="*70)
        print("PHASE 1: DATA PREPARATION")
        print("="*70)
        
        # Validate reference assembly exists
        run_with_env([
            "python3", "src/validate_reference.py", 
            "--config", str(args.config), "--row", str(args.row)
        ], "Validate reference assembly", SNAKE_ENV)
        
        # PHASE 2: FASTQ downloads (fastq-dl environment; TODO: fastq-dl env may lack pyyaml - consider adding deps to snakemake env and running from there, or use --skip-download when samples exist)
        print(">>> FASTQ downloads (fastq-dl environment; TODO: fastq-dl env may lack pyyaml - consider adding deps to snakemake env and running from there, or use --skip-download when samples exist)")
        if comparison_ids and not args.skip_download:
            print("\n" + "="*70)
            print("PHASE 2: FASTQ DOWNLOADS")
            print("="*70)
            
            run_with_env([
                "python3", "src/run_fastq_download.py",
                "--config", str(args.config),
                "--ids"] + comparison_ids,
                f"Download FASTQ for {len(comparison_ids)} samples", FASTQ_ENV
            )
        elif args.skip_download:
            print("\n>>> SKIPPING FASTQ downloads (--skip-download)")
        
        # Back to data preparation
        print("\n" + "="*70)
        print("PHASE 1b: DATASET GENERATION")
        print("="*70)
        
        # Generate dataset (limited for testing if needed)
        dataset_path = wd / "mgefinder_dataset.txt"
        if args.test_n is not None:
            run_with_env([
                "python3", "src/generate_mgefinder_dataset_limited.py",
                "--config", str(args.config),
                "--row", str(args.row),
                "--test-n", str(args.test_n),
                "--out", str(dataset_path)
            ], f"Generate LIMITED dataset ({args.test_n} samples)", SNAKE_ENV)
        else:
            run_with_env([
                "python3", "src/generate_mgefinder_dataset.py",
                "--config", str(args.config),
                "--row", str(args.row),
                "--out", str(dataset_path)
            ], f"Generate dataset", SNAKE_ENV)
        
        if args.fetch_only:
            print(f"\n>>> Fetch completed. Dataset: {dataset_path}")
            return
    
    # PHASE 3: Snakemake pipeline (snakemake environment)
    if not args.fetch_only:
        print("\n" + "="*70)
        print("PHASE 3: SNAKEMAKE PIPELINE")
        print("="*70)
        
        # Merge main config with genomes and write a single config so Snakemake gets full config when using --directory.
        # Use resolved paths so Snakefile looks for inputs/outputs in the same place as --directory (avoids path mismatch e.g. /home/dca36/rds vs /rds/project).
        merged = dict(config)
        merged["wd"] = str(wd)
        merged["data_dir"] = str(Path(config.get("data_dir", config["wd"])).resolve())
        merged["genomes"] = [ref_name]
        merged_yaml = Path("config/.mge_merged_config.yaml")
        merged_yaml.parent.mkdir(parents=True, exist_ok=True)
        with open(merged_yaml, "w") as f:
            yaml.dump(merged, f, default_flow_style=False, sort_keys=False)
        print(f">>> Created {merged_yaml} (config + genomes)")
        
        if not args.dry_run:
            check_first_inputs(wd, merged, ref_name)
        
        # Run Snakemake from working environment (single config file so all keys are present)
        config_abs = merged_yaml.resolve()
        snake_cmd = [
            "snakemake",
            "--configfile", str(config_abs),
            "--directory", str(wd),
            "-j", str(args.jobs),
            "all"
        ]
        
        if args.verbose:
            snake_cmd.extend(["--printshellcmds", "-p"])
        
        if args.dry_run:
            snake_cmd.append("-n")
        
        action = "DRY RUN Snakemake pipeline" if args.dry_run else "Run Snakemake pipeline"
        run_with_env(snake_cmd, action, MGEFINDER_ENV)
        
        if args.dry_run:
            print(f"\n>>> Dry run completed - no files were created")
            print(f"    To run for real: remove --dry-run flag")
        else:
            print(f"\n>>> Pipeline completed successfully!")
            print(f"    Results in: {wd}/results/{ref_name}/")


if __name__ == "__main__":
    main()
