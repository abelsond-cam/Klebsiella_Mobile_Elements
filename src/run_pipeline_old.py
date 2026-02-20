#!/usr/bin/env python3
"""Unified pipeline runner: fetch data + run MGEfinder pipeline."""

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


def run_command(cmd: list, description: str) -> None:
    """Run command and handle errors."""
    print(f">>> {description}")
    print(f"    Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, text=True)
        return result
    except subprocess.CalledProcessError as e:
        print(f"Error running {description}: {e}", file=sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(description="Run MGEfinder pipeline")
    parser.add_argument("--config", type=Path, default=Path("config/config.yaml"))
    parser.add_argument("--row", type=int, default=0, help="Row from reference_comparison_sets")
    parser.add_argument("--fetch-only", action="store_true", help="Only fetch data, don't run pipeline")
    parser.add_argument("--pipeline-only", action="store_true", help="Only run pipeline, skip fetch")
    parser.add_argument("--verbose", action="store_true", help="Use verbose Snakemake output")
    parser.add_argument("--jobs", "-j", type=int, default=1, help="Snakemake jobs")
    parser.add_argument("--test-n", type=int, help="TEST MODE: Limit to first N comparison samples (for login node testing)")
    parser.add_argument("--dry-run", action="store_true", help="Snakemake dry run (-n) - show what would be done")
    args = parser.parse_args()
    
    config = load_config(args.config)
    wd = Path(config["wd"]).resolve()
    ref_name, comparison_ids = get_reference_info(config, args.row, args.test_n)
    
    print(f">>> Pipeline setup: wd={wd}, reference={ref_name}")
    print(f">>> Comparison samples: {len(comparison_ids)} ({', '.join(comparison_ids[:3])}{'...' if len(comparison_ids) > 3 else ''})")
    
    # Fetch data phase
    if not args.pipeline_only:
        print("\n" + "="*60)
        print("FETCH DATA PHASE")
        print("="*60)
        
        # Validate reference assembly exists
        run_command([
            "python3", "src/validate_reference.py", 
            "--config", str(args.config), "--row", str(args.row)
        ], "Validate reference assembly")
        
        # Download FASTQ if needed
        if comparison_ids:
            run_command([
                "python3", "src/run_fastq_download.py",
                "--config", str(args.config),
                "--ids"] + comparison_ids,
                f"Download FASTQ for {len(comparison_ids)} samples"
            )
        
        # Generate dataset (limited for testing if needed)
        dataset_path = wd / "mgefinder_dataset.txt"
        if args.test_n is not None:
            # Use limited dataset generator for testing
            run_command([
                "python3", "src/generate_mgefinder_dataset_limited.py",
                "--config", str(args.config),
                "--row", str(args.row),
                "--test-n", str(args.test_n),
                "--out", str(dataset_path)
            ], f"Generate LIMITED dataset ({args.test_n} samples) -> {dataset_path}")
        else:
            run_command([
                "python3", "src/generate_mgefinder_dataset.py",
                "--config", str(args.config),
                "--row", str(args.row),
                "--out", str(dataset_path)
            ], f"Generate dataset -> {dataset_path}")
        
        if args.fetch_only:
            print(f"\n>>> Fetch completed. Dataset: {dataset_path}")
            return
    
    # Pipeline phase  
    if not args.fetch_only:
        print("\n" + "="*60)
        print("PIPELINE PHASE")
        print("="*60)
        
        # Create genomes config
        genomes_yaml = Path(".mge_genomes.yaml")
        genomes_yaml.write_text(f"genomes: [{ref_name}]\n")
        print(f">>> Created {genomes_yaml}")
        
        # Run Snakemake
        snake_cmd = [
            "snakemake",
            "--configfile", str(args.config),
            "--configfile", str(genomes_yaml),
            "--directory", str(wd),
            "-j", str(args.jobs),
            "all"
        ]
        
        if args.verbose:
            snake_cmd.extend(["--printshellcmds", "-p"])
        
        if args.dry_run:
            snake_cmd.append("-n")
        
        action = "DRY RUN Snakemake pipeline" if args.dry_run else "Run Snakemake pipeline"
        run_command(snake_cmd, action)
        
        if args.dry_run:
            print(f"\n>>> Dry run completed - no files were created")
            print(f"    To run for real: remove --dry-run flag")
        else:
            print(f"\n>>> Pipeline completed successfully!")
            print(f"    Results in: {wd}/results/{ref_name}/")


if __name__ == "__main__":
    main()