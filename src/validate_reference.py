#!/usr/bin/env python3
"""Validate that reference genome from TSV has matching assembly file."""

import argparse
import sys
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Validate reference genome has matching assembly"
    )
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--row", type=int, default=0)
    args = parser.parse_args()

    try:
        import yaml
    except ImportError:
        raise SystemExit("Error: PyYAML required")

    with open(args.config) as f:
        config = yaml.safe_load(f)

    data_dir = Path(config.get("data_dir", config["wd"])).resolve()
    assemblies_dir = data_dir / config["assemblies_dir"]
    refcomp_path = data_dir / config["reference_comparison_sets"]

    if not refcomp_path.exists():
        raise SystemExit(f"Error: reference_comparison_sets not found: {refcomp_path}")

    refcomp = pd.read_csv(refcomp_path, sep="\t")
    if args.row >= len(refcomp):
        raise SystemExit(f"Error: row {args.row} out of range")

    reference_name = refcomp.iloc[args.row]["reference_sample_name"]
    
    # Check if assembly exists (same logic as get_reference_path)
    for ext in [".fna", ".fa", ".fna.gz", ".fa.gz"]:
        assembly_path = assemblies_dir / f"{reference_name}{ext}"
        if assembly_path.exists():
            print(f"✓ Reference '{reference_name}' → {assembly_path}")
            return

    print(f"✗ ERROR: No assembly found for reference '{reference_name}'", file=sys.stderr)
    print(f"  Looked in: {assemblies_dir}", file=sys.stderr)
    print(f"  Expected: {reference_name}.{{fna,fa,fna.gz,fa.gz}}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()