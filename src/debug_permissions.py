#!/usr/bin/env python3
"""Debug file permissions issue."""

import argparse
import sys
from pathlib import Path

import pandas as pd
import yaml


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--row", type=int, default=0)
    args = parser.parse_args()
    
    with open(args.config) as f:
        config = yaml.safe_load(f)
    
    data_dir = Path(config.get("data_dir", config["wd"])).resolve()
    assemblies_dir = data_dir / config["assemblies_dir"]
    refcomp_path = data_dir / config["reference_comparison_sets"]
    
    print(f"Data dir: {data_dir}")
    print(f"Assemblies dir: {assemblies_dir}")
    print(f"Assemblies dir exists: {assemblies_dir.exists()}")
    
    if assemblies_dir.exists():
        try:
            files = list(assemblies_dir.glob("*.gz"))
            print(f"First few .gz files: {[f.name for f in files[:5]]}")
        except PermissionError as e:
            print(f"Cannot list assemblies dir: {e}")
    
    # Check reference comparison sets
    df = pd.read_csv(refcomp_path, sep="\t")
    row_data = df.iloc[args.row]
    comparison_ids = [s.strip() for s in str(row_data["mge_comparison_set"]).split(",") if s.strip()]
    
    print(f"\nFirst few comparison IDs: {comparison_ids[:3]}")
    
    # Try to load metadata and map first ID
    def load_metadata_local(metadata_path: Path):
        """Local copy of metadata loading."""
        df = pd.read_csv(metadata_path, sep="\t", low_memory=False)
        for col in ("Sample", "sample_accession"):
            if col not in df.columns:
                raise SystemExit(f"Error: metadata missing column '{col}'")
        return df
    
    metadata_path = data_dir / config["metadata_file"]
    meta = load_metadata_local(metadata_path)
    meta = meta.astype({"sample_accession": str})
    accession_to_sample = meta.set_index("sample_accession")["Sample"].to_dict()
    
    first_id = comparison_ids[0]
    sample_name = accession_to_sample.get(first_id)
    print(f"First ID {first_id} → Sample: {sample_name}")
    
    # Try to access the problematic file
    for ext in [".fna.gz", ".fa.gz", ".fna", ".fa"]:
        test_path = assemblies_dir / f"{sample_name}{ext}"
        print(f"Testing: {test_path}")
        try:
            exists = test_path.exists()
            print(f"  Exists: {exists}")
            if exists:
                stat = test_path.stat()
                print(f"  Size: {stat.st_size} bytes")
                break
        except PermissionError as e:
            print(f"  Permission error: {e}")
        except Exception as e:
            print(f"  Other error: {e}")


if __name__ == "__main__":
    main()