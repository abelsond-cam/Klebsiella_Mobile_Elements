#!/usr/bin/env python3
"""Generate limited mgefinder_dataset.txt for testing (wrapper around generate_mgefinder_dataset.py)."""

import argparse
import subprocess
import sys
import tempfile
from pathlib import Path

import pandas as pd
import yaml


def main():
    parser = argparse.ArgumentParser(
        description="Generate limited MGEfinder dataset for testing"
    )
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--row", type=int, default=0)
    parser.add_argument("--test-n", type=int, required=True, help="Limit to first N comparison samples")
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    
    with open(args.config) as f:
        config = yaml.safe_load(f)
    
    data_dir = Path(config.get("data_dir", config["wd"])).resolve()
    refcomp_path = data_dir / config["reference_comparison_sets"]
    
    # Read and modify reference comparison set
    df = pd.read_csv(refcomp_path, sep="\t")
    row_data = df.iloc[args.row].copy()
    
    # Limit comparison set
    comparison_ids = [s.strip() for s in str(row_data["mge_comparison_set"]).split(",") if s.strip()]
    limited_ids = comparison_ids[:args.test_n]
    row_data["mge_comparison_set"] = ",".join(limited_ids)
    
    print(f">>> Limited comparison set: {len(comparison_ids)} → {len(limited_ids)} samples")
    
    # Create temporary TSV with limited data
    with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as tmp:
        tmp_df = pd.DataFrame([row_data])
        tmp_df.to_csv(tmp.name, sep='\t', index=False)
        temp_tsv = tmp.name
    
    try:
        # Update config to point to temp TSV
        temp_config = config.copy()
        temp_config["reference_comparison_sets"] = temp_tsv
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.yaml', delete=False) as tmp:
            yaml.dump(temp_config, tmp)
            temp_config_path = tmp.name
        
        # Call original script with modified config
        cmd = [
            "python3", "src/generate_mgefinder_dataset.py",
            "--config", temp_config_path,
            "--row", "0",
            "--out", str(args.out)
        ]
        
        subprocess.run(cmd, check=True)
        
    finally:
        # Clean up temp files
        Path(temp_tsv).unlink(missing_ok=True)
        Path(temp_config_path).unlink(missing_ok=True)


if __name__ == "__main__":
    main()