#!/usr/bin/env python3
"""Generate mgefinder_dataset.txt from reference_comparison_sets.tsv and metadata.

Reads the selected row(s) from reference_comparison_sets.tsv (reference_sample_name,
mge_comparison_set), maps sample_accession -> Sample via metadata, discovers assembly
paths under assemblies_dir (supports .fna.gz, .fa.gz, .fna, .fa), and writes a
5-column TSV: data_dir, sample_id, sample_name, gff, contigs.

Usage:
  For first row only (n_set_to_run=1):
    python scripts/generate_mgefinder_dataset.py --config config/config.yaml --out mgefinder_dataset.txt
  For a specific row (when caller loops for n_set_to_run=-1):
    python scripts/generate_mgefinder_dataset.py --config config/config.yaml --row 0 --out mgefinder_dataset.txt
"""

import argparse
import sys
from pathlib import Path

import pandas as pd


# Assembly extensions to try (in order) when resolving path for a sample_name
ASSEMBLY_EXTENSIONS = (".fna.gz", ".fa.gz", ".fna", ".fa")


def discover_assembly_path(assemblies_dir: Path, sample_name: str) -> Path | None:
    """Return path to assembly file for sample_name, or None if not found."""
    for ext in ASSEMBLY_EXTENSIONS:
        p = assemblies_dir / f"{sample_name}{ext}"
        if p.exists():
            return p
    return None


def load_metadata(metadata_path: Path) -> pd.DataFrame:
    """Load metadata TSV; require Sample and sample_accession columns."""
    df = pd.read_csv(metadata_path, sep="\t", low_memory=False)
    for col in ("Sample", "sample_accession"):
        if col not in df.columns:
            raise SystemExit(f"Error: metadata missing column '{col}'")
    return df


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate mgefinder_dataset.txt from reference_comparison_sets.tsv"
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=Path("config/config.yaml"),
        help="Path to config YAML (default: config/config.yaml)",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("mgefinder_dataset.txt"),
        help="Output path for 5-column TSV (default: mgefinder_dataset.txt)",
    )
    parser.add_argument(
        "--row",
        type=int,
        default=0,
        metavar="N",
        help="Zero-based row index from reference_comparison_sets (default: 0)",
    )
    args = parser.parse_args()

    try:
        import yaml
    except ImportError:
        raise SystemExit("Error: PyYAML required. pip install pyaml")

    config_path = args.config
    if not config_path.exists():
        raise SystemExit(f"Error: config not found: {config_path}")

    with open(config_path) as f:
        config = yaml.safe_load(f)

    wd = Path(config["wd"]).resolve()
    assemblies_dir = wd / config["assemblies_dir"]
    fastq_dir = wd / config["fastq_dir"]
    metadata_path = wd / config["metadata_file"]
    refcomp_path = wd / config["reference_comparison_sets"]

    if not refcomp_path.exists():
        raise SystemExit(f"Error: reference_comparison_sets not found: {refcomp_path}")
    if not metadata_path.exists():
        raise SystemExit(f"Error: metadata not found: {metadata_path}")

    refcomp = pd.read_csv(refcomp_path, sep="\t")
    if "reference_sample_name" not in refcomp.columns or "mge_comparison_set" not in refcomp.columns:
        raise SystemExit(
            "Error: reference_comparison_sets must have columns reference_sample_name, mge_comparison_set"
        )

    if args.row < 0 or args.row >= len(refcomp):
        raise SystemExit(f"Error: --row {args.row} out of range [0, {len(refcomp)-1}]")

    row = refcomp.iloc[args.row]
    reference_sample_name = row["reference_sample_name"]
    mge_comparison_set = row["mge_comparison_set"]
    sample_ids = [s.strip() for s in str(mge_comparison_set).split(",") if s.strip()]

    meta = load_metadata(metadata_path)
    accession_to_sample = meta.set_index("sample_accession")["Sample"].to_dict()

    out_path = args.out if args.out.is_absolute() else wd / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows_out = []
    for sample_id in sample_ids:
        sample_name = accession_to_sample.get(sample_id)
        if sample_name is None or (isinstance(sample_name, float) and pd.isna(sample_name)):
            print(f"Warning: no Sample for sample_accession {sample_id}, skipping", file=sys.stderr)
            continue
        sample_name = str(sample_name).strip()

        contigs_path = discover_assembly_path(assemblies_dir, sample_name)
        if contigs_path is None:
            print(f"Warning: no assembly found for Sample {sample_name}, skipping", file=sys.stderr)
            continue

        data_dir = fastq_dir / sample_id
        rows_out.append({
            "data_dir": str(data_dir),
            "sample_id": sample_id,
            "sample_name": sample_name,
            "gff": ".",
            "contigs": str(contigs_path),
        })

    if not rows_out:
        raise SystemExit("Error: no valid rows to write (check metadata and assemblies)")

    out_df = pd.DataFrame(rows_out)
    out_df.to_csv(out_path, sep="\t", index=False, header=True)
    print(f"Wrote {len(rows_out)} rows to {out_path}")
    print(f"Reference for this run: {reference_sample_name}")


if __name__ == "__main__":
    main()
