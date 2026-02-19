#!/usr/bin/env python3
"""Create symlinks for assembly files: lookup by Sample, symlink name = original filename (Sample + extension)."""

import argparse
import os
import sys
from pathlib import Path

import pandas as pd


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Bulk symlink assembly files: lookup path by Sample, symlink name = original assembly filename"
    )
    parser.add_argument(
        "--assembly-list",
        default="/home/dca36/rds/rds-floto-bacterial-4k08a2yyQLw/david/raw/kpsc_assembly_files.txt",
        help="Path to file listing one assembly path per line (e.g. kpsc_assembly_files.txt)",
    )
    parser.add_argument(
        "--metadata",
        default="/home/dca36/rds/rds-floto-bacterial-4k08a2yyQLw/david/final/metadata_final_curated_slimmed.tsv",
        help="Path to metadata TSV with Sample and sample_accession columns",
    )
    parser.add_argument(
        "--out-dir",
        default="/home/dca36/rds/rds-floto-bacterial-4k08a2yyQLw/david/raw/assemblies",
        help="Directory where symlinks will be created (default: current directory)",
    )
    parser.add_argument(
        "--test-set",
        type=int,
        default=None,
        # "metavar" changes the name shown in the help message for this argument. "N" will appear as a placeholder for the argument value (number of links to create): e.g., "--test-set N".
        metavar="N",
        help="Create only N symlinks then stop (e.g. 100)",
    )
    parser.add_argument(
        "--path-prefix-old",
        default="/home/shgb2/rds",
        help="Path prefix in assembly list to replace (for converting to your mount)",
    )
    parser.add_argument(
        "--path-prefix-new",
        default="/home/dca36/rds",
        help="Path prefix to use for symlink targets (your equivalent path)",
    )
    args = parser.parse_args()

    # Resolve out_dir and create if needed
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {out_dir}")

    # Build path_by_sample from assembly list (convert path prefix for this lead if needed)
    path_by_sample: dict[str, str] = {}
    old_prefix = args.path_prefix_old
    new_prefix = args.path_prefix_new
    with open(args.assembly_list) as f:
        for line in f:
            path = line.strip()
            if not path:
                continue
            path = os.path.abspath(path)
            if old_prefix and new_prefix and path.startswith(old_prefix):
                path = new_prefix + path[len(old_prefix) :]
            basename = os.path.basename(path)
            for ext in (".fa.gz", ".fna.gz", ".fasta.gz"):
                if basename.endswith(ext):
                    sample_id = basename[: -len(ext)]
                    path_by_sample[sample_id] = path
                    break

    # Load metadata
    df = pd.read_csv(args.metadata, sep="\t", low_memory=False)
    for col in ("Sample", "sample_accession"):
        if col not in df.columns:
            print(f"Error: metadata missing column '{col}'", file=sys.stderr)
            sys.exit(1)

    # For debugging: assemblies present in kpsc list but not in metadata
    metadata_samples = {str(s).strip() for s in df["Sample"].dropna()}
    extra_assemblies = [sid for sid in path_by_sample if sid not in metadata_samples]

    created = 0
    skipped_existing = 0
    missing = 0
    matched_rows = 0
    matched_accessions: set[str] = set()
    matched_samples: set[str] = set()
    # For debugging: track which Samples share the same sample_accession
    samples_by_accession: dict[str, list[str]] = {}
    limit = args.test_set

    for _, row in df.iterrows():
        if limit is not None and created >= limit:
            break
        sample = row["Sample"]
        accession = row["sample_accession"]
        if pd.isna(sample) or pd.isna(accession):
            continue
        sample = str(sample).strip()
        accession = str(accession).strip()
        if sample not in path_by_sample:
            missing += 1
            continue

        matched_rows += 1
        matched_accessions.add(accession)
        matched_samples.add(sample)
        samples_by_accession.setdefault(accession, []).append(sample)

        target = path_by_sample[sample]
        link_path = out_dir / os.path.basename(target)
        if os.path.lexists(link_path):
            skipped_existing += 1
            continue
        os.symlink(target, link_path)
        created += 1

    print(f"Created: {created}, skipped (existing): {skipped_existing}, missing from assembly list: {missing}")
    print(f"Matched rows (Sample in assemblies): {matched_rows}")
    print(f"Unique Sample (symlink names): {len(matched_samples)}")

    # Show a few examples where multiple Samples share the same sample_accession
    dup_examples = [(acc, s_list) for acc, s_list in samples_by_accession.items() if len(s_list) > 1]
    if dup_examples:
        print(f"sample_accession values shared by >1 Sample: {len(dup_examples)} total")
        print("Example duplicates (up to 10):")
        for acc, s_list in dup_examples[:10]:
            print(f"  sample_accession={acc}")
            for s in s_list:
                print(f"    Sample={s}")

    # Show up to 10 example assembly entries that are not used (no matching metadata Sample)
    if extra_assemblies:
        print(f"Assemblies present in list but not in metadata Sample: {len(extra_assemblies)} total")
        print("Example (up to 10):")
        for sid in extra_assemblies[:10]:
            print(f"  {sid}\t{path_by_sample[sid]}")


if __name__ == "__main__":
    main()
