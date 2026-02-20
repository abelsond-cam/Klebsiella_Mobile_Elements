#!/usr/bin/env python3
"""Download FASTQ files using fastq-dl for selected genomes.

Accepts a list of sample_accessions via --ids-file (one per line) or --ids,
or falls back to select_genomes(n) from reference_comparison_sets. Writes
to fastq_dir/{sample_id}/{sample_id}_1.fastq.gz and _2.fastq.gz (fastq-dl
typically uses the accession as prefix when downloading to that dir).
"""

import argparse
import subprocess
import sys
from pathlib import Path

# No default FASTQ dir - require explicit --config or --fastq-dir


def check_fastq_exists(outdir: Path, sample_id: str) -> bool:
    """Return True if expected paired FASTQ files exist for sample_id."""
    if not outdir.exists():
        return False
    f1 = outdir / f"{sample_id}_1.fastq.gz"
    f2 = outdir / f"{sample_id}_2.fastq.gz"
    if f1.exists() and f2.exists():
        return True
    # Fallback: any common FASTQ in dir (for different naming from fastq-dl)
    for pattern in ("*.fastq.gz", "*.fq.gz", "*.fastq", "*.fq"):
        if list(outdir.glob(pattern)):
            return True
    return False


def download_fastq(
    accessions: list,
    output_base_dir: Path,
    skip_existing: bool = True,
) -> dict:
    """
    Download FASTQ files using fastq-dl.

    Args:
        accessions: List of SRA/sample accession IDs
        output_base_dir: Base directory; each acc gets output_base_dir/acc/
        skip_existing: If True, skip accessions that already have FASTQ files

    Returns:
        Dict with keys success, failed, skipped (lists of accession strings)
    """
    results = {"success": [], "failed": [], "skipped": []}

    for i, acc in enumerate(accessions, 1):
        print(f"\n[{i}/{len(accessions)}] Processing {acc}...")

        outdir = output_base_dir / acc
        outdir.mkdir(parents=True, exist_ok=True)

        if skip_existing and check_fastq_exists(outdir, acc):
            print(f"Skipping {acc} (FASTQ files already exist)")
            results["skipped"].append(acc)
            continue

        try:
            cmd = ["fastq-dl", "--accession", acc, "--outdir", str(outdir)]
            print(f"Command: {' '.join(cmd)}")
            subprocess.run(cmd, capture_output=True, text=True, check=True)
            print(f"Successfully downloaded {acc}")
            results["success"].append(acc)
        except subprocess.CalledProcessError as e:
            print(f"Failed to download {acc}: {e.stderr.strip()}")
            results["failed"].append(acc)
        except Exception as e:
            print(f"Unexpected error downloading {acc}: {e}")
            results["failed"].append(acc)

    return results


def get_accessions_from_args(args: argparse.Namespace) -> list:
    """Resolve list of sample_accessions from --ids-file, --ids, or select_genomes(--n)."""
    if args.ids_file is not None:
        path = Path(args.ids_file)
        if not path.exists():
            raise SystemExit(f"Error: --ids-file not found: {path}")
        return [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if args.ids:
        return list(args.ids)
    # Fallback: select_genomes(n) from reference_comparison
    from select_genomes_reference_comparison import select_genomes
    references, comparisons = select_genomes(args.n)
    accessions = [r["sample_accession"] for r in references] + [
        c["sample_accession"] for c in comparisons
    ]
    return accessions


def get_fastq_dir(args: argparse.Namespace) -> Path:
    """Resolve fastq output base dir from --config or --fastq-dir. Fail if neither works."""
    if args.fastq_dir is not None:
        return Path(args.fastq_dir).resolve()
    
    if args.config is not None and args.config.exists():
        try:
            import yaml
            with open(args.config) as f:
                cfg = yaml.safe_load(f)
            
            if "fastq_dir" not in cfg:
                raise SystemExit("Error: config missing 'fastq_dir' key")
            
            data_dir = Path(cfg.get("data_dir", cfg["wd"])).resolve()
            return data_dir / cfg["fastq_dir"]
        except Exception as e:
            raise SystemExit(f"Error: could not read fastq_dir from config: {e}")
    
    raise SystemExit("Error: must provide either --config (with fastq_dir) or --fastq-dir")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Download FASTQ files for MGEfinder analysis"
    )
    parser.add_argument(
        "--ids-file",
        type=Path,
        default=None,
        help="Path to file with one sample_accession per line (overrides --n and --ids)",
    )
    parser.add_argument(
        "--ids",
        nargs="+",
        default=[],
        help="Space-separated sample_accessions (overrides --n)",
    )
    parser.add_argument(
        "--n",
        type=int,
        default=10,
        help="Number of comparison genomes per KL subset (used only if neither --ids-file nor --ids)",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=None,
        help="Config YAML path; used to set fastq_dir from wd + fastq_dir",
    )
    parser.add_argument(
        "--fastq-dir",
        type=Path,
        default=None,
        help="Base directory for FASTQ output (overrides config)",
    )
    parser.add_argument(
        "--overwrite-existing",
        action="store_true",
        help="Re-download even if FASTQ files already exist",
    )
    args = parser.parse_args()

    accessions = get_accessions_from_args(args)
    if not accessions:
        raise SystemExit("Error: no accessions to download (provide --ids-file, --ids, or use --n)")

    fastq_dir = get_fastq_dir(args)
    skip_existing = not args.overwrite_existing
    mode_msg = "overwriting existing" if args.overwrite_existing else "skipping existing"
    print(f"===== FASTQ Download ({len(accessions)} accessions, {mode_msg}) =====")
    print(f"Output directory: {fastq_dir}")

    results = download_fastq(accessions, fastq_dir, skip_existing=skip_existing)

    print("\n===== Summary =====")
    print(f"Success: {len(results['success'])}, Skipped: {len(results['skipped'])}, Failed: {len(results['failed'])}")
    if results["failed"]:
        for acc in results["failed"]:
            print(f"  Failed: {acc}")
        sys.exit(1)
    print("All downloads completed successfully.")
