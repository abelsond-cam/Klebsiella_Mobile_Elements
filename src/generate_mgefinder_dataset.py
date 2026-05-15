#!/usr/bin/env python3
"""Generate mgefinder_dataset.txt, selecting samples by metadata Sublineage.

Primary (proper) mode: select samples where metadata Sublineage == <sublineage>,
optionally limited to the first N (validate_n / --limit). Assembly and GFF paths
are read from the metadata `assembly_file` / `gff_file` columns, resolved against
`metadata_path_base` (the RDS project root = parent of data_dir). The reference
genome assembly is resolved the same way from its own metadata row.

Legacy mode (back-compat): if --row is given without --sublineage (and config has
no `sublineage`), fall back to selecting from reference_comparison_sets.tsv.

Writes a 7-column TSV: data_dir, sample_id, sample_name, gff, contigs, fastq1, fastq2.

Usage:
  python src/generate_mgefinder_dataset.py --config config/config.yaml --out mgefinder_dataset.txt
  python src/generate_mgefinder_dataset.py --config config/config.yaml --sublineage SL39 --limit 1 --out ...
  python src/generate_mgefinder_dataset.py --config config/config.yaml --print-accessions   # for the downloader
"""

import argparse
import sys
from pathlib import Path
from typing import Optional

import pandas as pd


def discover_fastq_pair(fastq_dir: Path) -> Optional[tuple[Path, Path]]:
    """Find a FASTQ pair in fastq_dir. Looks for *_1.fastq.gz / *_2.fastq.gz (or .fq.gz).
    Returns (path1, path2) for the first pair found, or None if no pair found.
    Does not assume filenames (e.g. DRR* vs SAM*)."""
    for suffix in ("_1.fastq.gz", "_1.fq.gz"):
        for f1 in sorted(fastq_dir.glob(f"*{suffix}")):
            f2_name = f1.name.replace("_1.", "_2.", 1)
            f2 = f1.parent / f2_name
            if f2.exists():
                return (f1, f2)
    return None


def load_metadata(metadata_path: Path) -> pd.DataFrame:
    """Load metadata TSV; require the columns we depend on."""
    df = pd.read_csv(metadata_path, sep="\t", low_memory=False)
    required = ("Sample", "sample_accession", "Sublineage", "assembly_file", "gff_file")
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise SystemExit(f"Error: metadata missing column(s): {', '.join(missing)}")
    return df.astype({"sample_accession": str})


def resolve_under_base(base: Path, rel: object) -> Optional[Path]:
    """Resolve a metadata path column value against metadata_path_base. None if blank."""
    if rel is None or (isinstance(rel, float) and pd.isna(rel)):
        return None
    rel = str(rel).strip()
    if not rel or rel.lower() == "nan":
        return None
    return base / rel


def select_samples(meta: pd.DataFrame, config: dict, args) -> pd.DataFrame:
    """Return the metadata rows to run, in selection order.

    Sublineage mode (proper) when --sublineage / config['sublineage'] is set;
    otherwise legacy reference_comparison_sets row mode.
    """
    sublineage = args.sublineage or config.get("sublineage")
    if sublineage and not args.legacy_row:
        sel = meta[meta["Sublineage"].astype(str) == str(sublineage)].copy()
        if sel.empty:
            raise SystemExit(f"Error: no metadata rows with Sublineage == {sublineage}")
        limit = args.limit if args.limit is not None else config.get("validate_n", -1)
        if limit is not None and int(limit) > 0:
            sel = sel.head(int(limit))
        print(f">>> Sublineage {sublineage}: selected {len(sel)} sample(s) "
              f"(limit={limit})", file=sys.stderr)
        return sel

    # Legacy: reference_comparison_sets.tsv row
    data_dir = Path(config.get("data_dir", config["wd"])).resolve()
    refcomp_path = data_dir / config["reference_comparison_sets"]
    if not refcomp_path.exists():
        raise SystemExit(f"Error: reference_comparison_sets not found: {refcomp_path}")
    refcomp = pd.read_csv(refcomp_path, sep="\t")
    if args.row < 0 or args.row >= len(refcomp):
        raise SystemExit(f"Error: --row {args.row} out of range [0, {len(refcomp)-1}]")
    row = refcomp.iloc[args.row]
    sample_ids = [s.strip() for s in str(row["mge_comparison_set"]).split(",") if s.strip()]
    sel = meta[meta["sample_accession"].isin(sample_ids)].copy()
    # preserve the comparison-set order
    sel["__ord"] = sel["sample_accession"].map({a: i for i, a in enumerate(sample_ids)})
    sel = sel.sort_values("__ord").drop(columns="__ord")
    print(f">>> Legacy row {args.row}: {len(sel)}/{len(sample_ids)} samples found in metadata",
          file=sys.stderr)
    return sel


def get_reference_name(config: dict, args, meta: pd.DataFrame) -> str:
    """Reference genome name: config/CLI in sublineage mode, refcomp row in legacy mode."""
    ref = args.reference_sample_name or config.get("reference_sample_name")
    if ref:
        return str(ref)
    # No explicit reference. Only consult the legacy reference_comparison_sets
    # TSV when legacy row mode is explicitly requested — never silently default.
    if not args.legacy_row:
        raise SystemExit(
            "Error: no reference_sample_name (config or --reference-sample-name) "
            "and not in legacy --legacy-row mode. Refusing to guess a reference.")
    data_dir = Path(config.get("data_dir", config["wd"])).resolve()
    refcomp = pd.read_csv(data_dir / config["reference_comparison_sets"], sep="\t")
    return str(refcomp.iloc[args.row]["reference_sample_name"])


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, default=Path("config/config.yaml"))
    parser.add_argument("--wd", type=Path, default=None,
                        help="Override config 'wd' (must match run_pipeline --output-dir "
                             "so stream-fastq paths point at the run's 00.fastq).")
    parser.add_argument("--out", type=Path, default=Path("mgefinder_dataset.txt"))
    parser.add_argument("--row", type=int, default=0, metavar="N",
                        help="Legacy reference_comparison_sets row index")
    parser.add_argument("--legacy-row", action="store_true",
                        help="Force legacy reference_comparison_sets row mode")
    parser.add_argument("--sublineage", default=None,
                        help="Metadata Sublineage value to select (overrides config)")
    parser.add_argument("--limit", type=int, default=None,
                        help="Limit to first N selected samples (overrides config validate_n)")
    parser.add_argument("--reference-sample-name", dest="reference_sample_name",
                        default=None, help="Reference genome name (overrides config)")
    parser.add_argument("--ref-out", type=Path, default=None,
                        help="Write resolved reference assembly path here (default: <out dir>/mgefinder_reference_path.txt)")
    parser.add_argument("--stream-fastq", action="store_true",
                        help="Don't require FASTQ on disk; point fastq1/fastq2 at "
                             "wd/00.fastq/<id>_{1,2}.fastq.gz (produced+deleted by the "
                             "download_fastq Snakemake rule). Overrides config stream_fastq.")
    parser.add_argument("--print-accessions", action="store_true",
                        help="Print 'sample_accession<TAB>run_accession_used' for selected samples and exit")
    args = parser.parse_args()

    try:
        import yaml
    except ImportError:
        raise SystemExit("Error: PyYAML required. pip install pyyaml")

    if not args.config.exists():
        raise SystemExit(f"Error: config not found: {args.config}")
    with open(args.config) as f:
        config = yaml.safe_load(f)

    if args.wd is not None:
        config["wd"] = str(args.wd)
    data_dir = Path(config.get("data_dir", config["wd"])).resolve()
    wd = Path(config["wd"]).resolve()
    fastq_dir = data_dir / config["fastq_dir"]
    metadata_path = data_dir / config["metadata_file"]
    base = Path(config.get("metadata_path_base", str(data_dir.parent))).resolve()
    stream_fastq = args.stream_fastq or bool(config.get("stream_fastq", False))

    if not metadata_path.exists():
        raise SystemExit(f"Error: metadata not found: {metadata_path}")

    meta = load_metadata(metadata_path)
    selected = select_samples(meta, config, args)

    if args.print_accessions:
        run_col = "run_accession_used" if "run_accession_used" in selected.columns else None
        for _, r in selected.iterrows():
            run = str(r[run_col]).strip() if run_col else ""
            print(f"{r['sample_accession']}\t{run}")
        return

    # Resolve reference assembly path from its own metadata row
    ref_name = get_reference_name(config, args, meta)
    ref_rows = meta[meta["Sample"].astype(str) == ref_name]
    if ref_rows.empty:
        ref_rows = meta[meta["sample_accession"] == ref_name]
    if ref_rows.empty:
        raise SystemExit(f"Error: reference '{ref_name}' not found in metadata "
                          f"(by Sample or sample_accession)")
    ref_assembly = resolve_under_base(base, ref_rows.iloc[0]["assembly_file"])
    if ref_assembly is None or not ref_assembly.exists():
        raise SystemExit(f"Error: reference assembly not found on disk: {ref_assembly}")

    out_path = args.out if args.out.is_absolute() else wd / args.out
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ref_out = args.ref_out or (out_path.parent / "mgefinder_reference_path.txt")

    rows_out = []
    for _, r in selected.iterrows():
        sample_id = str(r["sample_accession"]).strip()
        sample_name = str(r["Sample"]).strip()

        # Skip read-less entries: RefSeq/assembly-only genomes (is_refseq=True,
        # GCF_* accession, no run_accession_used) have no short reads to align.
        # Including them breaks the bwa DAG (nothing produces their FASTQ).
        is_refseq = str(r.get("is_refseq", "")).strip().lower() in ("true", "1", "yes")
        run_used = str(r.get("run_accession_used", "")).strip()
        no_run = run_used == "" or run_used.lower() in ("nan", "na", "none")
        if is_refseq or no_run or sample_id.startswith("GCF_"):
            print(f"Skipping read-less sample {sample_id} "
                  f"(is_refseq={is_refseq}, run_accession_used={run_used or 'n/a'})",
                  file=sys.stderr)
            continue

        contigs = resolve_under_base(base, r["assembly_file"])
        if contigs is None or not contigs.exists():
            print(f"Warning: assembly missing for {sample_id}: {contigs}, skipping",
                  file=sys.stderr)
            continue
        gff = resolve_under_base(base, r["gff_file"])
        gff_str = str(gff) if gff is not None and gff.exists() else "."

        sample_fastq_dir = fastq_dir / sample_id
        if stream_fastq:
            # download_fastq rule will create (and Snakemake will temp()-delete) these.
            fastq1 = wd / "00.fastq" / f"{sample_id}_1.fastq.gz"
            fastq2 = wd / "00.fastq" / f"{sample_id}_2.fastq.gz"
        else:
            fastq_pair = discover_fastq_pair(sample_fastq_dir)
            if fastq_pair is None:
                run_used = str(r.get("run_accession_used", "")).strip()
                print(f"Warning: no FASTQ pair in {sample_fastq_dir} "
                      f"(run_accession_used={run_used or 'n/a'}); download then re-run. Skipping.",
                      file=sys.stderr)
                continue
            fastq1, fastq2 = fastq_pair
        rows_out.append({
            "data_dir": str(sample_fastq_dir),
            "sample_id": sample_id,
            "sample_name": sample_name,
            "gff": gff_str,
            "contigs": str(contigs),
            "fastq1": str(fastq1),
            "fastq2": str(fastq2),
        })

    if not rows_out:
        raise SystemExit("Error: no valid rows to write (check metadata paths and FASTQ downloads)")

    pd.DataFrame(rows_out).to_csv(out_path, sep="\t", index=False, header=True)
    ref_out.write_text(str(ref_assembly) + "\n")
    print(f"Wrote {len(rows_out)} rows to {out_path}")
    print(f"Reference for this run: {ref_name}")
    print(f"Reference assembly: {ref_assembly}  (recorded in {ref_out})")


if __name__ == "__main__":
    main()
