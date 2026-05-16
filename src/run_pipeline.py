#!/usr/bin/env python3
"""Pipeline runner: Snakemake from snakemake env; MGEfinder tools via conda in Snakefile."""

import argparse
import csv
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


def _stats(xs: list):
    """(mean, min, max) or NaNs if empty."""
    if not xs:
        nan = float("nan")
        return (nan, nan, nan)
    return (sum(xs) / len(xs), min(xs), max(xs))


def resolve_sl_reference(csv_path: Path, sample_accessions: list,
                         sublineage: str, level: str = "f"):
    """Per-SL reference from the panaroo best-reference CSV at a given level.

    `level` selects the ref_<level>/shared_<level> column pair: 'f' = ref_f
    (MGH 78578 universal basal), 'd' = ref_d (per-SL level-d). Joins CSV
    'Sample' == metadata sample_accession (the CSV 'run' column is unreliable
    for split SLs, e.g. SL258_part_*, and must NOT be used as the SL key).
    Enforces a single distinct ref_<level> across the matched samples
    (fail-fast — no silent modal pick). shared_d/shared_f vary per sample and
    are always summarised (mean, min, max) so the d-vs-f delta can be printed
    regardless of the chosen level.

    Returns (ref, shared_chosen_stats, shared_d_stats, shared_f_stats,
    n_matched, n_total), or (None, None, None, None, 0, n_total) if no SL
    sample is present in the CSV.
    """
    if not csv_path.exists():
        raise SystemExit(f"Error: best_reference_csv not found: {csv_path}")
    ref_col = f"ref_{level}"
    wanted = set(sample_accessions)
    ref_counts: dict = {}
    shared_d: list = []
    shared_f: list = []
    shared_chosen: list = []
    with open(csv_path, newline="") as fh:
        reader = csv.DictReader(fh)
        for col in ("Sample", ref_col, "shared_d", "shared_f",
                    f"shared_{level}"):
            if reader.fieldnames is None or col not in reader.fieldnames:
                raise SystemExit(f"Error: {csv_path} missing required column '{col}'")
        for row in reader:
            if row["Sample"] not in wanted:
                continue
            rv = (row[ref_col] or "").strip()
            if not rv:
                continue
            ref_counts[rv] = ref_counts.get(rv, 0) + 1
            for raw, acc in ((row.get("shared_d"), shared_d),
                             (row.get("shared_f"), shared_f),
                             (row.get(f"shared_{level}"), shared_chosen)):
                try:
                    acc.append(float(raw))
                except (TypeError, ValueError):
                    pass
    n_total = len(wanted)
    n_matched = sum(ref_counts.values())
    if n_matched == 0:
        return None, None, None, None, 0, n_total
    if len(ref_counts) > 1:
        dist = ", ".join(f"{k}={v}" for k, v in
                         sorted(ref_counts.items(), key=lambda kv: -kv[1]))
        raise SystemExit(
            f"Error: QC failed — sublineage {sublineage} has "
            f"{len(ref_counts)} distinct level-{level} references ({ref_col}) "
            f"among {n_matched} matched samples: {dist}. Expected exactly one; "
            f"refusing to pick one silently.")
    ref = next(iter(ref_counts))
    return (ref, _stats(shared_chosen), _stats(shared_d), _stats(shared_f),
            n_matched, n_total)


def get_reference_info(config: dict, args, config_path: Path) -> tuple:
    """Resolve (ref_name, comparison_ids) for the run.

    Sublineage mode (proper): ref_name from config['reference_sample_name'];
    sample accessions obtained from generate_mgefinder_dataset.py --print-accessions
    (single source of truth for selection). Legacy mode: reference_comparison_sets row.
    """
    sublineage = args.sublineage or config.get("sublineage")
    test_n = args.test_n

    if sublineage:
        ref_name = config.get("reference_sample_name") or args.reference_sample_name
        if not ref_name:
            raise SystemExit("Error: sublineage mode requires config 'reference_sample_name'")
        cmd = ["python3", "src/generate_mgefinder_dataset.py",
               "--config", str(config_path), "--sublineage", sublineage,
               "--print-accessions"]
        if test_n is not None and test_n > 0:
            cmd += ["--limit", str(test_n)]
        out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
        comparison_ids = [ln.split("\t")[0].strip()
                          for ln in out.splitlines() if ln.strip()]

        # Per-SL reference from the best-reference CSV at the chosen level
        # (default 'f' = ref_f = MGH 78578 basal; 'd' = ref_d level-d). The
        # --reference-level flag makes each job's reference DETERMINISTIC
        # regardless of later config/code changes. Resolved over the FULL SL
        # (independent of --test-n/validate_n) so the ref uniqueness QC and
        # the shared-content summary are meaningful.
        level = args.reference_level or config.get("reference_level", "f")
        csv_path = config.get("best_reference_csv")
        if csv_path:
            full_cmd = ["python3", "src/generate_mgefinder_dataset.py",
                        "--config", str(config_path), "--sublineage", sublineage,
                        "--print-accessions", "--limit", "-1"]
            full_out = subprocess.run(full_cmd, capture_output=True, text=True,
                                      check=True).stdout
            full_ids = [ln.split("\t")[0].strip()
                        for ln in full_out.splitlines() if ln.strip()]
            ref, sc, sd, sf, nm, nt = resolve_sl_reference(
                Path(csv_path), full_ids, sublineage, level)
            if ref is None:
                print(f">>> Warning: no {sublineage} samples in {csv_path}; "
                      f"falling back to basal reference {ref_name}",
                      file=sys.stderr)
            else:
                delta = sd[0] - sf[0]
                print(f">>> reference-level {level}: {sublineage} -> {ref} "
                      f"(matched {nm}/{nt}); shared_{level} mean {sc[0]:.0f} "
                      f"(range {sc[1]:.0f}-{sc[2]:.0f}). "
                      f"level-d vs level-f(MGH78578): shared_d mean {sd[0]:.0f} "
                      f"vs shared_f mean {sf[0]:.0f} (delta {delta:+.0f}"
                      f"{'; this run uses level-f = MGH78578 baseline' if level == 'f' else ''})")
                ref_name = ref

        print(f">>> Sublineage {sublineage}: {len(comparison_ids)} sample(s); "
              f"reference={ref_name}")
        return ref_name, comparison_ids

    # Legacy reference_comparison_sets row mode
    data_dir = Path(config.get("data_dir", config["wd"])).resolve()
    refcomp_path = data_dir / config["reference_comparison_sets"]
    if not refcomp_path.exists():
        raise SystemExit(f"Error: reference_comparison_sets not found: {refcomp_path}")
    df = pd.read_csv(refcomp_path, sep="\t")
    if args.row >= len(df):
        raise SystemExit(f"Error: row {args.row} out of range")
    row_data = df.iloc[args.row]
    ref_name = row_data["reference_sample_name"]
    comparison_ids = [s.strip() for s in str(row_data["mge_comparison_set"]).split(",") if s.strip()]
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


def _dataset_acc_to_name(dataset_path: Path) -> dict:
    """sample_accession (col2) -> sample_name (col3) from the 7-col dataset TSV."""
    m = {}
    with open(dataset_path) as f:
        next(f, None)
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if len(p) >= 3:
                m[p[1].strip()] = p[2].strip()
    return m


def chunk_targets(pool: Path, chunk_k: int, dataset_path: Path,
                   genome: str, mgefinder_dir: str, wd: Path,
                   dropped: set) -> tuple:
    """(target_paths, accessions) for chunk_k's 03.inferseq_* leaf outputs.

    Maps pool accessions -> sample_name via the dataset TSV; skips dropped or
    accessions absent from the dataset (avoids a bwa KeyError aborting the run).
    """
    chunk_file = pool / "chunks" / f"chunk_{chunk_k}.txt"
    if not chunk_file.exists():
        raise SystemExit(f"Error: chunk file not found: {chunk_file}")
    acc2name = _dataset_acc_to_name(dataset_path)
    targets, accs = [], []
    for acc in (l.strip() for l in chunk_file.read_text().splitlines() if l.strip()):
        if acc in dropped:
            continue
        name = acc2name.get(acc)
        if name is None:
            print(f">>> chunk {chunk_k}: {acc} not in dataset (filtered/dropped); "
                  f"skipping", file=sys.stderr)
            continue
        base = f"{mgefinder_dir}/{genome}/{name}/03.inferseq_{{}}.{name}.{genome}.tsv"
        for kind in ("assembly", "reference", "overlap"):
            targets.append(str(wd / base.format(kind)))
        accs.append((acc, name))
    return targets, accs


def mark_consumed(pool: Path, tag: str, accs: list, genome: str,
                  mgefinder_dir: str, wd: Path) -> None:
    """touch consumed/<tag>/<acc>.done for accs whose 3 inferseq TSVs now exist."""
    cdir = pool / "consumed" / tag
    cdir.mkdir(parents=True, exist_ok=True)
    done = 0
    for acc, name in accs:
        ok = all((wd / mgefinder_dir / genome / name /
                  f"03.inferseq_{k}.{name}.{genome}.tsv").exists()
                 for k in ("assembly", "reference", "overlap"))
        if ok:
            (cdir / f"{acc}.done").write_text("")
            done += 1
    print(f">>> consumed[{tag}]: {done}/{len(accs)} chunk samples produced "
          f"03.inferseq TSVs", file=sys.stderr)


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

    ref_path = merged.get("reference_assembly_path")
    if not ref_path or not Path(ref_path).exists():
        print(f"Error: reference assembly not found: {ref_path} "
              f"(resolved from metadata for reference '{ref_name}')", file=sys.stderr)
        sys.exit(1)

    if not Path(contigs).exists():
        print(f"Error: sample assembly not found: {contigs}", file=sys.stderr)
        sys.exit(1)

    if not has_fastq_cols:
        print("Error: dataset has no fastq1/fastq2 columns; re-run dataset generation to get 7-column format.", file=sys.stderr)
        sys.exit(1)
    if merged.get("stream_fastq"):
        # FASTQ are produced on demand by the download_fastq rule; not on disk yet.
        print(">>> stream_fastq: skipping FASTQ existence check (downloaded per-sample at run time)")
        return
    if merged.get("read_pool_mode"):
        # Pool reads arrive per chunk (producer-owned); the first dataset row
        # may be in a not-yet-downloaded chunk. Readiness is enforced by chunk
        # sentinels, not an up-front check.
        print(">>> read-pool: skipping FASTQ existence check (per-chunk readiness)")
        return
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
    parser.add_argument("--output-dir", type=Path, default=None,
                        help="Override config 'wd': all pipeline outputs (00.bam, "
                             "00.fastq, mgefinder, database, results) go here. "
                             "Use to isolate a run (e.g. .../mgefinder/SL39).")
    parser.add_argument("--row", type=int, default=0, help="Legacy: row from reference_comparison_sets")
    parser.add_argument("--sublineage", default=None, help="Metadata Sublineage to select (overrides config)")
    parser.add_argument("--reference-sample-name", dest="reference_sample_name", default=None,
                        help="Reference genome name (overrides config)")
    parser.add_argument("--fetch-only", action="store_true", help="Only fetch data, don't run pipeline")
    parser.add_argument("--pipeline-only", action="store_true", help="Only run pipeline, skip fetch")
    parser.add_argument("--verbose", action="store_true", help="Use verbose Snakemake output")
    parser.add_argument("--jobs", "-j", type=int, default=1, help="Snakemake jobs")
    parser.add_argument("--test-n", type=int, help="TEST MODE: Limit to first N comparison samples")
    parser.add_argument("--dry-run", action="store_true", help="Snakemake dry run (-n)")
    parser.add_argument("--skip-download", action="store_true", help="Skip FASTQ download step")
    parser.add_argument("--stream-fastq", action="store_true",
                        help="Per-sample download via Snakemake, deleted after bwa "
                             "(for whole-SL runs). Overrides config stream_fastq.")
    parser.add_argument("--reference-level", choices=["f", "d"], default=None,
                        help="Which best_reference_per_sample.csv level to use as "
                             "the per-SL reference: 'f' = ref_f (MGH 78578 basal, "
                             "default), 'd' = ref_d (per-SL level-d). Overrides "
                             "config reference_level. Makes the reference "
                             "deterministic per job.")
    parser.add_argument("--read-pool", type=Path, default=None,
                        help="Consume reads from a shared per-SL pool dir "
                             "(producer-owned). Forces non-stream, skips PHASE 2 "
                             "download; dataset fastq cols point at the pool.")
    parser.add_argument("--chunk", type=int, default=None,
                        help="Read-pool mode: process ONLY this chunk's accessions "
                             "(explicit per-sample 03.inferseq targets, no cohort "
                             "aggregation). Requires --read-pool + --consumer-tag.")
    parser.add_argument("--consumer-tag", default=None,
                        help="Read-pool consumer id (e.g. ref_f / ref_d) for the "
                             "pool consumed/<tag>/<acc>.done sentinels.")
    parser.add_argument("--aggregate-only", action="store_true",
                        help="Read-pool mode: regenerate dataset to the successful "
                             "set (--drop-file) and run ONLY the cohort aggregation "
                             "(all per-sample TSVs already exist).")
    parser.add_argument("--drop-file", type=Path, default=None,
                        help="Drop-list (producer dropped_samples.tsv) of permanently "
                             "undownloadable accessions; excluded at aggregation so "
                             "SAMPLES = successful set.")
    args = parser.parse_args()
    
    # Environment names (micromamba; no --use-conda so Snakemake never invokes conda)
    SNAKE_ENV = "snakemake"        # Data prep (validate, dataset generation)
    MGEFINDER_ENV = "mgefinder_env"  # Snakemake + all rule jobs (bwa, bowtie2, mgefinder); you manage this env yourself
    FASTQ_ENV = "fastq-dl"         # FASTQ downloads only
    
    config = load_config(args.config)
    if args.output_dir is not None:
        config["wd"] = str(args.output_dir)
        print(f">>> --output-dir override: wd = {args.output_dir}")
    wd = Path(config["wd"]).resolve()
    wd.mkdir(parents=True, exist_ok=True)

    read_pool = args.read_pool.resolve() if args.read_pool else None
    pool_mode = read_pool is not None
    if pool_mode and args.chunk is not None and not args.consumer_tag:
        raise SystemExit("Error: --chunk requires --consumer-tag (e.g. ref_f)")

    ref_name, comparison_ids = get_reference_info(config, args, args.config)
    
    print(f">>> Pipeline setup: wd={wd}, reference={ref_name}")
    print(f">>> Comparison samples: {len(comparison_ids)} ({', '.join(comparison_ids[:3])}{'...' if len(comparison_ids) > 3 else ''})")
    print(f">>> Data prep: {SNAKE_ENV}; Snakemake (no conda): {MGEFINDER_ENV}")
    
    # PHASE 1: Data preparation (snakemake environment)
    if not args.pipeline_only:
        print("\n" + "="*70)
        print("PHASE 1: DATA PREPARATION")
        print("="*70)
        
        stream_fastq = (False if pool_mode else
                        (args.stream_fastq or bool(config.get("stream_fastq", False))))

        # PHASE 2: FASTQ downloads. Skipped in stream_fastq mode (download_fastq
        # rule fetches on demand) AND in read-pool mode (the external producer
        # owns all downloads into the shared pool).
        if pool_mode:
            print(f"\n>>> read-pool: downloads owned by the producer "
                  f"({read_pool}); skipping in-pipeline download")
        elif stream_fastq:
            print("\n>>> stream_fastq: per-sample download handled by Snakemake "
                  "(reads deleted after bwa); skipping bulk download")
        elif comparison_ids and not args.skip_download:
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
        
        # Generate dataset (sublineage mode; --limit for small validation runs)
        dataset_path = wd / "mgefinder_dataset.txt"
        ref_path_file = wd / "mgefinder_reference_path.txt"
        gen_cmd = [
            "python3", "src/generate_mgefinder_dataset.py",
            "--config", str(args.config),
            "--row", str(args.row),
            "--out", str(dataset_path),
            "--ref-out", str(ref_path_file),
            # ref_name is the authoritative reference (basal, or per-SL ref_d
            # from resolve_sl_reference). Pass it so the generator resolves the
            # SAME genome's assembly — generate's --reference-sample-name wins
            # over config, closing any divergence with merged["genomes"].
            "--reference-sample-name", str(ref_name),
        ]
        if args.output_dir is not None:
            gen_cmd += ["--wd", str(wd)]
        sublineage = args.sublineage or config.get("sublineage")
        if sublineage:
            gen_cmd += ["--sublineage", sublineage]
        if stream_fastq:
            gen_cmd += ["--stream-fastq"]
        if pool_mode:
            gen_cmd += ["--read-pool", str(read_pool)]
        # Chunk runs use the FULL usable cohort (SAMPLES=full; explicit targets
        # drive scope). ONLY the final aggregation prunes to the successful set
        # so expand(sample=SAMPLES) cannot deadlock on a dropped accession.
        if args.aggregate_only and args.drop_file is not None:
            gen_cmd += ["--drop-file", str(args.drop_file)]
        if args.test_n is not None:
            gen_cmd += ["--limit", str(args.test_n)]
        run_with_env(gen_cmd, "Generate dataset", SNAKE_ENV)
        
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
        # --stream-fastq is a CLI flag; propagate the effective value into the
        # merged config so both check_first_inputs and the Snakefile see it
        # (config.yaml may still say stream_fastq: false).
        merged["stream_fastq"] = (False if pool_mode else
                                  bool(args.stream_fastq or config.get("stream_fastq", False)))
        merged["read_pool_mode"] = pool_mode
        ref_path_file = wd / "mgefinder_reference_path.txt"
        if ref_path_file.exists():
            merged["reference_assembly_path"] = ref_path_file.read_text().strip()
            print(f">>> Reference assembly: {merged['reference_assembly_path']}")
        else:
            print(f">>> Warning: {ref_path_file} not found; Snakefile will fall back to "
                  f"assemblies_dir glob for reference '{ref_name}'", file=sys.stderr)
        merged_yaml = Path("config/.mge_merged_config.yaml")
        merged_yaml.parent.mkdir(parents=True, exist_ok=True)
        with open(merged_yaml, "w") as f:
            yaml.dump(merged, f, default_flow_style=False, sort_keys=False)
        print(f">>> Created {merged_yaml} (config + genomes)")
        
        if not args.dry_run:
            check_first_inputs(wd, merged, ref_name)
        
        # Run Snakemake from working environment (single config file so all keys are present)
        config_abs = merged_yaml.resolve()
        # Cap concurrent FASTQ downloads (whole-SL streaming hammered ENA and
        # the mamba lock under -j 76, aborting runs). CPU rules still use all
        # --jobs cores; only download_fastq (resources: dl_slots=1) is throttled
        # to dl_slots. --keep-going so one bad sample doesn't waste the other
        # ~99% of in-progress work (resubmit then resumes via the DAG).
        mgefinder_dir = config.get("mgefinder_dir", "mgefinder")
        snake_base = ["snakemake", "--configfile", str(config_abs),
                      "--directory", str(wd), "-j", str(args.jobs)]
        chunk_accs = None
        if pool_mode and args.chunk is not None and not args.aggregate_only:
            # Per-chunk: explicit 03.inferseq leaf targets only — NEVER `all`,
            # so the expand(sample=SAMPLES) cohort barrier is not in the DAG.
            dropped = set()
            if args.drop_file is not None and args.drop_file.exists():
                dropped = {l.split("\t")[0].strip()
                           for l in args.drop_file.read_text().splitlines()
                           if l.strip() and not l.startswith("sample_accession")}
            targets, chunk_accs = chunk_targets(
                read_pool, args.chunk, dataset_path, ref_name,
                mgefinder_dir, wd, dropped)
            if not targets:
                raise SystemExit(f"Error: chunk {args.chunk} has no resolvable "
                                 f"targets (all dropped/filtered?)")
            snake_cmd = snake_base + ["--rerun-incomplete", "--keep-going"] + targets
        elif pool_mode and args.aggregate_only:
            # All per-sample TSVs already exist; SAMPLES = successful set
            # (dataset regenerated with --drop-file). Run the cohort barrier.
            snake_cmd = snake_base + ["--rerun-incomplete", "all"]
        else:
            # Legacy stream / non-stream whole-run (download_fastq in DAG).
            dl_slots = int(merged.get("download_concurrency", 12))
            snake_cmd = snake_base + ["--resources", f"dl_slots={dl_slots}",
                                      "--keep-going", "all"]

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
            if chunk_accs is not None:
                mark_consumed(read_pool, args.consumer_tag, chunk_accs,
                              ref_name, mgefinder_dir, wd)
                print(f"\n>>> Chunk {args.chunk} ({args.consumer_tag}) done; "
                      f"per-sample TSVs persisted in {wd}/{mgefinder_dir}/")
            else:
                print(f"\n>>> Pipeline completed successfully!")
                print(f"    Results in: {wd}/results/{ref_name}/")


if __name__ == "__main__":
    main()
