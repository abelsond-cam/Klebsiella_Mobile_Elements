#!/usr/bin/env python3
"""Shared per-SL read-pool PRODUCER: parallel, retried, chunked FASTQ download.

Decouples downloads from the MGEfinder compute DAG. Resolves the SL's usable
accessions via generate_mgefinder_dataset.py --print-accessions (single source
of truth — already filtered for read-less / missing-assembly / drop-file), then
downloads them in chunks with a thread pool and exponential-backoff retry,
into a pool shared by every reference-level consumer of that SL:

  <pool>/reads/<acc>_{1,2}.fastq.gz   <pool>/ready/<acc>.ok
  <pool>/chunks/chunk_<k>.txt         <pool>/chunks/chunk_<k>.READY
  <pool>/dropped_samples.tsv          <pool>/COMPLETE

Run INSIDE the fastq-dl micromamba env (the SLURM script activates it ONCE, so
`fastq-dl` is on PATH and there is no per-sample micromamba lock storm). The
accession listing shells out to the snakemake env (needs pandas/yaml).

Idempotent: skips accessions already in ready/ or dropped_samples.tsv, so a
requeued producer resumes. Permanently-failed accessions are dropped with a
loud warning; the SL run still completes (consumers exclude them via the same
drop file at aggregation).
"""

import argparse
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

_drop_lock = threading.Lock()


def list_accessions(config: Path, sublineage: str) -> list:
    """Usable accessions for the SL (post read-less/missing-assembly filters)."""
    cmd = ["micromamba", "run", "-n", "snakemake", "python3",
           "src/generate_mgefinder_dataset.py", "--config", str(config),
           "--sublineage", sublineage, "--limit", "-1", "--print-accessions"]
    out = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    accs = [ln.split("\t")[0].strip() for ln in out.splitlines() if ln.strip()]
    if not accs:
        raise SystemExit(f"Error: no usable accessions for sublineage {sublineage}")
    return accs


def chunks(seq: list, size: int) -> list:
    return [seq[i:i + size] for i in range(0, len(seq), size)]


def already_dropped(pool: Path) -> set:
    f = pool / "dropped_samples.tsv"
    if not f.exists():
        return set()
    return {ln.split("\t")[0].strip() for ln in f.read_text().splitlines()
            if ln.strip() and not ln.startswith("sample_accession")}


def record_drop(pool: Path, acc: str, attempts: int, err: str) -> None:
    line = f"{acc}\t{attempts}\t{err}\n"
    with _drop_lock:
        f = pool / "dropped_samples.tsv"
        if not f.exists():
            f.write_text("sample_accession\tattempts\tlast_error\n")
        with open(f, "a") as fh:
            fh.write(line)
    print(f">>> WARNING: permanently dropping {acc} after {attempts} "
          f"attempt(s): {err}", file=sys.stderr, flush=True)


def fetch_one(acc: str, pool: Path, max_retries: int, backoff_base: int) -> str:
    """Return 'ok' | 'skip' | 'drop'. Downloads acc into the pool with retry."""
    ready = pool / "ready" / f"{acc}.ok"
    r1 = pool / "reads" / f"{acc}_1.fastq.gz"
    r2 = pool / "reads" / f"{acc}_2.fastq.gz"
    if ready.exists() and r1.exists() and r2.exists():
        return "skip"
    tmp = pool / "reads" / f".dl_{acc}"
    last_err = ""
    for attempt in range(1, max_retries + 2):  # max_retries retries after try 1
        if tmp.exists():
            shutil.rmtree(tmp, ignore_errors=True)
        tmp.mkdir(parents=True, exist_ok=True)
        try:
            proc = subprocess.run(
                ["fastq-dl", "--accession", acc, "--outdir", str(tmp)],
                capture_output=True, text=True)
            if proc.returncode == 0:
                f1 = sorted(list(tmp.glob("*_1.fastq.gz")) + list(tmp.glob("*_1.fq.gz")))
                f2 = sorted(list(tmp.glob("*_2.fastq.gz")) + list(tmp.glob("*_2.fq.gz")))
                if f1 and f2:
                    shutil.move(str(f1[0]), str(r1))
                    shutil.move(str(f2[0]), str(r2))
                    shutil.rmtree(tmp, ignore_errors=True)
                    ready.write_text("")
                    return "ok"
                last_err = "no paired FASTQ produced"
            else:
                last_err = (proc.stderr or proc.stdout or "fastq-dl nonzero").strip().splitlines()[-1:]
                last_err = last_err[0] if last_err else "fastq-dl nonzero"
        except Exception as e:  # noqa: BLE001
            last_err = repr(e)
        if attempt <= max_retries:
            time.sleep(backoff_base * (2 ** (attempt - 1)))
    shutil.rmtree(tmp, ignore_errors=True)
    record_drop(pool, acc, max_retries + 1, last_err)
    return "drop"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", type=Path, default=Path("config/config.yaml"))
    ap.add_argument("--sublineage", required=True)
    ap.add_argument("--pool-dir", type=Path, required=True,
                    help="Per-SL pool dir, e.g. .../_readpool/<SL>")
    ap.add_argument("--chunk-size", type=int, default=100)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--max-retries", type=int, default=4)
    ap.add_argument("--backoff-base", type=int, default=20,
                    help="Seconds; sleep = backoff_base * 2**(attempt-1)")
    ap.add_argument("--print-chunks", action="store_true",
                    help="Write chunk files + print partition, then exit (no downloads)")
    args = ap.parse_args()

    pool = args.pool_dir
    for sub in ("reads", "ready", "chunks", "consumed"):
        (pool / sub).mkdir(parents=True, exist_ok=True)

    accs = list_accessions(args.config, args.sublineage)
    parts = chunks(accs, args.chunk_size)
    for k, part in enumerate(parts):
        (pool / "chunks" / f"chunk_{k}.txt").write_text("\n".join(part) + "\n")
    print(f">>> {args.sublineage}: {len(accs)} usable accessions in "
          f"{len(parts)} chunk(s) of <= {args.chunk_size}", flush=True)
    if args.print_chunks:
        return

    dropped = already_dropped(pool)
    for k, part in enumerate(parts):
        ready_flag = pool / "chunks" / f"chunk_{k}.READY"
        todo = [a for a in part if a not in dropped]
        n_ok = n_skip = n_drop = 0
        if todo:
            with ThreadPoolExecutor(max_workers=args.workers) as ex:
                for res in ex.map(lambda a: fetch_one(
                        a, pool, args.max_retries, args.backoff_base), todo):
                    n_ok += res == "ok"
                    n_skip += res == "skip"
                    n_drop += res == "drop"
        ready_flag.write_text("")
        print(f">>> chunk {k}/{len(parts) - 1}: ok={n_ok} skip={n_skip} "
              f"drop={n_drop} (READY)", flush=True)

    (pool / "COMPLETE").write_text("")
    n_drop_total = len(already_dropped(pool))
    print(f">>> {args.sublineage} producer COMPLETE: {len(accs)} accessions, "
          f"{n_drop_total} dropped (see {pool / 'dropped_samples.tsv'})",
          flush=True)


if __name__ == "__main__":
    main()
