# Streaming FASTQ for whole-SL runs (`stream_fastq`)

## Why

A full sublineage (e.g. SL39 ≈ 885 samples) at ~0.9 GB of reads per sample is
~750 GB if every FASTQ is downloaded up front and kept — it will blow the RDS
quota. The default (validation) mode expects reads already on disk and never
deletes them. `stream_fastq` mode downloads each sample's reads on demand and
deletes them as soon as they are no longer needed.

## How it works

FASTQ is only consumed by the `bwa` rule (read → reference alignment). After
`bwa`, every downstream step (`find`, `pair`, `inferseq*`, `make_database`,
`clusterseq`, `genotype`, `summarize`, `makefasta`) works off the BAM and small
per-sample TSVs — never the FASTQ. So:

- A `download_fastq` rule fetches one sample's reads (via `fastq-dl` in its own
  micromamba env) into `<wd>/00.fastq/<accession>_{1,2}.fastq.gz`.
- Those outputs are declared `temp()`, and `bwa` is their only consumer, so
  Snakemake **deletes each sample's FASTQ the moment `bwa` finishes**.
- With `-j 1` this is strictly serial: download → bwa → delete → next sample,
  so peak FASTQ disk ≈ **one sample**. With `-j N`, ≈ N samples in flight.

Aggregation is unchanged: `make_database`/`clusterseq`/`genotype` still gather
all per-sample inferseq/pair TSVs into one per-reference result.

## Enabling it

In `config/config.yaml` set `stream_fastq: true` (or pass `--stream-fastq` to
`run_pipeline.py`), and set `validate_n: -1` to run the whole sublineage. When
on, `run_pipeline.py` skips the bulk download phase and the pre-run FASTQ
existence check; `generate_mgefinder_dataset.py` points `fastq1/fastq2` at the
`<wd>/00.fastq/...` temp paths instead of requiring reads on disk.

## Caveat: compute-node internet

`download_fastq` runs `fastq-dl` inside the SLURM job, so the **compute node
must have outbound internet**. If CSD3 batch nodes cannot reach ENA, options
are: (a) run the whole pipeline somewhere with internet (e.g. locally — it is
CPU-only), or (b) pre-stage reads in small batches on the login node and run
with `stream_fastq: false`, deleting each batch after its job. Test with a
2-sample `stream_fastq` job before committing to the full SL39.

## Quick check (dry run, isolated)

```bash
micromamba run -n snakemake python3 src/generate_mgefinder_dataset.py \
  --config config/config.yaml --sublineage SL39 --limit 2 --stream-fastq \
  --out /tmp/st/mgefinder_dataset.txt --ref-out /tmp/st/ref.txt
# build a cfg with wd=/tmp/st, stream_fastq: true, genomes + reference_assembly_path,
# then: micromamba run -n mgefinder_env snakemake -s workflow/Snakefile \
#   --configfile /tmp/st/cfg.yaml --directory /tmp/st -j1 -n all
# Expect download_fastq jobs + "Would remove temporary output .../00.fastq/*_1.fastq.gz"
```
