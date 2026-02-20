micro# What Bowtie2 does in this pipeline (BAM → Bowtie2)

## What is Bowtie2?

**Bowtie2** is an aligner: it maps short sequencing reads to a reference genome (or assembly) and outputs alignments (SAM/BAM). Unlike BWA, which we use earlier in this pipeline to align **short reads to the reference**, Bowtie2 here is used only to **build indexes** of the reference and of each sample assembly. Those indexes are then used by **MGEfinder** to recover and extend insertion sequences.

## What Bowtie2 does in this pipeline

1. **Rule `index_genome_bowtie2`**  
   Builds a Bowtie2 index of the **reference genome** (e.g. `00.genome/GCF_000814305.1_ASM81430v1_genomic.fna`).  
   - **Input**: reference FASTA (from `copy_genome`).  
   - **Output**: `.bt2` index files under `00.genome/` (e.g. `*.fna.1.bt2`, `*.fna.2.bt2`, …).  
   - **Command**: `bowtie2-build <reference.fna> <reference.fna>`.

2. **Rule `index_assembly`**  
   Builds a Bowtie2 index of **each sample assembly** (the sample’s contigs, e.g. `00.assembly/SAMD00052619.fna`).  
   - **Input**: sample assembly FASTA (from `copy_assembly`).  
   - **Output**: `.bt2` index files under `00.assembly/` (e.g. `SAMD00052619.fna.1.bt2`, …).  
   - **Command**: `bowtie2-build <assembly.fna> <assembly.fna>`.

3. **How MGEfinder uses these indexes**  
   After **find** and **pair**, MGEfinder runs **inferseq** steps (e.g. `inferseq_reference`, `inferseq_assembly`) to recover insertion sequences. Those steps use the Bowtie2 indexes to align and extend sequences against the reference and the sample assembly. So:

   - **BAM** (from BWA + formatbam) → used by **find** / **pair** to detect MGE signals.  
   - **Bowtie2 indexes** (from `index_genome_bowtie2` and `index_assembly`) → used by **inferseq** to resolve and extend the inserted sequences.

So the flow is: **BAM** (short reads vs reference) → **find** / **pair** → then **Bowtie2 indexes** (reference + assembly) are used by **inferseq** to build the final insertion sequences and database.

## Where Bowtie2 should run (environment)

Bowtie2 is used only in rules that already use the **mgefinder** conda env (the Snakefile has a global `conda: "../envs/mgefinder.yaml"`). So **Bowtie2 should be installed in the same environment as MGEfinder and BWA**: **mgefinder_env**. That way:

- **Snakemake** runs in the **snakemake** env (driver and DAG).
- **Rules** (bwa, bwa_index, formatbam, index_genome_bowtie2, index_assembly, find, pair, inferseq, …) run inside **mgefinder_env** when Snakemake is called with **--use-conda**.

Ensure **mgefinder_env** contains at least: **bwa**, **bowtie2**, and **mgefinder**. The repo’s `envs/mgefinder.yaml` lists bwa and bowtie2; add MGEfinder to that env (e.g. `pip install mgefinder`). The pipeline is run with **--use-conda** so Snakemake activates that env for the rules and `bowtie2-build` is found. If you already have `mgefinder_env`, you can add the tools with:  
`micromamba install -n mgefinder_env -c bioconda -c conda-forge bowtie2 bwa`
