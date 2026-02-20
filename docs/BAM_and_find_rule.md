# What is a BAM file and what does it do in this pipeline?

## What is a BAM file?

**BAM** (Binary Alignment Map) is a compressed, indexed binary format for storing **alignments of sequencing reads to a reference genome**. It is the binary equivalent of SAM (Sequence Alignment Map).

- **Content**: For each short read (e.g. from your DRR FASTQs), it stores which position on the reference it aligns to, the alignment quality, and optional fields (flags, CIGAR, etc.).
- **Why it exists**: Alignments are needed to find “split” or “clipped” reads that span mobile genetic element (MGE) boundaries; MGEfinder uses these patterns to detect MGEs.

## What the BAM is doing in this pipeline

1. **Short reads (FASTQ)**  
   Paths to the paired FASTQs come from **mgefinder_dataset.txt**: the **fastq1** and **fastq2** columns list the actual files for each sample (e.g. DRR061421_1.fastq.gz and DRR061421_2.fastq.gz in that sample’s folder). The pipeline does not assume a fixed name like `sample_id_1.fastq.gz`.

2. **Alignment (rule `bwa`)**  
   `bwa mem` aligns those reads to the **reference genome** (e.g. GCF_000814305.1_ASM81430v1_genomic). It produces a **SAM** file (text alignments).

3. **Rule `formatbam`**  
   `mgefinder formatbam` converts that SAM into a **BAM** and ensures it’s in the form MGEfinder expects (sorted, indexed, etc.). The BAM is written under `00.bam/`.

4. **Rule `find` (and later steps)**  
   MGEfinder **find** (and then pair, inferseq, etc.) reads the BAM to detect MGE insertion signals from alignment patterns (clips, split reads). So the BAM is the **main input** to the MGEfinder “find” step.

So: **FASTQ → bwa (align to reference) → SAM → formatbam → BAM → MGEfinder find**.

The **sample name** in the BAM filename (e.g. SAMD00052619) is the **sample identifier** from your dataset (the row in mgefinder_dataset.txt). The **reads inside** that BAM come from the **fastq1** and **fastq2** paths listed for that row (e.g. DRR061421_1/2.fastq.gz). So the BAM name uses the sample id; the actual read filenames are whatever is in the dataset columns.

## Why you got “Missing input files for rule find”

Snakemake said the **input** to rule `find` is missing:  
`00.bam/SAMD00052619.GCF_000814305.1_ASM81430v1_genomic.bam`.

That file is supposed to be **produced by** rule `formatbam`. So either:

1. **Path mismatch**  
   The path used for the **output** of `formatbam` was not the same as the path used for the **input** of `find`. Then Snakemake does not see that `formatbam` creates the file `find` needs, so it treats that BAM as an external missing input and raises MissingInputException.

2. **Upstream not run**  
   If the path were correct, Snakemake would schedule `formatbam` (and before it `bwa`, etc.) to create the BAM. So the most likely cause here is (1): **output of formatbam** was given as a **relative** path (`00.bam/...`) while **input of find** was given as an **absolute** path (`join(BAM_DIR, ...)`). With `--directory` set, those can resolve to the same place on disk but still be different strings in the DAG, so Snakemake doesn’t link them.

The fix is to make **formatbam’s output** use the same path as **find’s input**: `join(BAM_DIR, "{sample}.{genome}.bam")` (with `temp()` if you want the file to be temporary). That way the DAG correctly sees that `formatbam` produces the file that `find` needs.
