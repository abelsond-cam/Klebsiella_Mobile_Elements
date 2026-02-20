# The Snakefile for MGEfinder which runs bwa, and indexes the bam file and  the assembly before running the standard MGEfinder pipeline. It also removes the bam file after finishing. This is the command I used to run it at the Sanger.
# snakemake -s mgefinder.end2end.snakefile --configfile /nfs/users/nfs_a/aw27/aw27/conda/mgefinder/lib/python3.8/site-packages/mgefinder-1.0.6-py3.8.egg/mgefinder/workflow/denovo.original.config.yml   --config memory=16000 wd=  --profile /nfs/users/nfs_a/aw27/.config/snakemake/lsf --default-resources mem_mb=16000 --restart-times 3 --rerun-incomplete
# It currently expects a table mgefinder_dataset.txt with the file locations which is parsed by Snakemake (at the begining of the pipeline) into dictionaries. If all the assemblies/reads are in the the same directory, this can be simplified. You'll have to re-download the reads yourself. Probably easier to explain on a quick call.

# MGEfinder end-to-end: BWA, formatbam, Bowtie2 index, MGEfinder find/pair/inferseq/makedatabase/clusterseq/genotype/summarize/makefasta.
# Entrypoint: workflow/Snakefile (includes this file). Run: snakemake -s workflow/Snakefile --configfile config/config.yaml --directory <WD> [targets]
# Dataset TSV: data_dir, sample_id, sample_name, gff, contigs [, fastq1, fastq2 ]. If fastq1/fastq2 columns present, bwa uses them; else assumes data_dir/sample_id_1.fastq.gz. SEGMENT comments mark parts for isolated testing.
# No conda: directive; run Snakemake with your env activated (e.g. micromamba activate mgefinder_env).

import os
from os.path import basename, join

WD = config.get("wd", "")
if WD:
    WD = os.path.abspath(WD)
DATA_DIR = config.get("data_dir", WD)
if DATA_DIR:
    DATA_DIR = os.path.abspath(DATA_DIR)

GENOME_DIR = join(WD, config.get("genome_dir", "00.genome"))
ASSEMBLY_DIR = join(WD, config.get('assembly_dir', "00.assembly"))
BAM_DIR = join(WD, config.get('bam_dir', "00.bam"))
BWA_DIR = join(WD, "bwa")
MUSTACHE_DIR = join(WD, config.get('mgefinder_dir', "mgefinder"))
DATABASE_DIR = join(WD, config.get('database_dir', "database"))
RESULTS_DIR = join(WD, config.get('results_dir', "results"))

# Debug: print where we read (DATA_DIR) vs where we write (WD and its subdirs)
print("MGEfinder paths: data_dir (read) =", DATA_DIR, "| wd (write) =", WD)
print("  subdirs under wd: genome_dir =", GENOME_DIR, "| assembly_dir =", ASSEMBLY_DIR, "| bam_dir =", BAM_DIR, "| bwa_dir =", BWA_DIR)

WC_genomes = glob_wildcards(join(GENOME_DIR, '{genome}.fna'))
WC_bam_samples = glob_wildcards(join(BAM_DIR, "{sample}.{genome}.bam"))

meta = join(WD, config.get("mgefinder_dataset", "mgefinder_dataset.txt"))

def _parse_meta_line(l):
    parts = l.strip().split('\t')
    return parts[0], parts[1], parts[2], parts[3], parts[4]

def get_assembly_dict():
    sample2filename = {}
    with open(meta) as f:
        f.readline() 
        for l in f:
            data_dir, sample_id, sample_name, gff, contigs = _parse_meta_line(l)
            sample2filename[sample_name] = contigs
    return sample2filename 

def get_data_dict():
    sample2filename = {}
    with open(meta) as f:
        f.readline() 
        for l in f:
            data_dir, sample_id, sample_name, gff, contigs = _parse_meta_line(l)
            sample2filename[sample_name] = data_dir 
    return sample2filename 

def get_sample_dict():
    sample2filename = {}
    with open(meta) as f:
        f.readline() 
        for l in f:
            data_dir, sample_id, sample_name, gff, contigs = _parse_meta_line(l)
            sample2filename[sample_name] = sample_id 
    return sample2filename 

def get_fastq1_dict():
    """sample_name -> path to first FASTQ. Uses fastq1 column if present (7-col), else data_dir/sample_id_1.fastq.gz (5-col)."""
    out = {}
    with open(meta) as f:
        f.readline()
        for l in f:
            parts = l.strip().split('\t')
            data_dir, sample_id, sample_name, gff, contigs = parts[0], parts[1], parts[2], parts[3], parts[4]
            if len(parts) >= 7:
                out[sample_name] = parts[5]
            else:
                out[sample_name] = data_dir + "/" + sample_id + "_1.fastq.gz"
    return out

def get_fastq2_dict():
    """sample_name -> path to second FASTQ. Uses fastq2 column if present (7-col), else data_dir/sample_id_2.fastq.gz (5-col)."""
    out = {}
    with open(meta) as f:
        f.readline()
        for l in f:
            parts = l.strip().split('\t')
            data_dir, sample_id, sample_name, gff, contigs = parts[0], parts[1], parts[2], parts[3], parts[4]
            if len(parts) >= 7:
                out[sample_name] = parts[6]
            else:
                out[sample_name] = data_dir + "/" + sample_id + "_2.fastq.gz"
    return out

def get_assembly():
    sample2filename = {}
    with open(meta) as f:
        f.readline() 
        for l in f:
            data_dir, sample_id, sample_name, gff, contigs = _parse_meta_line(l)
            sample2filename[sample_name] = contigs
    return sample2filename 

def get_reference_path(wildcards):
    """Resolve reference assembly path (supports .fna, .fa, .fna.gz, .fa.gz). Assemblies live under data_dir."""
    base = join(DATA_DIR, config["assemblies_dir"], wildcards.genome)
    for ext in [".fna", ".fa", ".fna.gz", ".fa.gz"]:
        p = base + ext
        if os.path.exists(p):
            return p
    return base + ".fna"

assembly_dict = get_assembly_dict()
data_dict = get_data_dict()
sample_dict = get_sample_dict()
fastq1_dict = get_fastq1_dict()
fastq2_dict = get_fastq2_dict()

def get_samples():
    samples = [] 
    with open(meta) as f:
        f.readline() 
        i = 0
        for l in f:
            #if i == 1:
            #    break
            i+=1
            data_dir, sample_id, sample_name, gff, contigs = _parse_meta_line(l)
            samples.append(sample_name)
    return samples 


SAMPLES = get_samples()
# Get genomes from config (set by run_pipeline.py in merged config)
try:
    GENOMES = config["genomes"]
    if not GENOMES:
        raise KeyError("genomes list is empty")
except KeyError:
    raise SystemExit("ERROR: 'genomes' not found in config. Run the pipeline via: make test-dry-skip-dl (or data/dry-run) so run_pipeline.py injects genomes.")


rule all:
    input:
        expand(join(RESULTS_DIR, "{genome}/04.makefasta.{genome}.all_seqs.fna"), genome=GENOMES)
    run:
        pass

# --- SEGMENT: Setup (copy_genome, copy_assembly) ---
# Reference from config (assemblies_dir + genome name); assemblies from dataset TSV (contigs column).
rule copy_genome:
    input:
        ref = lambda wc: get_reference_path(wc),
    output:
        join(GENOME_DIR, "{genome}.fna"),
    shell:
        """
        echo ">>> copy_genome: read {input.ref} write {output}"
        if [[ "{input.ref}" == *.gz ]]; then
            gunzip -c "{input.ref}" > "{output}"
        else
            cp "{input.ref}" "{output}"
        fi
        """

rule copy_assembly:
    input:
        assembly = lambda wc: assembly_dict[wc.sample],
    output:
        join(ASSEMBLY_DIR, "{sample}.fna"),
    shell:
        """
        echo ">>> copy_assembly: read {input.assembly} write {output}"
        if [[ "{input.assembly}" == *.gz ]]; then
            gunzip -c "{input.assembly}" > "{output}"
        else
            ln -sf "{input.assembly}" "{output}"
        fi
        """

# --- SEGMENT: BWA (bwa_index, bwa, formatbam) ---
rule formatbam:
    input:
        join(BWA_DIR, "{sample}.{genome}.bwa.sam")
    output:
        temp(join(BAM_DIR, "{sample}.{genome}.bam"))
    shell:
        """
        echo ">>> formatbam: read {input} write {output}"
        mgefinder formatbam  {input} {output}
        """


rule bwa_index:
    input:
        join(GENOME_DIR, "{genome}.fna")
    output:
        join(GENOME_DIR, "{genome}.fna.amb"),
    shell:
        """
        echo ">>> bwa_index: read {input} write {output}"
        bwa index {input} {output}
        """

rule bwa:
    input:
        index = join(GENOME_DIR, "{genome}.fna.amb"),
        genome = join(GENOME_DIR, "{genome}.fna"),
        fastq1 = lambda wc: fastq1_dict[wc.sample],
        fastq2 = lambda wc: fastq2_dict[wc.sample]
    params:
        #data_dir = lambda wc: data_dict[wc.sample],       
        #sample_id = lambda wc: sample_dict[wc.sample],
    output:
        sample=temp(join(BWA_DIR, "{sample}.{genome}.bwa.sam"))
    shell:
        """
        echo ">>> bwa: read {input.genome} {input.fastq1} {input.fastq2} write {output.sample}"
        bwa mem {input.genome} {input.fastq1} {input.fastq2}   -o  {output.sample}
        """

# --- SEGMENT: Index (index_genome_bowtie2, index_assembly) ---
rule index_genome_bowtie2:
    log: join(GENOME_DIR, "log/{genome}.index_bowtie2.log")
    benchmark: join(GENOME_DIR, "log/{genome}.index_bowtie2.benchmark.txt")
    input:
        join(GENOME_DIR, "{genome}.fna")
    output:
        one=join(GENOME_DIR, "{genome}.fna.1.bt2"),
        two=join(GENOME_DIR, "{genome}.fna.2.bt2"),
        three=join(GENOME_DIR, "{genome}.fna.3.bt2"),
        four=join(GENOME_DIR, "{genome}.fna.4.bt2"),
        revone=join(GENOME_DIR, "{genome}.fna.rev.1.bt2"),
        revtwo=join(GENOME_DIR, "{genome}.fna.rev.2.bt2")
    shell:
        """
        echo ">>> index_genome_bowtie2: indexing {wildcards.genome}"
        bowtie2-build {input} {input} 1> {log} 2> {log}.err || \
        (cat {log}.err; exit 1)
        """

rule index_assembly:
    log: join(ASSEMBLY_DIR, "log/{sample}.index_assembly.log")
    benchmark: join(ASSEMBLY_DIR, "log/{sample}.index_assembly.benchmark.txt")
    input:
        contigs=join(ASSEMBLY_DIR, "{sample}.fna")
    output:
        one=join(ASSEMBLY_DIR, "{sample}.fna.1.bt2"),
        two=join(ASSEMBLY_DIR, "{sample}.fna.2.bt2"),
        three=join(ASSEMBLY_DIR, "{sample}.fna.3.bt2"),
        four=join(ASSEMBLY_DIR, "{sample}.fna.4.bt2"),
        revone=join(ASSEMBLY_DIR, "{sample}.fna.rev.1.bt2"),
        revtwo=join(ASSEMBLY_DIR, "{sample}.fna.rev.2.bt2")
    shell:
        """
        echo ">>> index_assembly: indexing assembly {wildcards.sample}"
        bowtie2-build {input} {input} 1> {log} 2> {log}.err || \
        (cat {log}.err; exit 1)
        """

# --- SEGMENT: MGEfinder find / pair ---
rule find:
    log: join(MUSTACHE_DIR, "{genome}/{sample}/log/{sample}.{genome}.find.log")
    benchmark: join(MUSTACHE_DIR, "{genome}/{sample}/log/{sample}.{genome}.find.benchmark.txt")
    input:
        bam=join(BAM_DIR, "{sample}.{genome}.bam"),
    output:
        find=join(MUSTACHE_DIR, '{genome}/{sample}/01.find.{sample}.{genome}.tsv')
    params:
        sample='{sample}',
        minlen=config['find']['minlen'],
        mincount=config['find']['mincount'],
        minq=config['find']['minq'],
        minial=config['find']['minial'],
        mindist=config['find']['mindist'],
        minratio=config['find']['minratio'],
        maxir=config['find']['maxir'],
        lins=config['find']['lins'],
        mcc=config['find']['mcc'],
        check_bwa=config['find']['check_bwa_flag']
    shell:
        """
        echo ">>> find: read {input.bam} write {output.find}"
        mgefinder find -id {params.sample} -minlen {params.minlen} -mincount {params.mincount} -minq {params.minq} \
        -minial {params.minial} -mindist {params.mindist} -minratio {params.minratio} -maxir {params.maxir} -lins \
        {params.lins} -mcc {params.mcc} {params.check_bwa} {input.bam} -o {output.find} 1> {log} 2> {log}.err || \
        (cat {log}.err; exit 1)
        """


rule pair:
    log: join(MUSTACHE_DIR, "{genome}/{sample}/log/{sample}.{genome}.pair.log")
    benchmark: join(MUSTACHE_DIR, "{genome}/{sample}/log/{sample}.{genome}.pair.benchmark.txt")
    input:
        find=ancient(join(MUSTACHE_DIR, '{genome}/{sample}/01.find.{sample}.{genome}.tsv')),
        bam=join(BAM_DIR, "{sample}.{genome}.bam"),
        genome=join(GENOME_DIR, "{genome}.fna")
    output:
        pair=join(MUSTACHE_DIR, '{genome}/{sample}/02.pair.{sample}.{genome}.tsv')
    params:
        maxdr=config['pair']['maxdr'],
        minq=config['pair']['minq'],
        minial=config['pair']['minial'],
        maxjsp=config['pair']['maxjsp'],
        lins=config['pair']['lins']
    shell:
        """
        echo ">>> pair: MGEfinder pair {wildcards.sample} {wildcards.genome}"
        mgefinder pair -maxdr {params.maxdr} -minq {params.minq} -minial {params.minial} -maxjsp {params.maxjsp} \
        -lins {params.lins} {input.find} {input.bam} {input.genome} -o {output.pair} 1> {log} 2> {log}.err || \
        (cat {log}.err; exit 1)
        """

# --- SEGMENT: MGEfinder inferseq (assembly, reference, overlap) ---
rule inferseq_assembly:
    log: join(MUSTACHE_DIR, "{genome}/{sample}/log/{sample}.{genome}.inferseq_assembly.log")
    benchmark: join(MUSTACHE_DIR, "{genome}/{sample}/log/{sample}.{genome}.inferseq_assembly.benchmark.txt")
    input:
        pair=join(MUSTACHE_DIR, '{genome}/{sample}/02.pair.{sample}.{genome}.tsv'),
        bam=join(BAM_DIR, "{sample}.{genome}.bam"),
        recover_reference=join(GENOME_DIR, "{genome}.fna"),
        recover_assembly=join(ASSEMBLY_DIR, "{sample}.fna"),
        one_ref=join(GENOME_DIR, '{genome}.fna.1.bt2'),
        one_asm=join(ASSEMBLY_DIR, '{sample}.fna.1.bt2')
    output:
        recover=join(MUSTACHE_DIR, '{genome}/{sample}/03.inferseq_assembly.{sample}.{genome}.tsv')
    params:
        minident=config['inferseq_assembly']['minident'],
        maxclip=config['inferseq_assembly']['maxclip'],
        maxsize=config['inferseq_assembly']['maxsize'],
        minsize=config['inferseq_assembly']['minsize']
    shell:
        """
        echo ">>> inferseq_assembly: {wildcards.sample} {wildcards.genome}"
        mgefinder inferseq-assembly -minident {params.minident} -maxclip {params.maxclip} -maxsize {params.maxsize} \
        -minsize {params.minsize} {input.pair} {input.bam} {input.recover_assembly} {input.recover_reference} -o {output.recover} 1> {log} 2> {log}.err || \
        (cat {log}.err; exit 1)
        """


rule inferseq_reference:
    log: join(MUSTACHE_DIR, "{genome}/{sample}/log/{sample}.{genome}.inferseq_reference.log")
    benchmark: join(MUSTACHE_DIR, "{genome}/{sample}/log/{sample}.{genome}.inferseq_reference.benchmark.txt")
    input:
        pair=join(MUSTACHE_DIR, '{genome}/{sample}/02.pair.{sample}.{genome}.tsv'),
        recover_reference=join(GENOME_DIR, "{genome}.fna"),
        one_ref=join(GENOME_DIR, '{genome}.fna.1.bt2')
    output:
        recover=join(MUSTACHE_DIR, '{genome}/{sample}/03.inferseq_reference.{sample}.{genome}.tsv')
    params:
        minident=config['inferseq_reference']['minident'],
        maxclip=config['inferseq_reference']['maxclip'],
        maxsize=config['inferseq_reference']['maxsize'],
        minsize=config['inferseq_reference']['minsize']
    shell:
        """
        echo ">>> inferseq_reference: {wildcards.sample} {wildcards.genome}"
        mgefinder inferseq-reference -minident {params.minident} -maxclip {params.maxclip} -maxsize {params.maxsize} \
        -minsize {params.minsize} {input.pair} {input.recover_reference} -o {output.recover} 1> {log} 2> {log}.err || \
        (cat {log}.err; exit 1)
        """


rule inferseq_overlap:
    log: join(MUSTACHE_DIR, "{genome}/{sample}/log/{sample}.{genome}.inferseq_overlap.log")
    benchmark: join(MUSTACHE_DIR, "{genome}/{sample}/log/{sample}.{genome}.benchmark.txt")
    input:
        pair=join(MUSTACHE_DIR, '{genome}/{sample}/02.pair.{sample}.{genome}.tsv'),
    output:
        outfile=join(MUSTACHE_DIR, '{genome}/{sample}/03.inferseq_overlap.{sample}.{genome}.tsv')
    params:
        minscore=config['inferseq_overlap']['minscore'],
        minopi=config['inferseq_overlap']['minopi'],
        minsize=config['inferseq_overlap']['minsize']
    shell:
        """
        echo ">>> inferseq_overlap: {wildcards.sample} {wildcards.genome}"
        mgefinder inferseq-overlap -minscore {params.minscore} -minopi {params.minopi} -minsize {params.minsize} \
        {input.pair} -o {output.outfile} 1> {log} 2> {log}.err
        """

# --- SEGMENT: MGEfinder database (make_inferseq_file_path_list, make_database, inferseq_database, file lists) ---
rule make_inferseq_file_path_list:
    input:
        expand(join(MUSTACHE_DIR, '{{genome}}/{sample}/03.inferseq_assembly.{sample}.{{genome}}.tsv'), sample=SAMPLES),
        expand(join(MUSTACHE_DIR, '{{genome}}/{sample}/03.inferseq_reference.{sample}.{{genome}}.tsv'), sample=SAMPLES),
        expand(join(MUSTACHE_DIR, '{{genome}}/{sample}/03.inferseq_overlap.{sample}.{{genome}}.tsv'), sample=SAMPLES)
    output:
        outfile=join(MUSTACHE_DIR, '{genome}/{genome}.all_inferseq.txt')
    run:
        with open(output.outfile, 'w') as outfile:
            for f in input:
                outfile.write(f+'\n')

rule make_database:
    benchmark: join(DATABASE_DIR, '{genome}/{genome}.database.benchmark.txt')
    input:
        join(MUSTACHE_DIR, '{genome}/{genome}.all_inferseq.txt')
    output:
        db=join(DATABASE_DIR, '{genome}/{genome}.database.fna'),
        index=join(DATABASE_DIR, '{genome}/{genome}.database.fna.1.bt2')
    threads:
        16
    params:
        memory=config['memory'],
        outdir=join(DATABASE_DIR, '{genome}'),
        prefix='{genome}.database',
        minsize=config['makedatabase']['minsize'],
        maxsize=config['makedatabase']['maxsize'],
    shell:
        """
        echo ">>> make_database: building database for {wildcards.genome}"
        mgefinder makedatabase -minsize {params.minsize} -maxsize {params.maxsize} --threads {threads} \
        --memory {params.memory} -o {params.outdir} -p {params.prefix} {input}  --force
        """


rule inferseq_database:
    log: join(MUSTACHE_DIR, "{genome}/{sample}/log/{sample}.{genome}.inferseq_database.log")
    benchmark: join(MUSTACHE_DIR, "{genome}/{sample}/log/{sample}.{genome}.inferseq_database.benchmark.txt")
    input:
        pair=join(MUSTACHE_DIR, '{genome}/{sample}/02.pair.{sample}.{genome}.tsv'),
        database=join(DATABASE_DIR, '{genome}/{genome}.database.fna'),
        one=join(DATABASE_DIR, '{genome}/{genome}.database.fna.1.bt2')
    output:
        outfile=join(MUSTACHE_DIR, '{genome}/{sample}/04.inferseq_database.{sample}.{genome}.tsv')
    params:
        minident=config['inferseq_database']['minident'],
        maxclip=config['inferseq_database']['maxclip'],
        maxedgedist=config['inferseq_database']['maxedgedist']
    shell:
        """
        echo ">>> inferseq_database: {wildcards.sample} {wildcards.genome}"
        mgefinder inferseq-database -minident {params.minident} -maxclip {params.maxclip} -maxedgedist \
        {params.maxedgedist} {input.pair} {input.database} -o {output.outfile} 1> {log} 2> {log}.err
        """

rule make_inferseq_database_file_path_list:
    input:
        expand(join(MUSTACHE_DIR, '{{genome}}/{sample}/03.inferseq_assembly.{sample}.{{genome}}.tsv'), sample=SAMPLES),
        expand(join(MUSTACHE_DIR, '{{genome}}/{sample}/03.inferseq_reference.{sample}.{{genome}}.tsv'), sample=SAMPLES),
        expand(join(MUSTACHE_DIR, '{{genome}}/{sample}/03.inferseq_overlap.{sample}.{{genome}}.tsv'), sample=SAMPLES),
        expand(join(MUSTACHE_DIR, '{{genome}}/{sample}/04.inferseq_database.{sample}.{{genome}}.tsv'), sample=SAMPLES)
    output:
        outfile=join(MUSTACHE_DIR, '{genome}/{genome}.all_inferseq_database.txt')
    run:
        with open(output.outfile, 'w') as outfile:
            for f in input:
                outfile.write(f+'\n')

# --- SEGMENT: MGEfinder clusterseq, genotype, summarize, makefasta ---
rule clusterseq:
    input:
        join(MUSTACHE_DIR, '{genome}/{genome}.all_inferseq_database.txt')
    output:
        clusterseq=join(RESULTS_DIR, '{genome}/01.clusterseq.{genome}.tsv')
    threads:
        16
    params:
        memory=config['memory'],
        prefix='{genome}.database',
        minsize=config['clusterseq']['minsize'],
        maxsize=config['clusterseq']['maxsize'],
    shell:
        """
        echo ">>> clusterseq: clustering sequences for {wildcards.genome}"
        mgefinder clusterseq -minsize {params.minsize} -maxsize {params.maxsize} --threads {threads} \
        --memory {params.memory} {input} -o {output}
        """

rule make_pair_file_path_list:
    input:
        expand(join(MUSTACHE_DIR, '{{genome}}/{sample}/02.pair.{sample}.{{genome}}.tsv'), sample=SAMPLES)
    output:
        outfile=join(MUSTACHE_DIR, '{genome}/{genome}.all_pair.txt')
    run:
        with open(output.outfile, 'w') as outfile:
            for f in input:
                outfile.write(f+'\n')


rule genotype:
    log: join(RESULTS_DIR, "{genome}/log/{genome}.genotype.log")
    benchmark: join(RESULTS_DIR, "{genome}/log/{genome}.genotype.benchmark.txt")
    input:
        join(RESULTS_DIR, '{genome}/01.clusterseq.{genome}.tsv'),
        outfile=join(MUSTACHE_DIR, '{genome}/{genome}.all_pair.txt')
    output:
        genotype=join(RESULTS_DIR, '{genome}/02.genotype.{genome}.tsv')
    params:
        filter_clusters=config['genotype']['filter_clusters']
    shell:
        """
        if [ "{params.filter_clusters}" == "True" ]; then
            mgefinder genotype --filter-clusters-inferred-assembly {input} -o {output} 1> {log} 2> {log}.err || \
            (cat {log}.err; exit 1)
        else
            mgefinder genotype --no-filter-clusters-inferred-assembly {input} -o {output} 1> {log} 2> {log}.err || \
            (cat {log}.err; exit 1)
        fi
        """

rule summarize:
    log: join(RESULTS_DIR, "{genome}/log/{genome}.summarize.log")
    benchmark: join(RESULTS_DIR, "{genome}/log/{genome}.summarize.benchmark.txt")
    input:
        join(RESULTS_DIR, '{genome}/01.clusterseq.{genome}.tsv'),
        join(RESULTS_DIR, '{genome}/02.genotype.{genome}.tsv')
    output:
        outfile1=join(RESULTS_DIR, '{genome}/03.summarize.{genome}.clusters.tsv'),
        outfile2=join(RESULTS_DIR, '{genome}/03.summarize.{genome}.groups.tsv')
    params:
        prefix=join(RESULTS_DIR, '{genome}/03.summarize.{genome}')
    shell:
        """
        echo ">>> summarize: summarizing results for {wildcards.genome}"
        mgefinder summarize {input} -o {params.prefix} 1> {log} 2> {log}.err || \
        (cat {log}.err; exit 1)
        """


rule makefasta:
    log: join(RESULTS_DIR, "{genome}/log/{genome}.makefasta.log")
    benchmark: join(RESULTS_DIR, "{genome}/log/{genome}.makefasta.benchmark.txt")
    input:
        join(RESULTS_DIR, '{genome}/01.clusterseq.{genome}.tsv'),
        join(RESULTS_DIR, '{genome}/03.summarize.{genome}.clusters.tsv')
    output:
        outfile1=join(RESULTS_DIR, '{genome}/04.makefasta.{genome}.all_seqs.fna'),
        outfile2=join(RESULTS_DIR, '{genome}/04.makefasta.{genome}.repr_seqs.fna')
    params:
        prefix=join(RESULTS_DIR, '{genome}/04.makefasta.{genome}')
    shell:
        """
        echo ">>> makefasta: extracting fasta for {wildcards.genome}"
        mgefinder makefasta {input} -o {params.prefix} 1> {log} 2> {log}.err || \
        (cat {log}.err; exit 1)
        """
