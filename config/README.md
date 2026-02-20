# Config

- **config.yaml** – Single source of truth for the MGEfinder pipeline (paths, dirs, MGEfinder parameters). Edit this for your paths and settings.
  - **data_dir**: where the pipeline *reads* from (assemblies, fastq, metadata, reference_comparison_sets).
  - **wd**: where the pipeline *writes* all working and output files; subdirs (genome_dir, assembly_dir, bam_dir, mgefinder, database, results) are created under `wd`.
- **.mge_merged_config.yaml** – Generated at run time by `src/run_pipeline.py` (config.yaml + `genomes` from reference_comparison_sets). Do not edit; it is gitignored.
