# Envs

- **mgefinder.yaml** – Names the env used when **running Snakemake** (`mgefinder_env`). The pipeline does **not** use `--use-conda`: no conda is invoked by Snakemake. You manage `mgefinder_env` yourself (e.g. with micromamba). It must contain: **snakemake**, **bwa**, **bowtie2**, **mgefinder**. Data prep (validate, dataset generation) runs in the **snakemake** env; the Snakemake process and all rule jobs run in **mgefinder_env**.
