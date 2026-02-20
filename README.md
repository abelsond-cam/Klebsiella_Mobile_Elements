# Klebsiella_Mobile_Elements

Mobile Genetic Elements in Klebsiella, identification and isolation

## Project Organization

```
├── CITATION.cff       <- Contains metadata on how the project might eventually be published.
├── LICENSE
├── Makefile           <- Makefile with commands like `make data` or `make train`
├── README.md          <- The top-level README for developers using this project.
├── config             <- Configuration for the pipeline.
│   ├── config.yaml    <- Single pipeline config (paths, MGEfinder parameters). See config/README.md.
│   └── .mge_merged_config.yaml  <- Generated at run time (gitignored).
├── envs               <- Conda env specs for Snakemake rules (e.g. envs/mgefinder.yaml). See envs/README.md.
├── docs               <- A default Sphinx project; see sphinx-doc.org for details
├── img                <- A place to store images associated with the project/pipeline, e.g. a
│                         figure of the pipeline DAG.
├── notebooks          <- Jupyter or Rmd notebooks. Naming convention is a number (for ordering),
│                         the creator's initials, and a short `-` delimited description, e.g.
│                         `1.0-jqp-initial-data-exploration`.
├── references         <- Data dictionaries, manuals, and all other explanatory materials.
├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
│   └── figures        <- Generated graphics and figures to be used in reporting
├── resources          <- Place for data. By default excluded from the git repository.
│   ├── external       <- Data from third party sources.
│   └── raw_data       <- The original, immutable data dump.
├── results            <- Final output of the data processing pipeline. By default excluded from the git repository.
├── sandbox            <- A place to test scripts and ideas. By default excluded from the git repository.
├── scripts            <- A place for short shell or python scripts.
├── setup.py           <- Makes project pip installable (pip install -e .) so src can be imported
├── src                <- Source code for use in this project.
│   └── __init__.py    <- Makes src a Python module
├── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io
├── workflow           <- Main pipeline (Snakefile includes mgefinder.end2end.snakefile).
│   ├── rules          <- Optional .smk modules.
│   └── Snakefile      <- Entrypoint that includes mgefinder.end2end.snakefile.
└── workspace          <- Space for intermediate results in the pipeline. By default excluded from the git repository.
```

---

*Project inspired by the [cookiecutter data science project template](https://drivendata.github.io/cookiecutter-data-science/) and the [Snakemake workflow template](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html).*
