.PHONY: clean data data-legacy lint requirements fetch_data pipeline pipeline-verbose test-n1 dry-run test-dry test-dry-skip-dl submit-hpc recreate-symlinks test-mgefinder-env test_environment create_environment help

#################################################################################
# GLOBALS                                                                       #
#################################################################################

PROJECT_DIR := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
PROJECT_NAME = Klebsiella_Mobile_Elements
MGE_ENV = mgefinder_env
SNAKE_ENV = snakemake
FASTQ_ENV = fastq-dl

# Activate mgefinder_env (micromamba) - for MGEfinder pipeline only
define activate-mge
	eval "$$(micromamba shell hook --shell bash)" && micromamba activate $(MGE_ENV)
endef

# Activate snakemake env (micromamba) - for data prep and validation
define activate-snake
	eval "$$(micromamba shell hook --shell bash)" && micromamba activate $(SNAKE_ENV)
endef

# Activate fastq-dl env (micromamba) - for FASTQ downloads
define activate-fastq
	eval "$$(micromamba shell hook --shell bash)" && micromamba activate $(FASTQ_ENV)
endef

# use this function to load modules on the cluster
# usage in Makefile: $(call module-load,snakemake) && snakemake -j1 ...
define module-load
        eval `/usr/bin/modulecmd bash load $(1)`
endef


#################################################################################
# COMMANDS                                                                      #
#################################################################################

## Install Python Dependencies
requirements: test_environment
	#python3 -m pip install -U pip setuptools wheel
	#python3 -m pip install -r requirements.txt

## Make Dataset (multi-environment: snakemake for pipeline, mgefinder_env via conda)
data:
	$(call activate-snake) && python3 src/run_pipeline.py --config config/config.yaml

## Test with 1 sample only (safe for login node; skips FASTQ download - use if sample already present)
test-n1:
	$(call activate-snake) && python3 src/run_pipeline.py --config config/config.yaml --test-n 1 --verbose --skip-download

## Dry run: show what would be done (safe for login node)
dry-run:
	$(call activate-snake) && python3 src/run_pipeline.py --config config/config.yaml --dry-run --verbose

## Test dry run with 1 sample (safest for login node)
test-dry:
	$(call activate-snake) && python3 src/run_pipeline.py --config config/config.yaml --test-n 1 --dry-run --verbose

## Skip downloads if FASTQ already exists (useful for testing)
test-dry-skip-dl:
	$(call activate-snake) && python3 src/run_pipeline.py --config config/config.yaml --test-n 1 --dry-run --verbose --skip-download

## Submit full pipeline to HPC scheduler (adjust scripts/submit_hpc.sh first)
submit-hpc:
	sbatch scripts/submit_hpc.sh

## Recreate symlinks with correct paths (removes existing ones first)
recreate-symlinks:
	@echo ">>> Removing existing symlinks and recreating with correct paths"
	@ASSEMBLY_DIR=$$(python3 -c "import yaml; c=yaml.safe_load(open('config/config.yaml')); print(c.get('data_dir', c['wd']) + '/' + c['assemblies_dir'])") && \
	echo "Assembly dir: $$ASSEMBLY_DIR" && \
	cd "$$ASSEMBLY_DIR" && rm -f *.fa.gz *.fna.gz && \
	python3 $(PROJECT_DIR)/scripts/bulk_symlink_assemblies.py

## Test if mgefinder environment can run snakemake
test-mgefinder-env:
	@echo ">>> Testing if mgefinder_env can run snakemake..."
	$(call activate-mge) && snakemake --version

## Make Dataset (legacy bash version)
data-legacy:
	fetch_data
	pipeline

## Delete all compiled Python files
clean:
	find . -type f -name "*.py[co]" -delete
	find . -type d -name "__pycache__" -delete

## Lint using flake8 (requires mgefinder_env)
lint:
	$(call activate-mge) && flake8 src

## Fetch the raw data: from reference_comparison_sets.tsv get comparison IDs, check FASTQ, download if missing, generate mgefinder_dataset.txt
fetch_data:
	@$(call activate-snake) && \
	WD=$$(python3 -c "import yaml; print(yaml.safe_load(open('config/config.yaml'))['wd'])") && \
	DATA=$$(python3 -c "import yaml; c=yaml.safe_load(open('config/config.yaml')); print(c.get('data_dir', c['wd']))") && \
	echo ">>> Fetch data (wd=$$WD, data_dir=$$DATA)" && \
	ids=$$(python3 -c "import pandas as pd; df=pd.read_csv('$$DATA/processed/mgefinder/reference_comparison_sets.tsv', sep='\t'); row=df.iloc[0]; print(','.join([x.strip() for x in str(row['mge_comparison_set']).split(',')]))" 2>/dev/null) && \
	if [ -n "$$ids" ]; then echo "Checking FASTQ for comparison set..."; python3 src/run_fastq_download.py --config config/config.yaml --ids $$(echo "$$ids" | tr ',' ' ') 2>/dev/null || ./scripts/run_fastq_download.sh --config config/config.yaml --ids $$(echo "$$ids" | tr ',' ' '); fi && \
	python3 src/generate_mgefinder_dataset.py --config config/config.yaml --row 0 --out "$$WD/mgefinder_dataset.txt" && \
	echo ">>> Dataset written to $$WD/mgefinder_dataset.txt"

## Set up python interpreter environment (use micromamba; extend mgefinder_env manually)
create_environment:
	@echo ">>> Use: micromamba create -n $(MGE_ENV) ... and extend manually as needed"
	@echo ">>> Activate with: micromamba activate $(MGE_ENV)"

## Test mgefinder_env (the only environment we use after upgrade)
test_environment:
	@echo ">>> Testing mgefinder_env (our single environment for everything)..."
	$(call activate-mge) && python3 scripts/diagnose_env.py

#################################################################################
# PROJECT RULES                                                                 #
#################################################################################

## Execute the data pipeline: generate dataset for first row, then run Snakemake with config and WD as directory
pipeline:
	@$(call activate-snake) && \
	WD=$$(python3 -c "import yaml; print(yaml.safe_load(open('config/config.yaml'))['wd'])") && \
	DATA=$$(python3 -c "import yaml; c=yaml.safe_load(open('config/config.yaml')); print(c.get('data_dir', c['wd']))") && \
	REF=$$(python3 -c "import pandas as pd; df=pd.read_csv('$$DATA/processed/mgefinder/reference_comparison_sets.tsv', sep='\t'); print(df.iloc[0]['reference_sample_name'])") && \
	echo ">>> Running MGEfinder pipeline (wd=$$WD, reference=$$REF)" && \
	python3 src/validate_reference.py --config config/config.yaml --row 0 && \
	python3 src/generate_mgefinder_dataset.py --config config/config.yaml --row 0 --out "$$WD/mgefinder_dataset.txt" && \
	echo "genomes: [$$REF]" > $(PROJECT_DIR)/.mge_genomes.yaml && \
	snakemake --configfile config/config.yaml --configfile $(PROJECT_DIR)/.mge_genomes.yaml --directory "$$WD" -j 1 all

## Same as pipeline but with --printshellcmds and -p so you can verify each step as it runs
pipeline-verbose:
	@$(call activate-snake) && \
	WD=$$(python3 -c "import yaml; print(yaml.safe_load(open('config/config.yaml'))['wd'])") && \
	DATA=$$(python3 -c "import yaml; c=yaml.safe_load(open('config/config.yaml')); print(c.get('data_dir', c['wd']))") && \
	REF=$$(python3 -c "import pandas as pd; df=pd.read_csv('$$DATA/processed/mgefinder/reference_comparison_sets.tsv', sep='\t'); print(df.iloc[0]['reference_sample_name'])") && \
	echo ">>> Running MGEfinder pipeline VERBOSE (wd=$$WD, reference=$$REF)" && \
	python3 src/validate_reference.py --config config/config.yaml --row 0 && \
	python3 src/generate_mgefinder_dataset.py --config config/config.yaml --row 0 --out "$$WD/mgefinder_dataset.txt" && \
	echo "genomes: [$$REF]" > $(PROJECT_DIR)/.mge_genomes.yaml && \
	snakemake --configfile config/config.yaml --configfile $(PROJECT_DIR)/.mge_genomes.yaml --directory "$$WD" --printshellcmds -p -j 1 all

#################################################################################
# Self Documenting Commands                                                     #
#################################################################################

.DEFAULT_GOAL := help

# Inspired by <http://marmelab.com/blog/2016/02/29/auto-documented-makefile.html>
# sed script explained:
# /^##/:
# 	* save line in hold space
# 	* purge line
# 	* Loop:
# 		* append newline + line to hold space
# 		* go to next line
# 		* if line starts with doc comment, strip comment character off and loop
# 	* remove target prerequisites
# 	* append hold space (+ newline) to line
# 	* replace newline plus comments by `---`
# 	* print line
# Separate expressions are necessary because labels cannot be delimited by
# semicolon; see <http://stackoverflow.com/a/11799865/1968>
.PHONY: help
help:
	@echo "$$(tput bold)Available rules:$$(tput sgr0)"
	@echo
	@sed -n -e "/^## / { \
		h; \
		s/.*//; \
		:doc" \
		-e "H; \
		n; \
		s/^## //; \
		t doc" \
		-e "s/:.*//; \
		G; \
		s/\\n## /---/; \
		s/\\n/ /g; \
		p; \
	}" ${MAKEFILE_LIST} \
	| LC_ALL='C' sort --ignore-case \
	| awk -F '---' \
		-v ncol=$$(tput cols) \
		-v indent=19 \
		-v col_on="$$(tput setaf 6)" \
		-v col_off="$$(tput sgr0)" \
	'{ \
		printf "%s%*s%s ", col_on, -indent, $$1, col_off; \
		n = split($$2, words, " "); \
		line_length = ncol - indent; \
		for (i = 1; i <= n; i++) { \
			line_length -= length(words[i]) + 1; \
			if (line_length <= 0) { \
				line_length = ncol - indent - length(words[i]) - 1; \
				printf "\n%*s ", -indent, " "; \
			} \
			printf "%s ", words[i]; \
		} \
		printf "\n"; \
	}' \
	| more $(shell test $(shell uname) = Darwin && echo '--no-init --raw-control-chars')
