#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate longread

python -m snakemake -s ./snakemake/cosmos2023_c2c12.smk --configfile ./snakemake/cosmos2023_c2c12_config.yaml -j 8 \
--cluster 'sbatch -A seyedam_lab --mem {resources.mem_gb}G --cpus-per-task {threads} --output=$PWD/slurm_logs/slurm-%j.out --time=06:00:00' \
--latency-wait 120 \
#--dag | dot -Tpng > ruledag.png