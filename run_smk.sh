#!/bin/bash

source ~/miniconda3/etc/profile.d/conda.sh
conda activate longread

python -m snakemake -s ./snakemake/cosmos2023_c2c12.smk --configfile ./snakemake/cosmos2023_c2c12_config.yaml -j 8 \
--cluster 'sbatch -A COSMOS2023_LAB --mem {resources.mem_gb}G --cpus-per-task {threads} --output=$PWD/slurm_logs/slurm-%j.out' \
--latency-wait 120 \
#--dag | dot -Tpng > c2c12_diff_ruledag.png