# COSMOS 2023 - Long-read analysis of C2C12 Differentiation.

# Instructions
1. Clone this repository in your folder on the UCI HPC
2. Run the `cosmos2023_c2c12_set-up.sh` script 
	- install package manager `conda` and `mamba`
	- configure `conda` by adding channels
	- clone long-read packages `TranscriptClean` and `TALON`
	- create environment with `snakemake` and `TALON` installed
3. Run `cosmos2023_c2c12_smk.sh`
	- Map to genome using `minimap2`
	- Reverse flipped reads with custom script
	- Reference-based error correction using `TranscriptClean`
	- Quantify and categorize transcripts using `TALON`
