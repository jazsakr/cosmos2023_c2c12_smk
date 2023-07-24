### ---------------------------------------------------
### make folder for slurm logs
### ---------------------------------------------------

mkdir slurm_logs

### ---------------------------------------------------
### install package managers, conda and mamba
### ---------------------------------------------------

#wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
#bash Miniconda3-latest-Linux-x86_64.sh

### answer "yes" when it asks "Do you wish the installer to initialize Miniconda3 by running conda init?"

#source ~/.bashrc

#conda config --add channels bioconda
#conda config --add channels conda-forge
#conda update conda

#conda install mamba -n base -c conda-forge
#mamba init

### ---------------------------------------------------
### install packages for long-read analysis
### ---------------------------------------------------

git clone https://github.com/mortazavilab/TranscriptClean.git ./packages/TranscriptClean
git clone https://github.com/mortazavilab/TALON.git ./packages/TALON

### ---------------------------------------------------
### create environment and install TALON in it 
### ---------------------------------------------------

mamba create -c conda-forge -c bioconda -n longread snakemake=7.12.0 tabulate=0.8 graphviz minimap2 samtools pysam pyfaidx pyranges
### answer "y" when it asks "Confirm changes: [Y/n]""

source ~/miniconda3/etc/profile.d/conda.sh
conda activate longread

pip3 install ./packages/TALON