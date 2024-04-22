#!/bin/bash
#SBATCH --job-name=pip                         
#SBATCH --partition=cpu-med                                     
#SBATCH --time=08:00:00 # 8 hours                                
#SBATCH --mem=10G
#SBATCH -o slurm.pip.out
#SBATCH -e slurm.pip.err
#SBATCH --mail-user=olney.kimberly@mayo.edu

# source your bach profile to get your specific settings  
source $HOME/.bash_profile

module load python
conda activate Ecoli_pigs

# 1) get read information
#sh 01_sample_read_info.sh

# 2) create config
#python 02_create_pip_config.py

# 3) run snakemake - metaphlan alignment 
snakemake -s pip.Snakefile -j 3 --nolock --latency-wait 15 --rerun-incomplete --cluster "sbatch --ntasks 30 --partition=cpu-med --nodes 1 --mem=50G -t 36:00:00"


