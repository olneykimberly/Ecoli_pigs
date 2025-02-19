#!/bin/bash
#SBATCH --job-name=Ecoli_cellranger                                                      
#SBATCH --time=30:00:00 # 8 hours                                
#SBATCH --mem=10G
#SBATCH -o slurm.Ecoli_cellranger.out
#SBATCH -e slurm.Ecoli_cellranger.err

# source your bach profile to get your specific settings  
source $HOME/.bash_profile

#conda activate cellbender
conda activate cellbender_py38
# STEPS
# 0) Rename fastq files to standard format. 
# sh 00_rename_fastq_files.sh

# 1) get read information
# sh 01_get_read_info.sh

# 2) create config
# python 02_create_config.py

# 3) run snakemake - metaphlan alignment 
#snakemake -s Snakefile -j 8 --nolock --latency-wait 15 --rerun-incomplete --cluster "sbatch --ntasks 16 --nodes 1 --mem=50G -t 36:00:00"
snakemake -s Snakefile.v3 -j 8 --nolock --latency-wait 15 --rerun-incomplete --cluster "sbatch --partition=gpu-a100 --ntasks 1 --nodes 1 --mem=50G --gres=gpu:1 -t 8:00:00"

