#!/bin/bash
#SBATCH --job-name=Ecoli                         
#SBATCH --partition=cpu-short                                     
#SBATCH --time=8:00:00                               
#SBATCH --mem=2G
#SBATCH -n 10 # threaded 
#SBATCH -o slurm.Ecoli.out
#SBATCH -e slurm.Ecoli.err

# activate conda environment
source $HOME/.bash_profile
module load python
conda activate Ecoli_pigs

# change directory to where the bulk RNAseq Snakefile is located
CWD="scripts/"
cd $CWD

# run snakemake
snakemake -s 01_processing/Snakefile -j 52 --rerun-incomplete --latency-wait 10 --cluster "sbatch --ntasks 12 --partition=cpu-short --mem=32G -t 02:00:00"
