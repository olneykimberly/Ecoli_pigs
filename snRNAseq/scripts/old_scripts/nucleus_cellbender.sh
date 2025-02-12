#!/bin/bash                         
#SBATCH --partition=cpu-med                                     
#SBATCH --time=32:00:00 # 8 hours                                
#SBATCH --mem=50G
#SBATCH --mail-user=olney.kimberly@mayo.edu

# activate conda env
source $HOME/.bash_profile
module load python
conda activate cellbender_v0.3.0

# set variable
sample=$1

# print what sample you're on
echo "sample: $sample"

# change directory to your desired output folder
cd /research/labs/neurology/fryer/m239830/Ecoli_pigs/snRNAseq/cellbender/$sample

# run cellbender with cellranger count raw_feature_bc_matrix.h5 output from each sample
cellbender remove-background \
		   --input /research/labs/neurology/fryer/m239830/Ecoli_pigs/snRNAseq/pipseeker_output/$sample/raw_feature_bc_matrix.h5 \
		   --output /research/labs/neurology/fryer/m239830/Ecoli_pigs/snRNAseq/cellbender/$sample/$sample_cellbender.h5 \

touch $sample

