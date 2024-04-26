library(DropletUtils)
library(Seurat)

raw_matrix <- Read10X("/research/labs/neurology/fryer/m239830/Ecoli_pigs/snRNAseq/pipseeker_output/E2A/raw_matrix/")
write10xCounts("/research/labs/neurology/fryer/m239830/Ecoli_pigs/snRNAseq/pipseeker_output/E2A/raw_feature_bc_matrix.h5", raw_matrix, type = "HDF5")

raw_matrix <- Read10X("/research/labs/neurology/fryer/m239830/Ecoli_pigs/snRNAseq/pipseeker_output/E2B/raw_matrix/")
write10xCounts("/research/labs/neurology/fryer/m239830/Ecoli_pigs/snRNAseq/pipseeker_output/E2B/raw_feature_bc_matrix.h5", raw_matrix, type = "HDF5")