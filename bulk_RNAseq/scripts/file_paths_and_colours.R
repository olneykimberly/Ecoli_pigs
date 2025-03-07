#----------------- Libraries
.libPaths(c("/tgen_labs/jfryer/kolney/R/rstudio-4.3.0-4-with_modules.sif/libs", "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
.libPaths()
library(Matrix, lib.loc = "/usr/local/lib/R/site-library")
library(SeuratObject)
library(Signac)
library(Seurat) 
library(stringr)
library(ggplot2)
library(harmony)
library(remaCor)
library(gridExtra)
library(grid)
library(lattice)
library(R.utils)
library(SeuratWrappers)
library(Azimuth)
library(dittoSeq)
library(dplyr)
library(RColorBrewer)
library(DESeq2) # adds matrix
require(openxlsx)
library(ggrepel)
library(glmGamPoi)
library(devtools)
library(harmony)
library(DoubletFinder)
library(reshape2)
library(ggtree)
library(BiocParallel) 
library(edgeR)  
library(limma)  
library(ggrepel) 
library(ggplot2) 
library(gplots) 
library(grDevices)  
#library(philentropy) 
library(stringr) 
library(remaCor)
library(scales)
library(tximport)
library(tidyverse)
library(GenomicFeatures)
library(dplyr)
library(plyr)
library(gridExtra)
library(grid)
library(lattice)
library(data.table)
library(openxlsx)

#BiocManager::install("variancePartition", lib="/tgen_labs/jfryer/kolney/R/x86_64-pc-linux-gnu-library/4.3")
#BiocManager::install("variancePartition")
#devtools::install_github("DiseaseNeuroGenomics/variancePartition")
#BiocManager::install("limma", force = TRUE) 

#----------------- Define variables
tissue <- c("Brain") # Kidney or Brain
#control_color <- "gray29"
#Ecoli_color <- "green"
#High_Ecoli_color <- "darkgreen"
#Low_Ecoli_color <- "lightgreen"
#LPS_color <- "purple"
#output_dir <- c("") 
typeOfCount <- c("STAR.bamReadsPerGene.out.tab") 
pathToRef <- c("/tgen_labs/jfryer/projects/references/pig/ensembl_v7/")

#----------------- Functions
saveToPDF <- function(...) {
  d = dev.copy(pdf,...)
  dev.off(d)
}

#----------------- Data
# Both LPS and Ecoli pig project information combined into a single master metadata file 
metadata <- read.delim("/tgen_labs/jfryer/kolney/Ecoli_pigs/all_pigs_combined_metadata.tsv", header = TRUE, sep = "\t")
# Update path to star counts
metadata$path <- gsub("/research/labs/neurology/fryer/m239830/Ecoli_pigs/bulk_RNAseq/starAligned/", 
                   "/tgen_labs/jfryer/kolney/Ecoli_pigs/bulk_RNAseq/starAligned/", 
                   metadata$path)
# Remove pigs 9 & 13 from the LPS project, they died quickly and we're excluded 
metadata <- metadata[metadata$group != "LPS",]
metadata <- metadata[metadata$group != "Control",] # The saline samples from the LPS study

# We will only focus on the brain of saline and high dose samples
metadata <- metadata[metadata$condition_dose != "Low_dose_Ecoli",] # The saline samples from the LPS study
metadata <- metadata[metadata$flowcell != "HHKHJDRXY",] # The saline samples from the LPS study

# Remove samples that are excluded from the analysis. We only have single nuclues data for 4 saline and 4 Ecoli samples. 
# We will limit our analysis to those sample pigs. 
# Keep pigs: Saline 1, 2, 3, 6. Ecoli: 1, 4, 6, 8
brain_metadata <- metadata[metadata$tissue == "Brain", ]
keep_pigs <- c("S1", "S2", "S3", "S6", "E1", "E4", "E6", "E8")
brain_metadata <- brain_metadata[brain_metadata$pig_id %in% keep_pigs, ]
metadata <- brain_metadata
metadata$sample_ID <- paste0(metadata$sample, "_1_BR")
# clean up
rm(brain_metadata)

#------------------------------- Ecoli pig snRNAseq info
#sn_metadata <- read.delim("/tgen_labs/jfryer/projects/sepsis/pig/Ecoli/Ecoli_pig_snRNA_seq_info.txt", header = TRUE, sep = "\t")


#------------------------------- Ecoli and LPS seperate metadata files 
# read in Ecoli metadata
# Ecoli_meta <-
#   read.delim((
#     "/research/labs/neurology/fryer/projects/sepsis/pig/Ecoli/metadata.tsv"
#   ),
#   header = TRUE,
#   sep = "\t"
#   )
# Ecoli_meta <- Ecoli_meta[ -c(1)] # remove first column 
# Ecoli_meta$sample_name <- gsub("\\..*", "", Ecoli_meta$filename) # create sample_name column
# # create lane column 
# lane <- str_sub(Ecoli_meta$run_flowcell_lane,start=-1) 
# Ecoli_meta$lane <- paste0("L", lane)
# rm(lane)
# Ecoli_meta <- Ecoli_meta[Ecoli_meta$tissue == tissue, ]

# Read data LPS data
# read in metadata
# LPS_meta <-
#   read.delim((
#     "/research/labs/neurology/fryer/projects/sepsis/pig/LPS/metadata.tsv"
#   ),
#   header = TRUE,
#   sep = "\t"
#   )
# # subset for tissue 
# LPS_meta <- LPS_meta[LPS_meta$tissue == tissue, ]
# 
# # remove pigs 9 and 13
# LPS_meta <- LPS_meta[LPS_meta$pig_id != "9", ]
# LPS_meta <- LPS_meta[LPS_meta$pig_id != "13", ]
# LPS_meta <- LPS_meta[LPS_meta$blood_group != "BB", ]
# 
# # path to counts files
# LPS_count_files <-
#   file.path(paste0(
#     "/research/labs/neurology/fryer/m239830/Ecoli_pigs/bulk_RNAseq/LPS_pigs/starAligned/",
#     LPS_meta$featureCounts_name,"_",
#     typeOfCount
#   ))
# # add sample name to counts files
# names(LPS_count_files) <- paste0(LPS_meta$featureCounts_name)
# 
# 
# # sleuth and other tools requires path, sample and condition columns.
# # add this information to metadata
# LPS_meta$path <- LPS_count_files
# LPS_meta$sample <- LPS_meta$simplified_name
# LPS_meta$condition <- as.factor(LPS_meta$group)


