#----------------- Libraries
library(BiocParallel) 
library(edgeR)  
library(limma)  
library(ggrepel) 
library(ggplot2) 
library(gplots) 
library(grDevices)  
library(philentropy) 
library(stringr) 
library(variancePartition) 
library(tximport)
library(tidyverse)
library(GenomicFeatures)
#library(tximportData)
# library(wasabi)
library(sleuth)
library(dplyr)
library(plyr)
library(gridExtra)
library(grid)
library(lattice)
library(data.table)

#----------------- Define variables
tissue <- c("Brain") # Kidney or Brain
#control <- "Saline"
#treatment <- "Ecoli"
control_color <- "gray29"
Ecoli_color <- "green"
LPS_color <- "purple"
output_dir <- c("") 
typeOfCount <- c("STAR.bamReadsPerGene.out.tab") 
pathToRef <- c("/research/labs/neurology/fryer/projects/references/pig/ensembl_v7/")

#----------------- Data
# read in Ecoli metadata
Ecoli_meta <-
  read.delim((
    "/research/labs/neurology/fryer/projects/sepsis/pig/Ecoli/metadata.tsv"
  ),
  header = TRUE,
  sep = "\t"
  )
Ecoli_meta <- Ecoli_meta[ -c(1)] # remove first column 
Ecoli_meta$sample_name <- gsub("\\..*", "", Ecoli_meta$filename) # create sample_name column
# create lane column 
lane <- str_sub(Ecoli_meta$run_flowcell_lane,start=-1) 
Ecoli_meta$lane <- paste0("L", lane)
rm(lane)
Ecoli_meta <- Ecoli_meta[Ecoli_meta$tissue == tissue, ]

# path to counts files
Ecoli_count_files <-
  file.path(paste0(
    "../starAligned/",
    Ecoli_meta$sample_name,
    "_", Ecoli_meta$lane, "_",
    typeOfCount
  ))
# add sample name to counts files
names(Ecoli_count_files) <- paste0(Ecoli_meta$sample)


# sleuth and other tools requires path, sample and condition columns.
# add this information to metadata
Ecoli_meta$path <- Ecoli_count_files
Ecoli_meta$sample <- Ecoli_meta$pig_id
Ecoli_meta$condition <- as.factor(Ecoli_meta$group)


# Read data LPS data
# read in metadata
LPS_meta <-
  read.delim((
    "/research/labs/neurology/fryer/projects/sepsis/pig/LPS/metadata.tsv"
  ),
  header = TRUE,
  sep = "\t"
  )
# subset for tissue 
LPS_meta <- LPS_meta[LPS_meta$tissue == tissue, ]

# remove pigs 9 and 13
LPS_meta <- LPS_meta[LPS_meta$pig_id != "9", ]
LPS_meta <- LPS_meta[LPS_meta$pig_id != "13", ]
LPS_meta <- LPS_meta[LPS_meta$blood_group != "BB", ]

# path to counts files
LPS_count_files <-
  file.path(paste0(
    "/research/labs/neurology/fryer/m239830/LPS_pigs/bulkRNA/starAligned/",
    LPS_meta$featureCounts_name,"_",
    typeOfCount
  ))
# add sample name to counts files
names(LPS_count_files) <- paste0(LPS_meta$featureCounts_name)


# sleuth and other tools requires path, sample and condition columns.
# add this information to metadata
LPS_meta$path <- LPS_count_files
LPS_meta$sample <- LPS_meta$simplified_name
LPS_meta$condition <- as.factor(LPS_meta$group)


#----------------- Functions
saveToPDF <- function(...) {
  d = dev.copy(pdf,...)
  dev.off(d)
}