#--- libraries
library(Seurat)
library(stringr)
library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)
library(Azimuth)
library(dittoSeq)
library(dplyr)
library(RColorBrewer)

#--- variables
Ecoli_colors <- c("red")
sample_colors <- c("red")

#--- references and metadata
metadata <-
  read.delim(
    "/research/labs/neurology/fryer/m239830/Ecoli_pigs/metadata.tsv")

pathToRef = c("/research/labs/neurology/fryer/projects/references/pig/ensembl_v7/")
pathToRawData = c("/research/labs/neurology/fryer/projects/sepsis/pig/Ecoli/PIPseq_snRNAseq/2024_Brain/usftp21.novogene.com")


#--- functions 
saveToPDF <- function(...) {
  d = dev.copy(pdf, ...)
  dev.off(d)
}
