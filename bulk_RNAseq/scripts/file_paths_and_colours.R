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
library(gprofiler2)
library(ggpubr)
library(forcats)
library(stringr)
library(reshape2)

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
#rm(brain_metadata)

metadata_extra <- read.delim("/tgen_labs/jfryer/kolney/Ecoli_pigs/metadata_extra_clinical.txt", header = TRUE, sep = "\t")
metadata <- merge(metadata, metadata_extra, by = "pig_id")

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

#-----------------------
# function
publish_gostplot_intersect <- function (p, highlight_terms = NULL, filename = NULL, width = NA, 
                                        height = NA) 
{
  if (!("ggplot" %in% class(p))) {
    warning("Highlighting terms in a Manhattan plot is available for a ggplot object only.\nPlease set 'interactive = F' in the gostplot() function and try again.")
    return(NULL)
  }
  term_id <- logpval <- term_size_scaled <- id <- query <- p_value <- NULL
  if (!is.null(highlight_terms)) {
    if (is.data.frame(highlight_terms)) {
      message("The input 'highlight_terms' is a data.frame and therefore the column 'term_id' will be used for detection.")
      if ("term_id" %in% colnames(highlight_terms)) {
        highlight_terms <- highlight_terms$term_id
      }
      else {
        stop("No column named 'term_id'.")
      }
    }
    df <- p$data
    subdf <- base::subset(df, term_id %in% highlight_terms)
    if (nrow(subdf) == 0) {
      message("None of the term IDs in the 'highlight_terms' was found from the results.")
      return(p)
    }
    highlight_terms <- unique(highlight_terms)
    subdf$id <- match(subdf$term_id, highlight_terms)
    p <- p + ggplot2::geom_point(data = subdf, ggplot2::aes(x = order, 
                                                            y = logpval, size = term_size_scaled), pch = 21, 
                                 colour = "black")
    p <- p + ggplot2::geom_text(data = subdf, size = 4, colour = "white", 
                                ggplot2::aes(label = as.character(id), family = "mono", 
                                             fontface = "bold"), hjust = -1.2, vjust = -0.05) + 
      ggplot2::geom_text(data = subdf, size = 4, colour = "black", 
                         fontface = "bold", ggplot2::aes(label = as.character(id)), 
                         hjust = -1.2, vjust = -0.05)
    pseudo_gostres <- list(result = data.frame(subdf), meta = list(query_metadata = list(queries = sapply(unique(subdf$query), 
                                                                                                          function(x) NULL))))
    tb <- publish_gosttable(pseudo_gostres, highlight_terms = highlight_terms, 
                            use_colors = TRUE, show_columns = c("source", "term_name", 
                                                                "intersection_size"), filename = NULL, ggplot = FALSE)
    h <- grid::unit.c(grid::unit(1, "null"), sum(tb$heights) + 
                        grid::unit(3, "mm"))
    w <- grid::unit.c(grid::unit(1, "null"))
    tg <- gridExtra::grid.arrange(p, tb, ncol = 1, heights = h, 
                                  widths = w, newpage = TRUE)
    p <- ggplot2::ggplot() + ggplot2::annotation_custom(tg) + 
      ggplot2::geom_blank() + ggplot2::theme_void()
  }
  if (is.null(filename)) {
    return(p)
  }
  else {
    imgtype <- strsplit(basename(filename), split = "\\.")[[1]][-1]
    if (length(imgtype) == 0) {
      filename = paste0(filename, ".pdf")
    }
    if (tolower(imgtype) %in% c("png", "pdf", "jpeg", "tiff", 
                                "bmp")) {
      if (is.na(width)) {
        width = max(grDevices::dev.size()[1], 8)
      }
      if (is.na(height)) {
        height = max(grDevices::dev.size()[2], 6)
      }
      ggplot2::ggsave(filename = filename, plot = p, width = width, 
                      height = height, limitsize = F)
      message("The image is saved to ", filename)
      return(p)
    }
    else {
      stop("The given file format is not supported.\nPlease use one of the following extensions: .png, .pdf, .jpeg, .tiff, .bmp")
    }
  }
}


fromList <- function (input) {
  # Same as original fromList()...
  elements <- unique(unlist(input))
  data <- unlist(lapply(input, function(x) {
    x <- as.vector(match(elements, x))
  }))
  data[is.na(data)] <- as.integer(0)
  data[data != 0] <- as.integer(1)
  data <- data.frame(matrix(data, ncol = length(input), byrow = F))
  data <- data[which(rowSums(data) != 0), ]
  names(data) <- names(input)
  # ... Except now it conserves your original value names!
  row.names(data) <- elements
  return(data)
}

generate_heatmap <- function(go_term, gene_list, DEG_df, top_df) {
  # Filter top genes for this GO term
  top_genes <- top_df %>% filter(GO_term == go_term)
  gene_ids <- top_genes$gene_id
  
  # Get relevant DEGs
  DEG_subset <- DEG_df %>% filter(gene_id %in% gene_ids)
  DEG_merged <- merge(DEG_subset, top_genes, by = "gene_id", all.x = TRUE)
  
  DEG_merged$gene <- factor(DEG_merged$gene, levels = unique(top_genes$gene))
  DEG_merged$gene <- fct_rev(DEG_merged$gene)
  
  # Plot
  heatmap_plot <- ggplot(data = DEG_merged, aes(x = samples, y = gene)) +
    geom_tile(aes(fill = counts)) +
    facet_grid(~ Condition, scales = "free", switch = "both") +
    scale_fill_gradient2(
      low = "#FFFFCCFF", mid = "#FD8D3CFF", high = "#800026FF", midpoint = 2.5,
      space = "rgb", guide = "colourbar", breaks = c(-4, 0, 4, 8, 10),
      name = expression(log[2](CPM))
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 10),
      legend.position = "none",
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.margin = margin(0, 0.2, 0, 0.2, "cm"),
      panel.spacing = unit(0, 'lines'),
      plot.title = element_text(size = 12, vjust = -1, hjust = 0.5),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size = 10, margin = margin(r = -2)),
      axis.ticks.y = element_blank(), 
    ) +
    ggtitle(go_term)
  
  return(heatmap_plot)
}

addSmallLegend <- function(myPlot, pointSize = 6, textSize = 10, spaceLegend = .5) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

addSmallLegend_UMAP <- function(myPlot, pointSize = 6, textSize = 8, spaceLegend = .5) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
