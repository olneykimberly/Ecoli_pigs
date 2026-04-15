#----------------- Libraries
# R version 4.3 was used in for this analysis
# R/x86_64-pc-linux-gnu-library/4.3
.libPaths(c("/tgen_labs/jfryer/kolney/R/x86_64-pc-linux-gnu-library/4.3", "/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))

# Tidyverse & General Data Manipulation
library(tidyverse)      # Includes ggplot2, dplyr, stringr, forcats, etc.
library(data.table)
library(plyr)
library(reshape2)
library(R.utils)
library(car)

# Bioconductor & Bioinformatics
library(DESeq2)
library(edgeR)
library(limma)
library(glmGamPoi)
library(GenomicFeatures)
library(BiocParallel)
library(remaCor)
library(dittoSeq)
library(variancePartition)

# Visualization & Graphics
library(ggplot2)
library(dplyr)
library(forcats)
library(scales)
library(ggrepel)
library(gplots)
library(RColorBrewer)
library(gridExtra)
library(grid)
library(lattice)
library(grDevices)
library(purrr)
library(reshape)
library(ggforce)
library("readxl")
library(patchwork) 
library(ggpubr)
library(gridExtra)
library(ggplotify)
library(NatParksPalettes)
color.panel <- dittoColors()

# Utilities & Development
library(openxlsx)
library(devtools)

# Single nucleus 
library(SeuratObject)
library(Signac)
library(Seurat) 
library(harmony)
library(SeuratWrappers)
library(dittoSeq)
library(DoubletFinder)


#----------------- Create folder outputs
# Create directory structure if it doesn't exist
# This ensures that saveRDS() and saveToPDF() calls don't fail
required_dirs <- c("bulk_RNAseq/results", 
                   "bulk_RNAseq/rObjects",
                   "bulk_RNAseq/results/MDS",
                   "bulk_RNAseq/results/library",
                   "bulk_RNAseq/results/voom",
                   "bulk_RNAseq/results/variance",
                   "bulk_RNAseq/results/DEGs",
                   "bulk_RNAseq/results/PCA",
                   "bulk_RNAseq/results/volcano",
                   "bulk_RNAseq/results/gprofiler",
                   "bulk_RNAseq/rObjects/gene_tables", 
                   "snRNAseq/results/CellChat",
                   "snRNAseq/results/CellChat_Ecoli_vs_Saline",
                    "snRNAseq/results/DEGs",
                    "snRNAseq/results/DoubletFinder",
                    "snRNAseq/results/UMAP",
                    "snRNAseq/results/UpSet",
                    "snRNAseq/results/density",
                    "snRNAseq/results/dot_plot",
                    "snRNAseq/results/feature",
                    "snRNAseq/results/gprofiler",
                    "snRNAseq/results/markers",
                    "snRNAseq/results/metascape",
                    "snRNAseq/results/nuclei_count",
                    "snRNAseq/results/pca",
                    "snRNAseq/results/scatter",
                    "snRNAseq/results/top_transcripts",
                    "snRNAseq/results/tree",
                    "snRNAseq/results/violin",
                    "snRNAseq/results/volcanoes")

# Loop through and create them
lapply(required_dirs, function(x) if(!dir.exists(x)) dir.create(x, recursive = TRUE))

#----------------- Define variables
typeOfCount <- c("STAR.bamReadsPerGene.out.tab") 
pathToRef <- c("../../projects/references/pig/ensembl_v7/") # out of Ecoli, out of home, into projects
#----------------- Data
# Both LPS and Ecoli pig project information combined into a single master metadata file 
metadata <- read.delim("Ecoli_pigs_bulkRNAseq.tsv", header = TRUE, sep = "\t")
#----------------- Functions
saveToPDF <- function(...) {
  d = dev.copy(pdf,...)
  dev.off(d)
}

# for creating interscetion tables and retaining the gene name
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

my_seurat_theme <- function() {
  theme_classic() + # Start with a base theme
    theme(
    axis.title.x = element_text(size = 8),
    axis.text.x = element_text(size = 8),
    axis.title.y = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.title = element_text(size = 8), 
    legend.text = element_text(size = 8),  
    plot.title = element_text(size = 8), 
    theme(legend.position = "none")
    )
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
      guide = "colourbar", breaks = c(-4, 0, 4, 8, 10),
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
