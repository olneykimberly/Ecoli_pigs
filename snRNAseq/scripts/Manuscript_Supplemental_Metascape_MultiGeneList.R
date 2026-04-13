library(ggplot2)
library(reshape2)
library(dplyr)
library(forcats)
library(scales)
setwd("/tgen_labs/jfryer/kolney/Ecoli_pigs/bulk_RNAseq/scripts")

# --- 1. Setup and Themes ---

# Define the user-provided theme
my_seurat_theme <- function() {
  theme_classic() + 
    theme(
      axis.title.x = element_text(size = 8),
      axis.text.x = element_text(size = 8, angle = 270, vjust = 0.5, hjust = 0),
      axis.title.y = element_text(size = 8),
      axis.text.y = element_text(size = 8),
      legend.title = element_text(size = 8), 
      legend.text = element_text(size = 8),  
      plot.title = element_text(size = 8),
      legend.position = "right"
    )
}

# Helper to shrink legend size
addSmallLegend <- function(my_plot, pointSize = 3, textSize = 6, spaceLegend = .5) {
  my_plot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}

# Define the consistent cell type order
cell_type_order <- c("microglia","astrocyte", "oligodendrocyte", "opc", "endothelial", "mural", "VLMC")

# Function to process Metascape Heatmap CSVs
process_metascape_heatmap <- function(file_path, direction_suffix) {
  df <- read.csv(file_path)
  
  # Set Description as factor to preserve row order from Metascape
  df$Description <- factor(df$Description, levels = rev(df$Description))
  
  # Melt for ggplot
  df_melted <- melt(df, id.vars = "Description")
  
  # Clean column names
  df_melted$variable <- gsub(paste0("X_LogP_|_", direction_suffix), "", df_melted$variable)
  
  # Filter and set factor order
  df_melted <- df_melted %>% filter(variable %in% cell_type_order)
  df_melted$variable <- factor(df_melted$variable, levels = cell_type_order)
  
  # Logic to handle the "0 is grey" requirement
  df_melted$value <- as.numeric(df_melted$value)
  df_melted$plot_value <- ifelse(df_melted$value == 0, NA, df_melted$value)
  
  return(df_melted)
}

# --- 2. Upregulated Heatmap (A) ---

up_path <- "/tgen_labs/jfryer/kolney/Ecoli_pigs/snRNAseq/results/metascape/multi_input_up/Enrichment_heatmap/HeatmapSelectedGO.csv"
up_data <- process_metascape_heatmap(up_path, "up")

heatmap_up <- ggplot(up_data, aes(variable, Description, fill = plot_value)) + 
  geom_tile(color = "white") +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "", y = "", title = "Upregulated") +
  scale_fill_gradientn(
    colours = c("#800026FF", "#FD8D3CFF", "#FFFFCCFF"),
    name = expression(log[10](p)),
    limits = c(-25, 0), 
    oob = scales::squish, 
    na.value = "grey80" 
  ) + 
  my_seurat_theme()

heatmap_up <- addSmallLegend(heatmap_up)

# --- 3. Downregulated Heatmap (B) ---

down_path <- "/tgen_labs/jfryer/kolney/Ecoli_pigs/snRNAseq/results/metascape/multi_input_down/Enrichment_heatmap/HeatmapSelectedGO.csv"
down_data <- process_metascape_heatmap(down_path, "down")

heatmap_down <- ggplot(down_data, aes(variable, Description, fill = plot_value)) + 
  geom_tile(color = "white") +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "", y = "", title = "Downregulated") +
  scale_fill_gradientn(
    # Keeping same color scheme as UP for consistency
    colours = c("#800026FF", "#FD8D3CFF", "#FFFFCCFF"),
    name = expression(log[10](p)),
    
    # Adjusted limits for down if needed, or keep -25 for direct comparison
    limits = c(-25, 0), 
    
    oob = scales::squish, 
    na.value = "grey80" 
  ) + 
  my_seurat_theme()

heatmap_down <- addSmallLegend(heatmap_down)

# --- 4. Display/Print ---

print(heatmap_up)
print(heatmap_down)


supplementalFigure <-
  ggarrange(
    print(heatmap_up),
    print(heatmap_down),
    nrow = 2,
    labels = c("A","B"), 
    font.label = list(size = 10)
    )
supplementalFigure
path <- paste0("../results/manuscript_figures/Supplemental_figure_6_celltypes_metascape")
saveToPDF(paste0(path, ".pdf"), width = 7, height = 9)