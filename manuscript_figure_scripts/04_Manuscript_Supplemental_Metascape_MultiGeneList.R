source("bulk_RNAseq/scripts/file_paths_and_colours.R")
projectID <- "Ecoli_pigs_snRNA_object"


# --- 1. Setup and Themes ---
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

up_path <- "snRNAseq/results/metascape/multi_input_up/Enrichment_heatmap/HeatmapSelectedGO.csv"
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

down_path <- "snRNAseq/results/metascape/multi_input_down/Enrichment_heatmap/HeatmapSelectedGO.csv"
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
path <- paste0("manuscript_figures/Supplemental_figure_6_celltypes_metascape")
saveToPDF(paste0(path, ".pdf"), width = 7, height = 9)