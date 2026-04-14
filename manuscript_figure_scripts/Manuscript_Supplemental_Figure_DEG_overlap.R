# 1. Install/Load required for conversion if not present
library(ggplotify)
library(readxl)
library(dplyr)
library(purrr)
library(VennDiagram)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(readxl)
library(openxlsx)

getwd()

my_seurat_theme <- function() {
  theme_classic() + 
    theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_text(size = 8, angle = 45, hjust = 1),
      axis.title.y = element_text(size = 8),
      axis.text.y  = element_text(size = 8),
      legend.title = element_text(size = 8), 
      legend.text  = element_text(size = 8),  
      plot.title   = element_text(size = 8),
      legend.position = "none"
    )
}

# --- Bulk DEGs ---
bulk_data <- read_excel("/tgen_labs/jfryer/kolney/Ecoli_pigs/bulk_RNAseq/results/DEGs/DEGs.q1.00.xlsx")

bulk_deg_df <- bulk_data %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
  transmute(
    gene = gene_name,
    bulk_logFC = logFC,
    bulk_direction = ifelse(logFC > 0, "Up", "Down")
  ) %>%
  distinct()

bulk_degs <- bulk_deg_df$gene

# --- snRNA DEGs by cell type ---
sn_path   <- "/tgen_labs/jfryer/kolney/Ecoli_pigs/snRNAseq/results/DEGs/DESeq2_pseudobulk_exp_filter/Ecoli_vs_Saline.xlsx"
sn_sheets <- excel_sheets(sn_path)

sn_deg_df <- map_dfr(sn_sheets, function(sheet) {
  read_excel(sn_path, sheet = sheet) %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1) %>%
    transmute(
      CellType = sheet,
      gene = gene,
      sn_logFC = log2FoldChange,
      sn_direction = ifelse(log2FoldChange > 0, "Up", "Down")
    ) %>%
    distinct()
})

# --- Panel B data: overlap of each cell type with bulk ---
fraction_data <- sn_deg_df %>%
  inner_join(bulk_deg_df, by = "gene") %>%
  group_by(CellType) %>%
  summarise(
    n_shared = n_distinct(gene),
    pct_bulk = 100 * n_distinct(gene) / length(unique(bulk_degs)),
    pct_celltype = 100 * n_distinct(gene) / n_distinct(sn_deg_df$gene[sn_deg_df$CellType == first(CellType)]),
    n_same_direction = sum(sn_direction == bulk_direction, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    concordance = 100 * n_same_direction / n_shared
  ) %>%
  arrange(desc(n_shared))

fraction_data_concordant <- sn_deg_df %>%
  inner_join(bulk_deg_df, by = "gene") %>%
  filter(sn_direction == bulk_direction) %>%
  group_by(CellType) %>%
  summarise(
    n_concordant = n_distinct(gene),
    .groups = "drop"
  ) %>%
  arrange(desc(n_concordant))

panel_b_concordant <- ggplot(fraction_data_concordant,
                             aes(x = reorder(CellType, n_concordant), y = n_concordant)) +
  geom_col(fill = "firebrick3") +
  geom_text(aes(label = n_concordant), hjust = -0.1, size = 2.5) +
  coord_flip() +
  labs(
    title = "Concordant bulk DEGs recovered in each cell type",
    x = "Cell Type",
    y = "Number of concordant shared DEGs"
  ) +
  my_seurat_theme() +
  theme(axis.title.x = element_text(size = 8))

panel_b_concordant
# --- 5.5 Panel C: DEG Counts by Cell Type (Bidirectional) ---
# Process snRNA sheets to get Up/Down counts per cell type
deg_counts <- map_dfr(sn_sheets, function(sheet) {
  # Read and filter
  res <- read_excel(sn_path, sheet = sheet) %>%
    filter(padj < 0.05 & abs(log2FoldChange) >= 1)
  
  # Check if any DEGs exist to avoid type errors
  if (nrow(res) > 0) {
    res %>%
      mutate(Direction = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
      group_by(Direction) %>%
      tally() %>%
      mutate(CellType = sheet,
             Direction = as.character(Direction)) # Force character type
  } else {
    # Return an empty tibble with correct column types if no DEGs found
    tibble(Direction = character(), n = integer(), CellType = sheet)
  }
})

deg_counts$CellType <- gsub(" ", "", deg_counts$CellType)
# Complete the data: If a cell type has "Up" but no "Down" (or vice versa), 
# this ensures the plotting data is balanced.
deg_counts_plot <- deg_counts %>%
  complete(CellType, Direction = c("Up", "Down"), fill = list(n = 0)) %>%
  mutate(n_plot = ifelse(Direction == "Down", -n, n))


my_cell_order <- c("endothelial", "mural", "VLMC", "microglia", "astrocyte", "oligodendrocyte", "opc", 
                   "neuron")

# 2. Convert to factor
deg_counts_plot <- deg_counts_plot %>%
  mutate(CellType = factor(CellType, levels = rev(my_cell_order)))


panel_c <- ggplot(deg_counts_plot, aes(x = CellType, y = n_plot, fill = Direction)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n), 
            vjust = 0.5, 
            hjust = ifelse(deg_counts_plot$n_plot > 0, -0.2, 1.2),
            size = 2.) +
  scale_fill_manual(values = c("Down" = "blue", "Up" = "red")) +
  coord_flip() + 
  geom_hline(yintercept = 0, linetype = "solid", color = "black") +
  labs(title = "Number of DEGs within each cell type\nq < 0.05 & |log2FC| > 1",
       y = "", x = "") +
  my_seurat_theme() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    # --- Adjust Legend Box Size Here ---
    legend.key.size = unit(0.3, "cm"),      # Makes the boxes smaller
    legend.text = element_text(size = 7),    # Optional: shrink text to match
    legend.spacing.x = unit(0.1, "cm"),      # Keeps text close to the box
    # -----------------------------------
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    plot.title = element_text(size = 8),
    plot.margin = margin(.1, .35, .01, 0, "cm")

  ) +
  scale_y_continuous(limits = c(min(deg_counts_plot$n_plot) * 1.3, 
                                max(deg_counts_plot$n_plot) * 1.3))

panel_c

# --- 6. Combine and Save ---
Option1 <- ggarrange(
  panel_c, 
  panel_a, 
  panel_b,
  ncol = 3, # Side-by-side
  labels = c("C", "D", "E"),
  widths = c(1.15, .65, 1),
  font.label = list(size = 10)
)
Option1
# Custom save logic

path <- "../manuscript_figures/Supplemental_Figure_7_with_DEG_counts"
if(!dir.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE)
ggsave(paste0(path, ".pdf"), plot = Option1, width = 7, height = 3)

print("Figure saved and Dataframes (df_shared, df_unique_bulk) created.")

df_shared <- bulk_data %>% 
  filter(gene_name %in% intersect(bulk_degs, sn_union_degs)) %>%
  mutate(Comparison_Category = "Shared_Bulk_and_snRNA")

# B) UNIQUE TO BULK: Significant in bulk, but not found in any snRNA cell type
# Often represents signal from rare cell types or global low-level shifts
df_unique_bulk <- bulk_data %>% 
  filter(gene_name %in% setdiff(bulk_degs, sn_union_degs)) %>%
  mutate(Comparison_Category = "Unique_to_Bulk")

# C) UNIQUE TO snRNA: Significant in at least one cell type, but not in bulk
# These are your high-resolution, cell-specific signals
# We will create a combined dataframe of these from the snRNA results
# Load all snRNA tabs into a named list of dataframes
sn_deg_dfs <- map(sn_sheets, ~read_excel(sn_path, sheet = .x) %>% filter(padj < 0.05))
names(sn_deg_dfs) <- sn_sheets
df_unique_sn <- bind_rows(sn_deg_dfs, .id = "CellType") %>%
  filter(gene %in% setdiff(sn_union_degs, bulk_degs)) %>%
  mutate(Comparison_Category = "Unique_to_snRNA")

list_of_datasets <- list(
  "Shared"     = df_shared,
  "Unique_to_bulkRNA" = df_unique_bulk,
  "Unique_to_snRNA"   = df_unique_sn
)

# Append the raw filtered DEGs for each cell type as individual tabs
list_of_datasets <- c(list_of_datasets, sn_deg_dfs)

write.xlsx(list_of_datasets, file = "Sepsis_Pig_Transcriptome_Full_Comparison.xlsx")
