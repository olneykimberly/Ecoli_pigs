knitr::opts_knit$set(root.dir = ".")
source("bulk_RNAseq/scripts/file_paths_and_colours.R")

# --- 2. Data Processing ---
bulk_data <- read_excel("bulk_RNAseq/results/DEGs/DEGs.q1.00.xlsx")

# Make bulk DEG table unique by gene so all downstream counts are gene-level
bulk_deg_df <- bulk_data %>%
  mutate(
    gene_name = as.character(gene_name),
    logFC = as.numeric(logFC),
    adj.P.Val = as.numeric(adj.P.Val)
  ) %>%
  filter(!is.na(gene_name), gene_name != "") %>%
  filter(!is.na(adj.P.Val), !is.na(logFC)) %>%
  filter(adj.P.Val < 0.05, abs(logFC) > 1) %>%
  transmute(
    gene = gene_name,
    bulk_logFC = logFC,
    bulk_direction = ifelse(logFC > 0, "Up", "Down")
  ) %>%
  distinct(gene, .keep_all = TRUE)

bulk_degs <- unique(bulk_deg_df$gene)

sn_path   <- "snRNAseq/results/DEGs/DESeq2_pseudobulk_exp_filter/Ecoli_vs_Saline.xlsx"
sn_sheets <- excel_sheets(sn_path)

# Keep one row per CellType-gene pair
sn_deg_df <- map_dfr(sn_sheets, function(sheet) {
  dat <- read_excel(sn_path, sheet = sheet) %>%
    mutate(
      gene = as.character(gene),
      log2FoldChange = suppressWarnings(as.numeric(log2FoldChange)),
      padj = suppressWarnings(as.numeric(padj))
    ) %>%
    filter(!is.na(gene), gene != "") %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    filter(padj < 0.05, abs(log2FoldChange) > 1)
  
  # If no DEGs in this sheet, return an empty tibble with fixed column types
  if (nrow(dat) == 0) {
    return(tibble(
      CellType = character(),
      gene = character(),
      sn_logFC = numeric(),
      sn_direction = character()
    ))
  }
  
  dat %>%
    transmute(
      CellType = as.character(sheet),
      gene = gene,
      sn_logFC = log2FoldChange,
      sn_direction = case_when(
        log2FoldChange > 0 ~ "Up",
        log2FoldChange < 0 ~ "Down",
        TRUE ~ NA_character_
      )
    ) %>%
    distinct(CellType, gene, .keep_all = TRUE)
})

sn_deg_list <- split(sn_deg_df$gene, sn_deg_df$CellType)
sn_union_degs <- unique(sn_deg_df$gene)

# --- 3. Gene Comparison Dataframes ---
df_shared <- data.frame(gene = intersect(bulk_degs, sn_union_degs), category = "Shared")
df_unique_bulk <- data.frame(gene = setdiff(bulk_degs, sn_union_degs), category = "Unique_Bulk")

# --- 4. Panel A: Venn Diagram ---
panel_a <- as.ggplot(expression(
  draw.pairwise.venn(
    area1 = length(unique(bulk_degs)),
    area2 = length(unique(sn_union_degs)),
    cross.area = length(intersect(unique(bulk_degs), unique(sn_union_degs))),
    category = c("Bulk", "snRNA"),
    fill = c("skyblue", "mediumorchid"),
    alpha = 0.5,
    lty = "blank",
    cex = 0.75,
    cat.cex = 0.7,
    cat.pos = c(-20, 20),
    margin = 0.15
  )
))
panel_a
# --- 5. NEW Panel B: Dominant cell type assignment (counts) ---
# For each bulk DEG, assign it to the cell type with the largest |sn_logFC|
# among same-direction overlaps. If no same-direction snRNA DEG exists,
# label it "Not recovered".

bulk_sn_overlap <- bulk_deg_df %>%
  inner_join(sn_deg_df, by = "gene") %>%
  mutate(same_direction = bulk_direction == sn_direction)

dominant_assignments <- bulk_sn_overlap %>%
  filter(same_direction) %>%
  mutate(abs_sn_logFC = abs(sn_logFC)) %>%
  arrange(gene, desc(abs_sn_logFC), CellType) %>%
  group_by(gene) %>%
  slice(1) %>%
  ungroup() %>%
  select(gene, CellType)

dominant_all <- bulk_deg_df %>%
  select(gene) %>%
  distinct(gene) %>%
  left_join(dominant_assignments, by = "gene") %>%
  mutate(CellType = ifelse(is.na(CellType), "Not recovered", CellType))

dominant_summary <- dominant_all %>%
  count(CellType, name = "n_bulk_genes")

# Put "Not recovered" at the end of the plot
celltype_order <- dominant_summary %>%
  filter(CellType != "Not recovered") %>%
  arrange(desc(n_bulk_genes)) %>%
  pull(CellType)

celltype_order <- c(celltype_order, "Not recovered")

NR_genes <- subset(dominant_all, CellType == "Not recovered") 
df_unique_bulk$gene

setdiff(NR_genes$gene, df_unique_bulk$gene)
dominant_summary <- dominant_summary %>%
  mutate(CellType = factor(CellType, levels = celltype_order))

panel_b <- ggplot(dominant_summary,
                  aes(x = n_bulk_genes , y = CellType)) +
  geom_col(fill = "indianred") +
  geom_text(aes(label = n_bulk_genes), vjust = -.15, size = 2.5) +
  coord_flip() +
  labs(
    title = "Dominant cell type\ndriving bulk DEG signal",
    x = "Number of bulk DEGs",
    y = "Cell type"
  ) +
  my_seurat_theme() +
  theme(
    axis.title.x = element_text(size = 8)
  )

panel_b

# --- 5.5 Panel C: DEG Counts by Cell Type (Bidirectional) ---
deg_counts <- map_dfr(sn_sheets, function(sheet) {
  res <- read_excel(sn_path, sheet = sheet) %>%
    mutate(
      log2FoldChange = as.numeric(log2FoldChange),
      padj = as.numeric(padj)
    ) %>%
    filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    filter(padj < 0.05 & abs(log2FoldChange) >= 1)
  
  if (nrow(res) > 0) {
    res %>%
      mutate(Direction = ifelse(log2FoldChange > 0, "Up", "Down")) %>%
      group_by(Direction) %>%
      tally() %>%
      mutate(CellType = sheet,
             Direction = as.character(Direction))
  } else {
    tibble(Direction = character(), n = integer(), CellType = sheet)
  }
})

deg_counts$CellType <- gsub(" ", "", deg_counts$CellType)

deg_counts_plot <- deg_counts %>%
  complete(CellType, Direction = c("Up", "Down"), fill = list(n = 0)) %>%
  mutate(n_plot = ifelse(Direction == "Down", -n, n))

my_cell_order <- c("endothelial", "mural", "VLMC", "microglia", "astrocyte", "oligodendrocyte", "opc", 
                   "neuron")

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
    legend.key.size = unit(0.3, "cm"),
    legend.text = element_text(size = 7),
    legend.spacing.x = unit(0.1, "cm"),
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
  ncol = 3,
  labels = c("C", "D", "E"),
  widths = c(1.15, .65, 1),
  font.label = list(size = 10)
)

Option1

path <- "../manuscript_figures/Supplemental_Figure_7_with_DEG_counts"
if(!dir.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE)
ggsave(paste0(path, ".pdf"), plot = Option1, width = 7, height = 3)

print("Figure saved and Dataframes (df_shared, df_unique_bulk) created.")

# --- 7. Output comparison tables ---
df_shared <- bulk_data %>% 
  filter(gene_name %in% intersect(bulk_degs, sn_union_degs)) %>%
  mutate(Comparison_Category = "Shared_Bulk_and_snRNA")

df_unique_bulk <- bulk_data %>% 
  filter(gene_name %in% setdiff(bulk_degs, sn_union_degs)) %>%
  mutate(Comparison_Category = "Unique_to_Bulk")

sn_deg_dfs <- map(sn_sheets, ~read_excel(sn_path, sheet = .x) %>% filter(padj < 0.05))
names(sn_deg_dfs) <- sn_sheets

df_unique_sn <- bind_rows(sn_deg_dfs, .id = "CellType") %>%
  filter(gene %in% setdiff(sn_union_degs, bulk_degs)) %>%
  mutate(Comparison_Category = "Unique_to_snRNA")

list_of_datasets <- list(
  "Shared" = df_shared,
  "Unique_to_bulkRNA" = df_unique_bulk,
  "Unique_to_snRNA" = df_unique_sn
)

list_of_datasets <- c(list_of_datasets, sn_deg_dfs)

write.xlsx(list_of_datasets, file = "manuscript_figures/Sepsis_Pig_Transcriptome_Full_Comparison.xlsx")