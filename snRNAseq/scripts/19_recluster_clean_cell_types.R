knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/Ecoli_pigs/snRNAseq/scripts/")
setwd("/tgen_labs/jfryer/kolney/Ecoli_pigs/snRNAseq/scripts/")


source(here::here("/tgen_labs/jfryer/kolney/Ecoli_pigs/bulk_RNAseq/scripts/", "file_paths_and_colours.R"))
projectID <- "pigs_cellbender_after_recluster_harm_int_noise_removed_after_annotation_FINAL"
color.panel <- dittoColors()

dataObject <- readRDS(paste0("../rObjects/",projectID,"_annotated.rds"))

# Is PLP1 expression showing up slightly in all cell types, primarily driven by sample S1?
library(patchwork)
# UMAP showing the expression of select features
umap_feature <-
  FeaturePlot(dataObject,
              reduction = "umap",
              split.by = "sample",
              features = c("PLP1"))+
  plot_layout(nrow = 2, ncol = 4) 
umap_feature

#----------------------
# Re-look at the neuron cluster annotation. Maybe the neurons should be annotated into two different types of neuron annotations.
dataObject.neurons <- subset(dataObject, cells = WhichCells(dataObject, idents = "neuron"))
# neuron makers
neuron_markers <- c("AIF1",
                    "RBFOX3", 
                    "SNAP25",
                    "SYT1")

interneuron_markers <- c("GAD2",
                         "GAD1")

FeaturePlot(dataObject.neurons,
            reduction = "umap",
           # split.by = "sample",
            features = c("GAD2")) #+ plot_layout(nrow = 2, ncol = 4) 
help(FeaturePlot)
dittoDimPlot(object = dataObject.neurons,
                           var = "group",
                           reduction.use = "umap",
                           do.label = TRUE,
                           order = c("randomize"),
                           labels.highlight = TRUE)


dataObject.neurons <- RunAzimuth(dataObject.neurons, reference = "humancortexref")
# currently, the object has two layers in the RNA assay: counts, and data
# DimPlot(dataObject, group.by = "predicted.subclass", label = TRUE, label.size = 3) + NoLegend()
dittoDimPlot(object = dataObject.neurons,
                           var="predicted.subclass",
                           reduction.use = "umap",
                           do.label = FALSE,
                           labels.highlight = FALSE)



#-----------------------------------------------
# inspect
dataObject
Layers(dataObject)
# List of cell types to analyze
cell_types <- c(unique(dataObject$individual_clusters))

# Split by cell type 
object_list <- SplitObject(dataObject, split.by = "individual_clusters")

## ----neuron-------------------------------------------------------------------------------------------------------
dataObject.neuron <- SCTransform(object_list$neuron, verbose = FALSE)

# run PCA on the merged object
dataObject.neuron<- RunPCA(object = dataObject.neuron)
Idents(dataObject.neuron) <- "sample"

# Determine the K-nearest neighbor graph
dataObject.neuron <- FindNeighbors(object = dataObject.neuron,
                                   assay = "SCT",
                                   reduction = "pca",
                                   dims = 1:15)
# Run UMAP
dataObject.neuron <- RunUMAP(dataObject.neuron,
                             dims = 1:15,
                             reduction = "pca",
                             n.components = 3)

# Determine the clusters for various resolutions
dataObject.neuron <- FindClusters(object = dataObject.neuron,
                                  algorithm = 1, # 1= Louvain
                                  resolution = 0.7)

Idents(dataObject.neuron) <- dataObject.neuron$SCT_snn_res.0.7
dataObject.neuron$seurat_clusters <- dataObject.neuron$SCT_snn_res.0.7
ditto_umap <- dittoDimPlot(object = dataObject.neuron,
                           var = "seurat_clusters",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
ditto_umap
ditto_umap <- dittoDimPlot(object = dataObject.neuron,
                           var = "sample",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
ditto_umap

ditto_umap <- dittoDimPlot(object = dataObject.neuron,
                           var = "sample",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
ditto_umap

FeaturePlot(dataObject.neuron,
            reduction = "umap",
            #split.by = "sample",
            features = c("RBFOX3"))

FeaturePlot(dataObject.neuron,
            reduction = "umap",
            #split.by = "sample",
            features = c("GAD1"))

markers.to.plot <-
  c(
    "AIF1",
    "RBFOX3", 
    "SNAP25",
    "SYT1",
    "GAD1",
    "GAD2"
  )


dot_ind <- DotPlot(dataObject.neuron,
                   features = markers.to.plot, 
                   # split.by = "strain", 
                   cluster.idents = TRUE,
                   dot.scale = 8) + RotatedAxis()
dot_ind

## ----microglia-------------------------------------------------------------------------------------------------------
dataObject.microglia <- SCTransform(object_list$microglia, verbose = FALSE)

# run PCA on the merged object
dataObject.microglia<- RunPCA(object = dataObject.microglia)
Idents(dataObject.microglia) <- "sample"

# Determine the K-nearest neighbor graph
dataObject.microglia <- FindNeighbors(object = dataObject.microglia,
                                        assay = "SCT",
                                        reduction = "pca",
                                        dims = 1:15)
# Run UMAP
dataObject.microglia <- RunUMAP(dataObject.microglia,
                                  dims = 1:15,
                                  reduction = "pca",
                                  n.components = 3)

# Determine the clusters for various resolutions
dataObject.microglia <- FindClusters(object = dataObject.microglia,
                                       algorithm = 1, # 1= Louvain
                                       resolution = 0.7)

Idents(dataObject.microglia) <- dataObject.microglia$SCT_snn_res.0.7
dataObject.microglia$seurat_clusters <- dataObject.microglia$SCT_snn_res.0.7
ditto_umap <- dittoDimPlot(object = dataObject.microglia,
                           var = "group",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
ditto_umap
dataObject.neuron <- RunAzimuth(dataObject.neuron, reference = "humancortexref")
# currently, the object has two layers in the RNA assay: counts, and data
# DimPlot(dataObject, group.by = "predicted.subclass", label = TRUE, label.size = 3) + NoLegend()
dittoDimPlot(object = dataObject.neuron,
             var="predicted.subclass",
             reduction.use = "umap",
             do.label = FALSE,
             labels.highlight = FALSE)
## ----endothelial-------------------------------------------------------------------------------------------------------
# transform
dataObject.endothelial <- SCTransform(object_list$endothelial, verbose = FALSE)

# run PCA on the merged object
dataObject.endothelial <- RunPCA(object = dataObject.endothelial)
Idents(dataObject.endothelial) <- "sample"

# Determine the K-nearest neighbor graph
dataObject.endothelial <- FindNeighbors(object = dataObject.endothelial,
                                      assay = "SCT",
                                      reduction = "pca",
                                      dims = 1:15)
# Run UMAP
dataObject.endothelial <- RunUMAP(dataObject.endothelial,
                                dims = 1:15,
                                reduction = "pca",
                                n.components = 3)

# Determine the clusters for various resolutions
dataObject.endothelial <- FindClusters(object = dataObject.endothelial,
                                     algorithm = 1, # 1= Louvain
                                     resolution = 0.7)

Idents(dataObject.endothelial) <- dataObject.endothelial$SCT_snn_res.0.7
dataObject.endothelial$seurat_clusters <- dataObject.endothelial$SCT_snn_res.0.7
ditto_umap <- dittoDimPlot(object = dataObject.endothelial,
                           var = "sample",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
ditto_umap
## ----DE-------------------------------------------------------------------------------------------------------

gene_info <- read.delim("../rObjects/v7_ensembl_protein_coding_genes.txt")
protein_coding_genes <- subset(gene_info, gene_biotype == "protein_coding")

info <- (dataObject@meta.data)

output_dir <- "../results/DEGs/DESeq2_pseudobulk_exp_filter"
# List of cell types to analyze
# Pseudo-bulk the counts based on sample, group, and cell type
dataObject.pseudo <- AggregateExpression(
  dataObject.microglia, 
  assays = "RNA", # DESeq works with raw counts
  features = protein_coding_genes$gene_name,
  return.seurat = TRUE, 
  group.by = c("sample", "group", "individual_clusters")
)

# Function to calculate pct.1 and pct.2
calculate_percentages <- function(counts, group_labels) {
  percent_expressing <- function(group_counts) {
    apply(group_counts, 1, function(x) sum(x > 0) / length(x))
  }
  pct_list <- lapply(levels(group_labels), function(label) {
    percent_expressing(counts[, group_labels == label])
  })
  pct_df <- do.call(cbind, pct_list)
  colnames(pct_df) <- paste0("percent_", levels(group_labels))
  return(pct_df)
}


dataObject.pseudo$celltype.group <- paste(dataObject.pseudo$individual_clusters, dataObject.pseudo$group, sep = "_")
Idents(dataObject.pseudo) <- "celltype.group"
current_cells <- WhichCells(dataObject.pseudo, idents = paste("microglia", c("Ecoli", "Saline"), sep = "_"))
pseudo_bulk_counts <- GetAssayData(dataObject.pseudo, slot = "counts")[, current_cells]
gene_means <- rowMeans(pseudo_bulk_counts)
gene_means_summary <- summary(gene_means)
  # Filter genes with mean expression greater than 1st Qu. 
  filtered_genes <- names(gene_means[gene_means > gene_means_summary[[2]]])
  filtered_counts <- pseudo_bulk_counts[filtered_genes, ]
  pseudo_bulk_meta <- dataObject.pseudo@meta.data[current_cells, ]
  
  
  dds <- DESeqDataSetFromMatrix(
    countData = filtered_counts,
    colData = pseudo_bulk_meta,
    design = ~ group
  )
  
  dds <- DESeq(dds)
  res <- results(dds, contrast = c("group", "Saline", "Ecoli"))
  # Calculate pct.1 and pct.2
  pct_data <- calculate_percentages(filtered_counts, factor(pseudo_bulk_meta$group))
  # Combine results
  res$gene <- rownames(res)
  res <- as.data.frame(res)
  DEGs <- cbind(res[, c("gene", "baseMean", "log2FoldChange","lfcSE", "stat", "pvalue", "padj")])
  # Output
  group_output_path <- file.path(output_dir, paste0("microglia_Ecoli_vs_Saline_comparison_pseudobulk_v2.txt"))
  write.table(DEGs, file = group_output_path, sep = "\t", quote = FALSE, row.names = FALSE)
  rm(dds, res, DEGs) # clean up 
  
rm(dataObject.pseudo)

