## ----setup, include=FALSE--------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/Ecoli_pigs/snRNAseq/scripts/")
setwd("/tgen_labs/jfryer/kolney/Ecoli_pigs/snRNAseq/scripts/")

## ----echo=FALSE, message=FALSE---------------------------------------------------------------------------------------
source(here::here("/tgen_labs/jfryer/kolney/Ecoli_pigs/bulk_RNAseq/scripts/", "file_paths_and_colours.R"))
projectID <- "pigs_cellbender"
color.panel <- dittoColors()


## ----read_object-----------------------------------------------------------------------------------------------------
dataObject.dirty <- readRDS(paste0("../rObjects/",projectID,"_dirty.rds"))

# Visualize the number of cell counts per sample
data <- as.data.frame(table(dataObject.dirty$Sample_ID))
colnames(data) <- c("Sample_ID","frequency")

ncells <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) + 
  geom_col() +
  theme_classic() +
  geom_text(aes(label = frequency), 
            position=position_dodge(width=0.9), 
            #vjust=-0.25, 
            hjust = -.025,
            angle = 90) +
  #scale_fill_manual(values = sample_colors) + 
  scale_y_continuous(breaks = seq(0,60000, by = 5000), limits = c(0,60000)) +
  ggtitle("Removed dirty: nuclei per sample") +
  theme(legend.position =  "none") + 
  theme(axis.text.x = element_text(angle = 45, hjust=1))
path <- paste0("../results/nuclei_count/",projectID, 
               "_cells_per_sample_post_QC_and_doublets_removed_unannotated_after_recluster_dirty_nuclei")
ncells
saveToPDF(paste0(path, ".pdf"), width = 10, height = 6)

# mean cell count 
mean(data$frequency)

# median cell count 
median(data$frequency)

# Inspect 
dataObject.dirty
Layers(dataObject.dirty[["RNA"]]) # Inspect layers
dataObject.dirty[["RNA"]] <- JoinLayers(dataObject.dirty[["RNA"]])
# Inspect 
Layers(dataObject.dirty[["RNA"]])

# Re-normalizing data and finding clusters
dataObject.dirty <- SCTransform(dataObject.dirty, verbose = TRUE, conserve.memory = TRUE)
help(SCTransform)
all.genes <- rownames(dataObject.dirty)
# Inspect 
dataObject.dirty

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dataObject.dirty), 10)
top10

# run PCA on the merged object
dataObject.dirty <- RunPCA(object = dataObject.dirty)
Idents(dataObject.dirty) <- "Sample_ID"

# Determine the K-nearest neighbor graph
dataObject.dirty <- FindNeighbors(object = dataObject.dirty, 
                                  assay = "SCT", 
                                  reduction = "pca",
                                  dims = 1:15)
# Run UMAP
dataObject.dirty <- RunUMAP(dataObject.dirty,
                            dims = 1:15,
                            reduction = "pca",
                            n.components = 3) 

# Determine the clusters for various resolutions
dataObject.dirty <- FindClusters(object = dataObject.dirty,
                                 algorithm = 1, # 1= Louvain
                                 resolution = 0.9)

Idents(dataObject.dirty) <- dataObject.dirty$SCT_snn_res.0.9
dataObject.dirty$seurat_clusters <- dataObject.dirty$SCT_snn_res.0.9

# Save 
saveRDS(dataObject.dirty, paste0("../rObjects/",projectID,"_merged_reclusters_dirty_SCTransform.rds"))


ditto_umap <- dittoDimPlot(object = dataObject.dirty,
                           var = "seurat_clusters",
                           dim.1 = 1, 
                           dim.2 = 2,
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
ditto_umap
pdf(
  paste0(
    "../results/UMAP/unannotated/",
    projectID,
    "_individual_clusters_dirty_UMAP.pdf"
  ),
  width = 9,
  height = 7
)
ditto_umap
dev.off()

ditto_umap_anno <- dittoDimPlot(object = dataObject.dirty,
                                var = "individual_clusters",
                                dim.1 = 1, 
                                dim.2 = 2,
                                reduction.use = "umap",
                                do.label = TRUE,
                                labels.highlight = TRUE)
ditto_umap_anno


pdf(
  paste0(
    "../results/UMAP/annotated/",
    projectID,
    "_individual_clusters_dirty_UMAP.pdf"
  ),
  width = 9,
  height = 7
)
ditto_umap_anno
dev.off()

markers.to.plot <-
  c(
    "CLU", 
    "GFAP", 
    "AQP4", 
    "GJA1",
    "CLDN5",
    "ADGRF5",
    "FLT1",
    "COL1A1",
    "COL1A2",
    "DCN",
    "HEXB",
    "C1QA",
    "C1QB",
    "C1QC",
    "TMEM119",
    "ITGAM",
    "TYROBP",
    "P2RY12",
    "AIF1",
    "RBFOX1",
    "RBFOX3", 
    "SNAP25",
    "SYT1",
    "GAD1",
    "GAD2",
    "PLP1",
    "MBP", 
    "MOG", 
    "OLIG1",
    "PDGFRA",
    "VCAN",
    "TNR",
    "ACTA2",
    "VTN"
  )

dot_ind_celltype <- DotPlot(dataObject.dirty,
                            features = markers.to.plot, 
                            cluster.idents = FALSE,
                            dot.scale = 8) + RotatedAxis()
dot_ind_celltype
pdf(
  paste0(
    "../results/dot_plot/",
    projectID,
    "_individual_clusters_dirty.pdf"
  ),width = 10, height = 4)
dot_ind_celltype
dev.off()


