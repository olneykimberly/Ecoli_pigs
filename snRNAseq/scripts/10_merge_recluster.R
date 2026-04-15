## ----setup, include=FALSE--------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/Ecoli_pigs/snRNAseq/scripts/")
setwd("/tgen_labs/jfryer/kolney/Ecoli_pigs/snRNAseq/scripts/")

## ----echo=FALSE, message=FALSE---------------------------------------------------------------------------------------
source(here::here("/tgen_labs/jfryer/kolney/Ecoli_pigs/bulk_RNAseq/scripts/", "file_paths_and_colours.R"))
projectID <- "pigs_cellbender"
color.panel <- dittoColors()


## ----read_object-----------------------------------------------------------------------------------------------------
# read object
# dataObject.astrocyte <- readRDS(paste0("../rObjects/", projectID, "_astrocyte.rds"))
# dataObject.oligodendrocyte <- readRDS(paste0("../rObjects/",projectID,"_oligodendrocyte.rds"))
# dataObject.polydendrocyte<- readRDS(paste0("../rObjects/",projectID,"_polydendrocyte.rds"))
# dataObject.endothelial <- readRDS(paste0("../rObjects/",projectID,"_endothelial.rds"))
# dataObject.fibroblast <- readRDS(paste0("../rObjects/",projectID,"_fibroblast.rds"))
# dataObject.mural<- readRDS(paste0("../rObjects/",projectID,"_mural.rds"))
# dataObject.microglia <- readRDS(paste0("../rObjects/",projectID,"_microglia.rds"))
# dataObject.neuron <- readRDS(paste0("../rObjects/",projectID,"_neuron.rds"))
# dataObject.ependyma <- readRDS(paste0("../rObjects/",projectID,"_ependyma.rds"))
# 
# ## ----clean_objects-----------------------------------------------------------------------------------------------------
# dataObject.astrocyte.clean <- subset(dataObject.astrocyte, cells = WhichCells(dataObject.astrocyte, idents = setdiff(unique(dataObject.astrocyte$seurat_clusters), c(0, 4, 19, 15, 13, 20))))
# dataObject.astrocyte.dirty <- subset(dataObject.astrocyte, cells = WhichCells(dataObject.astrocyte, idents = c(0, 4, 19, 15, 13, 20)))
# rm(dataObject.astrocyte)
# 
# dataObject.oligodendrocyte.clean <- subset(dataObject.oligodendrocyte, cells = WhichCells(dataObject.oligodendrocyte, idents = setdiff(unique(dataObject.oligodendrocyte$seurat_clusters), c(7, 14))))
# dataObject.oligodendrocyte.dirty <- subset(dataObject.oligodendrocyte, cells = WhichCells(dataObject.oligodendrocyte, idents = c(7, 14)))
# rm(dataObject.oligodendrocyte)
# 
# dataObject.polydendrocyte.clean <- subset(dataObject.polydendrocyte, cells = WhichCells(dataObject.polydendrocyte, idents = setdiff(unique(dataObject.polydendrocyte$seurat_clusters), c(1, 2, 4, 6, 8, 10, 11, 12, 13, 15, 16, 18, 19, 20, 21))))
# dataObject.polydendrocyte.dirty <- subset(dataObject.polydendrocyte, cells = WhichCells(dataObject.polydendrocyte, idents = c(1, 2, 4, 6, 8, 10, 11, 12, 13, 15, 16, 18, 19, 20, 21)))
# rm(dataObject.polydendrocyte)
# 
# dataObject.endothelial.clean <- subset(dataObject.endothelial, cells = WhichCells(dataObject.endothelial, idents = setdiff(unique(dataObject.endothelial$seurat_clusters), c(16, 19, 8, 0, 20, 15, 11, 4, 21))))
# dataObject.endothelial.dirty <- subset(dataObject.endothelial, cells = WhichCells(dataObject.endothelial, idents = c(16, 19, 8, 0, 20, 15, 11, 4, 21)))
# rm(dataObject.endothelial)
# 
# dataObject.fibroblast.clean <- subset(dataObject.fibroblast, cells = WhichCells(dataObject.fibroblast, idents = setdiff(unique(dataObject.fibroblast$seurat_clusters), c(0, 2, 3, 5, 6, 7, 9, 11, 12, 13, 15))))
# dataObject.fibroblast.dirty <- subset(dataObject.fibroblast, cells = WhichCells(dataObject.fibroblast, idents = c(0, 2, 3, 5, 6, 7, 9, 11, 12, 13, 15)))
# rm(dataObject.fibroblast)
# 
# dataObject.mural.clean <- subset(dataObject.mural, cells = WhichCells(dataObject.mural, idents = setdiff(unique(dataObject.mural$seurat_clusters), c(0, 8, 4, 10, 13, 11, 6))))
# dataObject.mural.dirty <- subset(dataObject.mural, cells = WhichCells(dataObject.mural, idents = c(0, 8, 4, 10, 13, 11, 6)))
# rm(dataObject.mural)
# 
# dataObject.microglia.clean <- subset(dataObject.microglia, cells = WhichCells(dataObject.microglia, idents = setdiff(unique(dataObject.microglia$seurat_clusters), c(13, 17, 9, 16, 15))))
# dataObject.microglia.dirty <- subset(dataObject.microglia, cells = WhichCells(dataObject.microglia, idents =c(13, 17, 9, 16, 15)))
# rm(dataObject.microglia)
# 
# dataObject.neuron.clean <- subset(dataObject.neuron, cells = WhichCells(dataObject.neuron, idents = setdiff(unique(dataObject.neuron$seurat_clusters), c(0, 2, 3, 5, 6, 7, 8, 10, 11, 13, 17, 18, 19))))
# dataObject.neuron.dirty <- subset(dataObject.neuron, cells = WhichCells(dataObject.neuron, idents = c(0, 2, 3, 5, 6, 7, 8, 10, 11, 13, 17, 18, 19)))
# rm(dataObject.neuron)
# 
# dataObject.ependyma.clean <- subset(dataObject.ependyma, cells = WhichCells(dataObject.ependyma, idents = setdiff(unique(dataObject.ependyma$seurat_clusters), c(9, 7))))
# dataObject.ependyma.dirty <- subset(dataObject.ependyma, cells = WhichCells(dataObject.ependyma, idents = c(9, 7)))
# rm(dataObject.ependyma)
# 
# ## ----merge_objects-----------------------------------------------------------------------------------------------------
# dataObject.clean <- merge(x = dataObject.astrocyte.clean, y = c(dataObject.oligodendrocyte.clean, dataObject.polydendrocyte.clean, dataObject.endothelial.clean, dataObject.fibroblast.clean, dataObject.mural.clean, dataObject.microglia.clean, dataObject.neuron.clean, dataObject.ependyma.clean))
# saveRDS(dataObject.clean, paste0("../rObjects/",projectID,"_clean.rds"))
# # Inspect clean
# dataObject.clean
# 
# # Inspect dirty 
# dataObject.dirty <- merge(x = dataObject.astrocyte.dirty, y = c(dataObject.oligodendrocyte.dirty, dataObject.polydendrocyte.dirty, dataObject.endothelial.dirty, dataObject.fibroblast.dirty, dataObject.mural.dirty, dataObject.microglia.dirty, dataObject.neuron.dirty, dataObject.ependyma.dirty))
# saveRDS(dataObject.dirty, paste0("../rObjects/",projectID,"_dirty.rds"))
# 
# # Visualize the number of cell counts per sample
# data <- as.data.frame(table(dataObject.clean$Sample_ID))
# colnames(data) <- c("Sample_ID","frequency")
# 
# ncells <- ggplot(data, aes(x = Sample_ID, y = frequency, fill = Sample_ID)) + 
#   geom_col() +
#   theme_classic() +
#   geom_text(aes(label = frequency), 
#             position=position_dodge(width=0.9), 
#             #vjust=-0.25, 
#             hjust = -.025,
#             angle = 90) +
#   #scale_fill_manual(values = sample_colors) + 
#   scale_y_continuous(breaks = seq(0,60000, by = 5000), limits = c(0,60000)) +
#   ggtitle("Filtered: nuclei per sample") +
#   theme(legend.position =  "none") + 
#   theme(axis.text.x = element_text(angle = 45, hjust=1))
# path <- paste0("../results/nuclei_count/",projectID, 
#                "_cells_per_sample_post_QC_and_doublets_removed_unannotated_after_recluster")
# ncells
# saveToPDF(paste0(path, ".pdf"), width = 10, height = 6)
# 
# # mean cell count 
# mean(data$frequency)
# 
# # median cell count 
# median(data$frequency)
# 
# # Inspect 
# dataObject.clean
# Layers(dataObject.clean[["RNA"]]) # Inspect layers
# dataObject.clean[["RNA"]] <- JoinLayers(dataObject.clean[["RNA"]])
# # Inspect 
# Layers(dataObject.clean[["RNA"]])
# saveRDS(dataObject.clean, paste0("../rObjects/",projectID,"_merged_reclusters.rds"))
dataObject.clean <- readRDS(paste0("../rObjects/",projectID,"_merged_reclusters.rds"))

# Re-normalizing data and finding clusters
dataObject.clean <- SCTransform(dataObject.clean, verbose = FALSE, conserve.memory = TRUE)
all.genes <- rownames(dataObject.clean)
# Inspect 
dataObject.clean

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(dataObject.clean), 10)

# run PCA on the merged object
dataObject.clean <- RunPCA(object = dataObject.clean)
Idents(dataObject.clean) <- "Sample_ID"

# Determine the K-nearest neighbor graph
dataObject.clean <- FindNeighbors(object = dataObject.clean, 
                                  assay = "SCT", 
                                  reduction = "pca",
                                  dims = 1:15)
# Run UMAP
dataObject.clean <- RunUMAP(dataObject.clean,
                            dims = 1:15,
                            reduction = "pca",
                            n.components = 3) 

# Determine the clusters for various resolutions
dataObject.clean <- FindClusters(object = dataObject.clean,
                                 algorithm = 1, # 1= Louvain
                                 resolution = 0.9)

Idents(dataObject.clean) <- dataObject.clean$SCT_snn_res.0.9
dataObject.clean$seurat_clusters <- dataObject.clean$SCT_snn_res.0.9

# Save 
saveRDS(dataObject.clean, paste0("../rObjects/",projectID,"_merged_reclusters_clean_SCTransform.rds"))


ditto_umap <- dittoDimPlot(object = dataObject.clean,
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
    "_individual_clusters_clean_UMAP.pdf"
  ),
  width = 9,
  height = 7
)
ditto_umap
dev.off()

ditto_umap_anno <- dittoDimPlot(object = dataObject.clean,
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
    "_individual_clusters_clean_UMAP.pdf"
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


dot_ind_celltype <- DotPlot(dataObject.clean,
                            features = markers.to.plot, 
                            cluster.idents = FALSE,
                            dot.scale = 8) + RotatedAxis()
dot_ind_celltype
pdf(
  paste0(
    "../results/dot_plot/",
    projectID,
    "_individual_clusters_clean.pdf"
  ),width = 10, height = 4)
dot_ind_celltype
dev.off()


