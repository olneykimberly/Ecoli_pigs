## ----setup, include=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/Ecoli_pigs/snRNAseq/scripts/")
setwd("/tgen_labs/jfryer/kolney/Ecoli_pigs/snRNAseq/scripts/")

## ----echo=FALSE, message=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
source(here::here("/tgen_labs/jfryer/kolney/Ecoli_pigs/bulk_RNAseq/scripts/", "file_paths_and_colours.R"))
projectID <- "pigs_cellbender"
color.panel <- dittoColors()


# ----read_object--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
dataObject.singlets <- readRDS(paste0("../rObjects/",projectID,"_singlets.rds"))
dataObject.singlets
dataObject <- dataObject.singlets


## ----reprocess----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# transform
dataObject <- SCTransform(dataObject, verbose = FALSE)

# run PCA on the merged object
dataObject <- RunPCA(object = dataObject)
Idents(dataObject) <- "Sample_ID"
dataObject # inspect

# re-join layers
dataObject[["RNA"]] <- JoinLayers(dataObject[["RNA"]])
dataObject

# Determine the K-nearest neighbor graph
dataObject <- FindNeighbors(object = dataObject,
                                     assay = "SCT", # set as default after SCTransform
                                     reduction = "pca", # pca, harmony
                                     dims = 1:15)
# Run UMAP
dataObject <- RunUMAP(dataObject,
                               dims = 1:15,
                               reduction = "pca",
                               n.components = 3) # set to 3 to use with VR

# Determine the clusters for various resolutions
dataObject <- FindClusters(object = dataObject,
                                 algorithm = 1, # 1= Louvain
                                 resolution = seq(0.2,1,by=0.2))

saveRDS(dataObject, paste0("../rObjects/",projectID,"_unannotated_doublets_removed.rds"))
dataObject <- readRDS(paste0("../rObjects/",projectID,"_unannotated_doublets_removed.rds"))

## ----umap_noDoublets----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Idents(dataObject) <- dataObject$SCT_snn_res.0.4
dataObject$seurat_clusters <- dataObject$SCT_snn_res.0.4

## ----save_umap_noDoublets-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# save
path <- paste0("../results/UMAP/unannotated/",projectID,
               "_UMAP_unannotated_doublets_removed_SCTres0.2")
umap0.2 <- DimPlot(dataObject,
        group.by = "SCT_snn_res.0.2",
        label = TRUE)
umap0.2
saveToPDF(paste0(path, ".pdf"), width = 7, height = 6.6)

# save
path <- paste0("../results/UMAP/unannotated/",projectID,
               "_UMAP_unannotated_doublets_removed_SCTres0.4")
umap0.4 <- DimPlot(dataObject,
                   group.by = "SCT_snn_res.0.4",
                   label = TRUE)
umap0.4
saveToPDF(paste0(path, ".pdf"), width = 7, height = 6.6)

## ----Harmony integration-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(45)
dataObject.integrated <- IntegrateLayers(
  object = dataObject, method = HarmonyIntegration,
  normalization.method = "SCT",
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)
# Determine the K-nearest neighbor graph
dataObject.integrated <- FindNeighbors(object = dataObject.integrated, 
                                       reduction = "harmony", # pca, harmony 
                                       dims = 1:30)

# Determine the clusters for various resolutions
dataObject.integrated <- FindClusters(object = dataObject.integrated, resolution = 0.2)
dataObject.integrated <- RunUMAP(dataObject.integrated, reduction = "harmony", dims = 1:30)
saveRDS(dataObject.integrated, paste0("../rObjects/",projectID,"_unannotated_doublets_removed_harmony_int.rds"))
dataObject.integrated <- readRDS(paste0("../rObjects/",projectID,"_unannotated_doublets_removed_harmony_int.rds"))


p1 <- DimPlot(
  dataObject.integrated,
  reduction = "harmony",
  group.by = c("Sample_ID"),
  combine = FALSE, label.size = 2
)
p1
path <- paste0("../results/pca/",projectID,
               "_harmony_integration")
p1
saveToPDF(paste0(path, ".pdf"), width = 7, height = 6.6)

# umap - sample 
path <- paste0("../results/UMAP/unannotated/",projectID,
               "_UMAP_unannotated_doublets_removed_harmony_intergration_label_sample")
umap_harm_sample <- DimPlot(dataObject.integrated,
                            group.by = "Sample_ID")
umap_harm_sample
saveToPDF(paste0(path, ".pdf"), width = 7, height = 6.6)

# umap - group 
path <- paste0("../results/UMAP/unannotated/",projectID,
               "_UMAP_unannotated_doublets_removed_harmony_intergration_label_group")
umap_harm_group <- DimPlot(dataObject.integrated,
                            group.by = "group")
umap_harm_group
saveToPDF(paste0(path, ".pdf"), width = 7, height = 6.6)

## ----nuclei_per_cluster_noDoublets--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Idents(dataObject.integrated) <- dataObject.integrated$seurat_clusters
sample_ncells <- FetchData(dataObject.integrated, 
                     vars = c("ident", "Sample_ID")) %>%
  dplyr::count(ident,Sample_ID) %>%
  tidyr::spread(ident, n)
write.table(sample_ncells, 
            paste0("../results/nuclei_count/",
                   projectID, 
                   "_nuclei_per_cluster_doublets_removed.txt"),
            quote = FALSE, sep = "\t")
sample_ncells


## ----dot_individual-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
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

dot_ind <- DotPlot(dataObject.integrated,
                   features = markers.to.plot, 
                   cluster.idents = TRUE,
                   dot.scale = 8) + RotatedAxis()
dot_ind


## ----save_dot_individual, echo=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pdf(
  paste0(
    "../results/dot_plot/",
    projectID,
    "_clusters_DotPlot_doublets_removed_harmony_integration.pdf"
  ),
  width = 14,
  height = 10
)
dot_ind
dev.off()

## -----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


