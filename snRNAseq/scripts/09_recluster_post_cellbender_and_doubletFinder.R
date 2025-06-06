## ----setup, include=FALSE--------------------------------------------------------------------------------------------
knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/Ecoli_pigs/snRNAseq/scripts/")
setwd("/tgen_labs/jfryer/kolney/Ecoli_pigs/snRNAseq/scripts/")

## ----echo=FALSE, message=FALSE---------------------------------------------------------------------------------------
source(here::here("/tgen_labs/jfryer/kolney/Ecoli_pigs/bulk_RNAseq/scripts/", "file_paths_and_colours.R"))
projectID <- "pigs_cellbender"
color.panel <- dittoColors()
set.seed(28)

## ----read_object-----------------------------------------------------------------------------------------------------
# read object
dataObject.annotated <- readRDS(paste0("../rObjects/",projectID,"_annotated_before_recluster.rds"))
# inspect
dataObject.annotated
# Split by cell type 
object_list <- SplitObject(dataObject.annotated, split.by = "individual_clusters")

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

# cell types
# mural, endothelial, ependyma, oligodendrocyte, fibroblast, polydendrocytes, astrocyte, neuron, microglia  

## ----astrocyte-------------------------------------------------------------------------------------------------------
# transform
dataObject.astrocyte <- SCTransform(object_list$astrocyte, verbose = FALSE)

# run PCA on the merged object
dataObject.astrocyte <- RunPCA(object = dataObject.astrocyte)
Idents(dataObject.astrocyte) <- "sample"

# Determine the K-nearest neighbor graph
dataObject.astrocyte <- FindNeighbors(object = dataObject.astrocyte,
                                      assay = "SCT",
                                      reduction = "pca",
                                      dims = 1:15)
# Run UMAP
dataObject.astrocyte <- RunUMAP(dataObject.astrocyte,
                                dims = 1:15,
                                reduction = "pca",
                                n.components = 3)

# Determine the clusters for various resolutions
dataObject.astrocyte <- FindClusters(object = dataObject.astrocyte,
                                     algorithm = 1, # 1= Louvain
                                     resolution = 0.7)

Idents(dataObject.astrocyte) <- dataObject.astrocyte$SCT_snn_res.0.7
dataObject.astrocyte$seurat_clusters <- dataObject.astrocyte$SCT_snn_res.0.7
ditto_umap <- dittoDimPlot(object = dataObject.astrocyte,
                           var = "seurat_clusters",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
ditto_umap
pdf(
  paste0("../results/UMAP/reclusters/",
         projectID, "_recluster_astrocyte.pdf"
  ),
  width = 7,
  height = 6
)
ditto_umap
dev.off()


dot_ind <- DotPlot(dataObject.astrocyte,
                   features = markers.to.plot,
                   cluster.idents = FALSE,
                   dot.scale = 8) + RotatedAxis()
pdf(
  paste0("../results/dot_plot/reclusters/",
         projectID, "_recluster_astrocyte.pdf"
  ),
  width = 14,
  height = 7
)
dot_ind
dev.off()

saveRDS(dataObject.astrocyte, paste0("../rObjects/",projectID,"_astrocyte.rds"))


## ----oligodendrocyte-------------------------------------------------------------------------------------------------
# transform
dataObject.oligodendrocyte <- SCTransform(object_list$oligodendrocyte, verbose = FALSE)

# run PCA on the merged object
dataObject.oligodendrocyte <- RunPCA(object = dataObject.oligodendrocyte)
Idents(dataObject.oligodendrocyte) <- "sample"

# Determine the K-nearest neighbor graph
dataObject.oligodendrocyte <- FindNeighbors(object = dataObject.oligodendrocyte,
                                            assay = "SCT",
                                            reduction = "pca",
                                            dims = 1:15)
# Run UMAP
dataObject.oligodendrocyte <- RunUMAP(dataObject.oligodendrocyte,
                                      dims = 1:15,
                                      reduction = "pca",
                                      n.components = 3)

# Determine the clusters for various resolutions
dataObject.oligodendrocyte <- FindClusters(object = dataObject.oligodendrocyte,
                                           algorithm = 1, # 1= Louvain
                                           resolution = 0.7)

Idents(dataObject.oligodendrocyte) <- dataObject.oligodendrocyte$SCT_snn_res.0.7
dataObject.oligodendrocyte$seurat_clusters <- dataObject.oligodendrocyte$SCT_snn_res.0.7
ditto_umap <- dittoDimPlot(object = dataObject.oligodendrocyte,
                           var = "seurat_clusters",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
ditto_umap
pdf(
  paste0("../results/UMAP/reclusters/",
         projectID, "_recluster_oligodendrocyte.pdf"
  ),
  width = 7,
  height = 6
)
ditto_umap
dev.off()

dot_ind <- DotPlot(dataObject.oligodendrocyte,
                   features = markers.to.plot,

                   cluster.idents = FALSE,
                   dot.scale = 8) + RotatedAxis()
dot_ind
pdf(
  paste0("../results/dot_plot/reclusters/",
         projectID, "_recluster_oligodendrocyte.pdf"
  ),
  width = 14,
  height = 7
)
dot_ind
dev.off()

saveRDS(dataObject.oligodendrocyte, paste0("../rObjects/",projectID,"_oligodendrocyte.rds"))


## ----polydendrocyte--------------------------------------------------------------------------------------------------
# transform
dataObject.polydendrocyte <- SCTransform(object_list$polydendrocyte, verbose = FALSE)

# run PCA on the merged object
dataObject.polydendrocyte <- RunPCA(object = dataObject.polydendrocyte)
Idents(dataObject.polydendrocyte) <- "sample"

# Determine the K-nearest neighbor graph
dataObject.polydendrocyte <- FindNeighbors(object = dataObject.polydendrocyte,
                                           assay = "SCT",
                                           reduction = "pca",
                                           dims = 1:15)
# Run UMAP
dataObject.polydendrocyte <- RunUMAP(dataObject.polydendrocyte,
                                     dims = 1:15,
                                     reduction = "pca",
                                     n.components = 3)

# Determine the clusters for various resolutions
dataObject.polydendrocyte <- FindClusters(object = dataObject.polydendrocyte,
                                          algorithm = 1, # 1= Louvain
                                          resolution = 0.7)

Idents(dataObject.polydendrocyte) <- dataObject.polydendrocyte$SCT_snn_res.0.7
dataObject.polydendrocyte$seurat_clusters <- dataObject.polydendrocyte$SCT_snn_res.0.7
ditto_umap <- dittoDimPlot(object = dataObject.polydendrocyte,
                           var = "seurat_clusters",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
ditto_umap
pdf(
  paste0("../results/UMAP/reclusters/",
         projectID, "_recluster_polydendrocyte.pdf"
  ),
  width = 7,
  height = 6
)
ditto_umap
dev.off()


dot_ind <- DotPlot(dataObject.polydendrocyte,
                   features = markers.to.plot,

                   cluster.idents = FALSE,
                   dot.scale = 8) + RotatedAxis()
dot_ind
pdf(
  paste0("../results/dot_plot/reclusters/",
         projectID, "_recluster_polydendrocyte.pdf"
  ),
  width = 14,
  height = 7
)
dot_ind
dev.off()

saveRDS(dataObject.polydendrocyte, paste0("../rObjects/",projectID,"_polydendrocyte.rds"))

## ----endothelial-----------------------------------------------------------------------------------------------------
# transform
dataObject.endothelial <- SCTransform(object_list$endothelial, verbose = FALSE)
# run PCA on the merged object
dataObject.endothelial <- RunPCA(object = dataObject.endothelial)
dataObject.endothelial
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
                           var = "seurat_clusters",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
ditto_umap
pdf(
  paste0("../results/UMAP/reclusters/",
         projectID, "_recluster_endothelial.pdf"
  ),
  width = 7,
  height = 6
)
ditto_umap
dev.off()

dot_ind <- DotPlot(dataObject.endothelial,
                    features = markers.to.plot, 
                    cluster.idents = TRUE,
                    dot.scale = 8) + RotatedAxis()
pdf(
   paste0("../results/dot_plot/reclusters/",
          projectID, "_recluster_endothelial.pdf"
   ),
   width = 14,
   height = 7
 )
dot_ind
dev.off()

saveRDS(dataObject.endothelial, paste0("../rObjects/",projectID,"_endothelial.rds"))


## ----fibroblast------------------------------------------------------------------------------------------------------
# transform
dataObject.fibroblast <- SCTransform(object_list$fibroblast, verbose = FALSE)

# run PCA on the merged object
dataObject.fibroblast <- RunPCA(object = dataObject.fibroblast)
Idents(dataObject.fibroblast) <- "sample"

# Determine the K-nearest neighbor graph
dataObject.fibroblast <- FindNeighbors(object = dataObject.fibroblast, 
                                       assay = "SCT", 
                                       reduction = "pca",
                                       dims = 1:15)
# Run UMAP
dataObject.fibroblast <- RunUMAP(dataObject.fibroblast,
                                 dims = 1:15,
                                 reduction = "pca",
                                 n.components = 3) 

# Determine the clusters for various resolutions
dataObject.fibroblast <- FindClusters(object = dataObject.fibroblast,
                                      algorithm = 1, # 1= Louvain
                                      resolution = 0.7)

Idents(dataObject.fibroblast) <- dataObject.fibroblast$SCT_snn_res.0.7
dataObject.fibroblast$seurat_clusters <- dataObject.fibroblast$SCT_snn_res.0.7
ditto_umap <- dittoDimPlot(object = dataObject.fibroblast,
                           var = "seurat_clusters",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
ditto_umap
pdf(
  paste0("../results/UMAP/reclusters/",
         projectID, "_recluster_fibroblast.pdf"
  ),
  width = 7,
  height = 6
)
ditto_umap
dev.off()

dot_ind <- DotPlot(dataObject.fibroblast,
                   features = markers.to.plot, 
                   cluster.idents = FALSE,
                   dot.scale = 8) + RotatedAxis()
dot_ind
pdf(
  paste0("../results/dot_plot/reclusters/",
         projectID, "_recluster_fibroblast.pdf"
  ),
  width = 14,
  height = 7
)
dot_ind
dev.off()

saveRDS(dataObject.fibroblast, paste0("../rObjects/",projectID,"_fibroblast.rds"))

## ----mural-----------------------------------------------------------------------------------------------------------
# transform
dataObject.mural <- SCTransform(object_list$mural, verbose = FALSE)

# run PCA on the merged object
dataObject.mural <- RunPCA(object = dataObject.mural)
Idents(dataObject.mural) <- "sample"

# Determine the K-nearest neighbor graph
dataObject.mural <- FindNeighbors(object = dataObject.mural, 
                                  assay = "SCT", 
                                  reduction = "pca",
                                  dims = 1:15)
# Run UMAP
dataObject.mural <- RunUMAP(dataObject.mural,
                            dims = 1:15,
                            reduction = "pca",
                            n.components = 3) 

# Determine the clusters for various resolutions
dataObject.mural <- FindClusters(object = dataObject.mural,
                                 algorithm = 1, # 1= Louvain
                                 resolution = 0.7)

Idents(dataObject.mural) <- dataObject.mural$SCT_snn_res.0.7
dataObject.mural$seurat_clusters <- dataObject.mural$SCT_snn_res.0.7
ditto_umap <- dittoDimPlot(object = dataObject.mural,
                           var = "seurat_clusters",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
ditto_umap
pdf(
  paste0("../results/UMAP/reclusters/",
         projectID, "_recluster_mural.pdf"
  ),
  width = 7,
  height = 6
)
ditto_umap
dev.off()


dot_ind <- DotPlot(dataObject.mural,
                   features = markers.to.plot, 
                   cluster.idents = FALSE,
                   dot.scale = 8) + RotatedAxis()
dot_ind
pdf(
  paste0("../results/dot_plot/reclusters/",
         projectID, "_recluster_mural.pdf"
  ),
  width = 14,
  height = 7
)
dot_ind
dev.off()

saveRDS(dataObject.mural, paste0("../rObjects/",projectID,"_mural.rds"))

## ----microglia-------------------------------------------------------------------------------------------------------
# transform
dataObject.microglia <- SCTransform(object_list$microglia, verbose = FALSE)

# run PCA on the merged object
dataObject.microglia <- RunPCA(object = dataObject.microglia)
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
                           var = "seurat_clusters",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
ditto_umap
pdf(
  paste0("../results/UMAP/reclusters/",
         projectID, "_recluster_microglia.pdf"
  ),
  width = 7,
  height = 6
)
ditto_umap
dev.off()


dot_ind <- DotPlot(dataObject.microglia,
                   features = markers.to.plot, 
                   
                   cluster.idents = TRUE,
                   dot.scale = 8) + RotatedAxis()
dot_ind
pdf(
  paste0("../results/dot_plot/reclusters/",
         projectID, "_recluster_microglia.pdf"
  ),
  width = 14,
  height = 7
)
dot_ind
dev.off()

saveRDS(dataObject.microglia, paste0("../rObjects/",projectID,"_microglia.rds"))

## ----neuron----------------------------------------------------------------------------------------------------------
# transform
dataObject.neuron <- SCTransform(object_list$neuron, verbose = FALSE)

# run PCA on the merged object
dataObject.neuron <- RunPCA(object = dataObject.neuron)
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
pdf(
  paste0("../results/UMAP/reclusters/",
         projectID, "_recluster_neuron.pdf"
  ),
  width = 7,
  height = 6
)
ditto_umap
dev.off()


dot_ind <- DotPlot(dataObject.neuron,
                   features = markers.to.plot, 
                   cluster.idents = FALSE,
                   dot.scale = 8) + RotatedAxis()
dot_ind
pdf(
  paste0("../results/dot_plot/reclusters/",
         projectID, "_recluster_neuron.pdf"
  ),
  width = 14,
  height = 7
)
dot_ind
dev.off()

saveRDS(dataObject.neuron, paste0("../rObjects/",projectID,"_neuron.rds"))

## ----ependyma----------------------------------------------------------------------------------------------------------
# transform
dataObject.ependyma <- SCTransform(object_list$ependyma, verbose = FALSE)

# run PCA on the merged object
dataObject.ependyma <- RunPCA(object = dataObject.ependyma)
Idents(dataObject.ependyma) <- "sample"

# Determine the K-nearest neighbor graph
dataObject.ependyma <- FindNeighbors(object = dataObject.ependyma, 
                                     assay = "SCT", 
                                     reduction = "pca",
                                     dims = 1:15)
# Run UMAP
dataObject.ependyma <- RunUMAP(dataObject.ependyma,
                               dims = 1:15,
                               reduction = "pca",
                               n.components = 3) 

# Determine the clusters for various resolutions
dataObject.ependyma <- FindClusters(object = dataObject.ependyma,
                                    algorithm = 1, # 1= Louvain
                                    resolution = 0.7)

Idents(dataObject.ependyma) <- dataObject.ependyma$SCT_snn_res.0.7
dataObject.ependyma$seurat_clusters <- dataObject.ependyma$SCT_snn_res.0.7
ditto_umap <- dittoDimPlot(object = dataObject.ependyma,
                           var = "seurat_clusters",
                           reduction.use = "umap",
                           do.label = TRUE,
                           labels.highlight = TRUE)
ditto_umap
pdf(
  paste0("../results/UMAP/reclusters/",
         projectID, "_recluster_ependyma.pdf"
  ),
  width = 7,
  height = 6
)
ditto_umap
dev.off()


dot_ind <- DotPlot(dataObject.ependyma,
                   features = markers.to.plot, 
                   cluster.idents = FALSE,
                   dot.scale = 8) + RotatedAxis()
dot_ind
pdf(
  paste0("../results/dot_plot/reclusters/",
         projectID, "_recluster_ependyma.pdf"
  ),
  width = 14,
  height = 7
)
dot_ind
dev.off()

saveRDS(dataObject.ependyma, paste0("../rObjects/",projectID,"_ependyma.rds"))
