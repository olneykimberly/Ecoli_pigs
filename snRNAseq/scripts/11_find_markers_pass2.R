knitr::opts_knit$set(root.dir = "/tgen_labs/jfryer/kolney/Ecoli_pigs/snRNAseq/scripts/")
setwd("/tgen_labs/jfryer/kolney/Ecoli_pigs/snRNAseq/scripts/")

source(here::here("/tgen_labs/jfryer/kolney/Ecoli_pigs/bulk_RNAseq/scripts/", "file_paths_and_colours.R"))
projectID <- "pigs_cellbender"
color.panel <- dittoColors()

# read object
dataObject <- readRDS(paste0("../rObjects/",projectID,"_after_recluster_harm_int.rds"))
# inspect
dataObject
DefaultAssay(dataObject) 
Layers(dataObject[["RNA"]])

# set the default assay to RNA
DefaultAssay(dataObject) <- "RNA"
dataObject$seurat_clusters <- dataObject$SCT_snn_res.0.9
Idents(dataObject) <- "seurat_clusters"
dataObject <- NormalizeData(dataObject)
# inspect
dataObject


## Markers per cluster
markers <- SeuratWrappers::RunPrestoAll(
  object = dataObject,
  assay = "RNA",
  slot = "counts",
  only.pos = FALSE
)
write.table(markers, 
            paste0("../results/markers/", projectID, "_after_recluster_harm_int_markers.tsv"),
            quote = FALSE,
            row.names = FALSE)
saveRDS(markers, paste0("../rObjects/", projectID, "_after_recluster_harm_int_markers.rds"))
markers <- readRDS(paste0("../rObjects/", projectID,"_after_recluster_harm_int_markers.rds"))

# rearrange to order by cluster & filter to only include log2FC > 1 & FDR < 0.05
all.markers.strict <- markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1 & p_val_adj < 0.05)

saveRDS(all.markers.strict, paste0("../rObjects/", projectID,"_after_recluster_harm_int_markers_log2FC1_q0.01.rds"))
all.markers.strict <- readRDS(paste0("../rObjects/", projectID,"_after_recluster_harm_int_markers_log2FC1_q0.01.rds"))