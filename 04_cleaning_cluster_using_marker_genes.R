library(Seurat)
library(ggplot2)
library(tidyverse)
library(SeuratObject)


plot_top10_markers <- function(seurat_obj) {
  
  # Step 1: Join layers (safe for Seurat v5)
  seurat_obj <- SeuratObject::JoinLayers(seurat_obj)
  
  # Step 2: Find markers
  markers <- FindAllMarkers(seurat_obj,
                            only.pos = TRUE,
                            min.pct = 0.25,
                            logfc.threshold = 0.25)
  
  # Step 3: Get top 10 markers per cluster
  top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
  
  # Step 4: Unique feature list
  unique_features <- unique(top10$gene)
  
  # Step 5: Plot
  p <- DotPlot(seurat_obj, features = unique_features, group.by = "seurat_clusters") +
    coord_flip() +
    theme(axis.text.y = element_text(size = 10)) +
    ggtitle("Top 10 Markers per Cluster")
  
  return(p)
}

pdf("04_top10_markers.pdf", height = 12, width = 12)
sc_ctrl_regr_cc <- readRDS(file="04_Ctrl/04_Ctrl_sc_mouse_filtered_merged_seurat_cc.rds")
plot_top10_markers(sc_ctrl_regr_cc)
rm(sc_ctrl_regr_cc)
gc()

sc_racl_regr_cc <- readRDS(file="04_RACL/04_RACL_sc_mouse_filtered_merged_seurat_cc.rds")
plot_top10_markers(sc_racl_regr_cc)
rm(sc_racl_regr_cc)
gc()


sc_nacl_regr_cc <- readRDS(file="04_NACL/04_NACL_sc_mouse_filtered_merged_seurat_cc.rds")
plot_top10_markers(sc_nacl_regr_cc)
rm(sc_nacl_regr_cc)
gc()

sc_xal_regr_cc <- readRDS(file="04_XAL/04_XAL_sc_mouse_filtered_merged_seurat_cc.rds")
plot_top10_markers(sc_xal_regr_cc)
rm(sc_xal_regr_cc)
gc()

dev.off()

# Decided to remove cluster 4 and 5 from Ctrl 
sc_ctrl_regr_cc <- readRDS(file="04_Ctrl/04_Ctrl_sc_mouse_filtered_merged_seurat_cc.rds")

clusters_to_remove <- c(4, 5)
cells_to_remove <- WhichCells(sc_ctrl_regr_cc, idents = c(4,5))
sc_ctrl_regr_cc <- subset(sc_ctrl_regr_cc, 
                          cells = setdiff(Cells(sc_ctrl_regr_cc), cells_to_remove))
sc_ctrl_regr_cc <- SeuratObject::JoinLayers(sc_ctrl_regr_cc)
sc_ctrl_regr_cc <- CellCycleScoring(
  sc_ctrl_regr_cc,
  s.features = cc.genes$s.genes,
  g2m.features = cc.genes$g2m.genes,
  set.ident = FALSE
)

sc_ctrl_regr_cc <- RunPCA(sc_ctrl_regr_cc, verbose = TRUE)
sc_ctrl_regr_cc <- FindNeighbors(sc_ctrl_regr_cc, dims = 1:30)
sc_ctrl_regr_cc <- FindClusters(sc_ctrl_regr_cc, resolution = 0.5)
sc_ctrl_regr_cc <- RunUMAP(sc_ctrl_regr_cc, dims = 1:30)
