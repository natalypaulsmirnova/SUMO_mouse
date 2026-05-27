# Load required libraries
library(Seurat)
library(CellChat)
library(dplyr)
library(ggplot2)


# Function to check cell type availability and log missing types
check_cell_types <- function(obj, condition, celltypes_to_keep) {
  available_types <- levels(Idents(obj))
  missing_types <- setdiff(celltypes_to_keep, available_types)
  
  cat("  Available cell types in", condition, ":", length(available_types), "\n")
  cat("  Cell types:", paste(available_types, collapse = ", "), "\n")
  
  if (length(missing_types) > 0) {
    cat("  WARNING: Missing cell types in", condition, ":", paste(missing_types, collapse = ", "), "\n")
  }
  
  # Return only the cell types that are actually present
  present_types <- intersect(celltypes_to_keep, available_types)
  return(present_types)
}

conditions <- c("ctrl", "nacl", "racl", "xal")
min_cells_per_group <- 10

seurat_obj <- readRDS("601_manual_labeling_marioni_filtering/601_racl_filtered.rds")

seurat_obj <- JoinLayers(seurat_obj)

Idents(seurat_obj) <- "time"
seurat_obj <- subset(seurat_obj, subset = time == "120h")

metadata <- seurat_obj@meta.data

metadata <- metadata %>% mutate(cellchat_collapse = case_when(
  manual_marioni_label == "Parietal endoderm" ~ "ExEnd",
  manual_marioni_label == "Visceral endoderm" ~ "ExEnd",
  manual_marioni_label == "Venous endothelium" ~ "Endothelium",
  manual_marioni_label == "Embryo proper endothelium" ~ "Endothelium",
  manual_marioni_label == "Foregut" ~ "Gut",
  manual_marioni_label == "Gut tube" ~ "Gut",
  manual_marioni_label == "Midgut" ~ "Gut",
  manual_marioni_label == "Hindgut" ~ "Gut",
  manual_marioni_label == "Presomitic mesoderm" ~ "Somite",
  manual_marioni_label == "Somitic mesoderm" ~ "Somite",
  manual_marioni_label == "Posterior somitic mesoderm" ~ "Somite",
  manual_marioni_label == "Posterior somitic tissues" ~ "Somite",
  manual_marioni_label == "Anterior somitic tissues" ~ "Somite",
  manual_marioni_label == "Dermomyotome" ~ "Somite",
  manual_marioni_label == "Sclerotome" ~ "Somite",
  manual_marioni_label == "Endotome" ~ "Somite",
  manual_marioni_label == "Neural tube" ~ "Neural",
  manual_marioni_label == "Spinal cord progenitors" ~ "Neural",
  manual_marioni_label == "Dorsal spinal cord progenitors" ~ "Neural",
  manual_marioni_label == "Hindbrain neural progenitors" ~ "Neural",
  manual_marioni_label == "Hindbrain floor plate" ~ "Neural",
  manual_marioni_label == "Dorsal midbrain neurons" ~ "Neural",
  manual_marioni_label == "Dorsal hindbrain progenitors" ~ "Neural",
  manual_marioni_label == "Dorsal spinal cord progenitors" ~ "Neural",
  manual_marioni_label == "Ventral hindbrain progenitors" ~ "Neural",
  .default = manual_marioni_label
))

seurat_obj <- AddMetaData(seurat_obj, metadata = metadata)
Idents(seurat_obj) <- "cellchat_collapse"
seurat_obj <- subset(seurat_obj, 
                     idents = names(which(table(Idents(seurat_obj)) >= min_cells_per_group)))
Idents(seurat_obj) <- droplevels(Idents(seurat_obj))
DefaultAssay(seurat_obj) <- "RNA"

expr <- GetAssayData(seurat_obj, assay = "RNA", layer = "data")
meta <- data.frame(labels = Idents(seurat_obj), row.names = colnames(seurat_obj))

# Create CellChat object with expression data and metadata
cellchat <- createCellChat(object = expr, meta = meta, group.by = "labels")
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = min_cells_per_group)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

group_size <- as.numeric(table(cellchat@idents))
pair_pathways <- subsetCommunication(cellchat)




