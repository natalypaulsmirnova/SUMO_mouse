library(Seurat)
library(CellChat)
library(patchwork)
library(tidyverse)
library(future)

obj <- readRDS("601_manual_labeling_marioni_filtering/601_nacl_filtered.rds")                     
celltype_col <- "manual_predicted.id" 
assay_to_use <- "RNA"
min_cells_per_group <- 10                       # <- drop extremely rare groups to stabilize probabilities
n_workers <- max(1, parallel::detectCores() - 6)

obj <- JoinLayers(obj)
#data.input <- seurat_object[["RNA"]]$data # normalized data matrix

# Put cell-type labels in Idents if needed
if (!is.null(celltype_col)) {
  stopifnot(celltype_col %in% colnames(obj@meta.data))
  Idents(obj) <- obj[[celltype_col, drop = TRUE]]
}
table(Idents(obj)) 

keep_ids <- names(which(table(Idents(obj)) >= min_cells_per_group))
obj <- subset(obj, idents = keep_ids)
Idents(obj) <- droplevels(Idents(obj))

DefaultAssay(obj) <- assay_to_use
expr <- GetAssayData(obj, assay = assay_to_use, slot = "data")  # log-normalized (or SCT corrected)
meta <- data.frame(labels = Idents(obj), row.names = colnames(obj))
unique(meta$labels)[1:5]

cellchat <- createCellChat(object = expr, meta = meta, group.by = "labels")

cellchat@DB <- CellChatDB

cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

plan("default")

cellchat <- computeCommunProb(cellchat, type = "triMean")
cellchat <- filterCommunication(cellchat, min.cells = min_cells_per_group)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellChat@idents))
group_size <- as.numeric(table(cellchat@idents))

netVisual_circle(cellchat@net$count, vertex.weight = group_size,
                 weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = group_size,
                 weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction strength")

cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming") 
netAnalysis_signalingRole_network(cellchat, signaling = pathway_show)

features <- c("Wnt3", "Bmp4", "Fgf8", "Notch1", "Dll1", "Jag1")  # edit for your system
DotPlot(obj, features = features, group.by = if (is.null(celltype_col)) NULL else celltype_col) + RotatedAxis()

#####################################################################

############################################################
objs <- SplitObject(obj, split.by = "time")

cellchat_list <- lapply(objs, function(o) {
  Idents(o) <- if (is.null(celltype_col)) Idents(o) else o[[celltype_col, drop = TRUE]]
  o <- subset(o, idents = names(which(table(Idents(o)) >= min_cells_per_group)))
  expr <- GetAssayData(o, assay = assay_to_use, slot = "data")
  meta <- data.frame(labels = Idents(o), row.names = colnames(o))
  cc <- createCellChat(expr, meta = meta, group.by = "labels")
  cc@DB <- CellChatDB.mouse
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc)
  cc <- filterCommunication(cc, min.cells = min_cells_per_group)
  cc <- computeCommunProbPathway(cc)
  aggregateNet(cc)
})

# Merge and compare
cellchat_merged <- mergeCellChat(cellchat_list, add.names = names(objs))
# Global comparison of interaction strength/counts between conditions
compareInteractions(cellchat_merged, show.legend = TRUE, group = c(3,4))  # e.g., 48h vs 72h
netVisual_diffInteraction(cellchat_merged)
rankNet(cellchat_merged, mode = "comparison")
############################################################


pathway_show <- c("ncWNT", "BMP", "FGF", "NOTCH")  # adjust as present
for (p in pathway_show) {
  if (p %in% names(cellchat@netP$pathways)) {
    print(netVisual_aggregate(cellchat, signaling = p, layout = "circle"))
    print(netVisual_bubble(cellchat, sources.use = NULL, targets.use = NULL, signaling = p))
  }
}

cellChat <- identifyCommunicationPatterns(cellChat, pattern = "outgoing", k = 5)

netVisual_aggregate(cellChat, "ncWNT", vertex.receiver = 20:29, layout="hierarchy")

netVisual_hierarchy2(cellChat, "ncWNT")
