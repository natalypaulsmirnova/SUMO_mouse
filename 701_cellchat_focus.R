# =============================================================================
# CellChat Analysis Script - Focused on Specific Cell Types
# =============================================================================
# This script performs cell-cell communication analysis using CellChat
# on a subset of specific cell types from a single-cell RNA-seq dataset

# Load required libraries
library(Seurat)
library(CellChat)
library(dplyr)

# =============================================================================
# Data Loading and Configuration
# =============================================================================

# Load the pre-processed Seurat object (RACL condition with manual cell type labels)
obj <- readRDS("601_manual_labeling_marioni_filtering/601_racl_filtered.rds")                     

# Define parameters for the analysis
celltype_col <- "manual_marioni_label"  # Column name containing cell type annotations
assay_to_use <- "RNA"                   # Assay to use for expression data
min_cells_per_group <- 10               # Minimum cells per group to stabilize probabilities

# Define specific cell types to focus on for communication analysis
# These represent key developmental stages in early mouse embryogenesis
celltypes_to_keep <- c("Epiblast", "Primitive Streak", "Anterior Primitive Streak",
                       "Visceral endoderm", "ExE endoderm", "Parietal endoderm", "Caudal epiblast",
                       "Ectoderm")

# =============================================================================
# Data Preprocessing
# =============================================================================

# Join layers if the object has multiple layers (e.g., from SCTransform)
obj <- JoinLayers(obj)

# Set cell type labels as the active identity in the Seurat object
if (!is.null(celltype_col)) {
  stopifnot(celltype_col %in% colnames(obj@meta.data))
  Idents(obj) <- obj[[celltype_col, drop = TRUE]]
}

# Ensure Idents are properly set
Idents(obj) <- Idents(obj)

# =============================================================================
# Cell Type Filtering and Quality Control
# =============================================================================

# Subset to only include the specified cell types of interest
obj_sub <- subset(obj, idents = celltypes_to_keep)

# Remove the original object to free up memory
rm(obj)

# Remove unused factor levels from the cell type identities
Idents(obj_sub) <- droplevels(Idents(obj_sub))

# Further filter to remove cell types with fewer than min_cells_per_group cells
# This helps stabilize probability calculations in CellChat
obj_sub <- subset(obj_sub, idents = names(which(table(Idents(obj_sub)) >= min_cells_per_group)))

# =============================================================================
# Prepare Data for CellChat Analysis
# =============================================================================

# Set the active assay
DefaultAssay(obj_sub) <- assay_to_use

# Extract normalized expression data
expr <- GetAssayData(obj_sub, assay = assay_to_use, slot = "data")

# Create metadata dataframe with cell type labels
meta <- data.frame(labels = Idents(obj_sub), row.names = colnames(obj_sub))

# =============================================================================
# CellChat Analysis Pipeline
# =============================================================================

# Create CellChat object with expression data and metadata
cellchat <- createCellChat(object = expr, meta = meta, group.by = "labels")

# Set the database to mouse (CellChatDB.mouse contains mouse ligand-receptor pairs)
cellchat@DB <- CellChatDB.mouse

# Subset the database to only include genes expressed in the dataset
cellchat <- subsetData(cellchat)

# Identify over-expressed genes in each cell group
cellchat <- identifyOverExpressedGenes(cellchat)

# Identify over-expressed ligand-receptor interactions
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute communication probabilities using triMean method
# triMean is robust to outliers and provides stable probability estimates
cellchat <- computeCommunProb(cellchat, type = "triMean")

# Filter out communications with too few cells
cellchat <- filterCommunication(cellchat, min.cells = min_cells_per_group)

# Compute communication probabilities at pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate the communication network
cellchat <- aggregateNet(cellchat)

# Compute network centrality measures for pathway analysis
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")

# =============================================================================
# Visualization and Results
# =============================================================================

# Calculate group sizes for visualization
group_size <- as.numeric(table(cellchat@idents))

# Create circular network plot showing interaction strength
# Vertex size represents the number of cells in each group
netVisual_circle(cellchat@net$weight, vertex.weight = group_size, title.name = "Interaction strength")

# =============================================================================
# Data Exploration
# =============================================================================

# Check the cell types present in the analysis
levels(cellchat@idents)

# Check the dimensions of the communication network
dim(cellchat@net$weight)

# Extract detailed communication information
pair_pathways <- subsetCommunication(cellchat)
