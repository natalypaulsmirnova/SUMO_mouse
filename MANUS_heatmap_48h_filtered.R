library(Seurat)
library(ggplot2)
library(viridis)
library(pheatmap)
library(dplyr)

####
# Heatmap of selected genes at 48h
# Data: 601_48h_filtered.rds
# Genes on X-axis, Conditions on Y-axis
# Color scheme: viridis
####

# Load the data
obj <- readRDS("601_manual_labeling_marioni_filtering_time_based/601_48h_filtered.rds")

# Set default assay to RNA
DefaultAssay(obj) <- "RNA"

# Merge layers in the RNA assay (has 20 layers that need to be merged)
cat("Merging layers in RNA assay...\n")
obj <- JoinLayers(obj, assay = "RNA")
cat("Layers merged successfully.\n")

# Define genes of interest in the exact order desired for the heatmap
genes <- c("Otx2", "Utf1", "Rasgrp2", "Dppa5a", "Fgf5", "Dnmt3a", "Sox2", 
           "Nodal", "Nanog", "Epha1", "Eomes", "Gsc", "Wnt3", "Mixl1", 
           "Wnt3a", "Nkx1-2", "T")

# Check which genes are available in the RNA assay
# Handle potential naming differences (e.g., Nkx1-2 might be Nkx1.2 in R)
available_genes <- rownames(obj[["RNA"]])
genes_found <- c()
genes_not_found <- c()

for (gene in genes) {
  # Try exact match first
  if (gene %in% available_genes) {
    genes_found <- c(genes_found, gene)
  } else {
    # Try case-insensitive match
    gene_match <- grep(paste0("^", gene, "$"), available_genes, ignore.case = TRUE, value = TRUE)
    if (length(gene_match) > 0) {
      genes_found <- c(genes_found, gene_match[1])
      cat("Found", gene, "as", gene_match[1], "\n")
    } else {
      # Try with dot instead of dash (Nkx1-2 -> Nkx1.2)
      gene_alt <- gsub("-", ".", gene)
      if (gene_alt %in% available_genes) {
        genes_found <- c(genes_found, gene_alt)
        cat("Found", gene, "as", gene_alt, "\n")
      } else {
        genes_not_found <- c(genes_not_found, gene)
      }
    }
  }
}

if (length(genes_not_found) > 0) {
  warning("The following genes were not found: ", paste(genes_not_found, collapse = ", "))
}

if (length(genes_found) == 0) {
  stop("No genes found in the object. Please check gene names.")
}

cat("Using", length(genes_found), "genes for heatmap\n")

# Get conditions (assuming they're stored in a "group" metadata column)
# Check available metadata columns
if (!"group" %in% colnames(obj@meta.data)) {
  # Try alternative column names
  if ("condition" %in% colnames(obj@meta.data)) {
    obj$group <- obj$condition
  } else if ("Group" %in% colnames(obj@meta.data)) {
    obj$group <- obj$Group
  } else {
    stop("Could not find condition/group column in metadata. Available columns: ", 
         paste(colnames(obj@meta.data), collapse = ", "))
  }
}

# Get unique conditions
conditions_found <- unique(obj$group)
cat("Conditions found in data:", paste(conditions_found, collapse = ", "), "\n")

# Define desired order: Ctrl, RACL, NACL, XAL
condition_order <- c("Ctrl", "RACL", "NACL", "XAL")

# Define label mapping: actual name -> display label
condition_labels <- c("Ctrl" = "Control",
                      "RACL" = "RACL",
                      "NACL" = "NACLB",
                      "XAL" = "XAL")

# Match found conditions to desired order (case-insensitive matching)
conditions_ordered <- c()
for (cond in condition_order) {
  # Try exact match first
  if (cond %in% conditions_found) {
    conditions_ordered <- c(conditions_ordered, cond)
  } else {
    # Try case-insensitive match
    match <- grep(paste0("^", cond, "$"), conditions_found, ignore.case = TRUE, value = TRUE)
    if (length(match) > 0) {
      conditions_ordered <- c(conditions_ordered, match[1])
      # Update label mapping to use actual condition name
      condition_labels[match[1]] <- condition_labels[cond]
    }
  }
}

# If we couldn't match all, use what we found
if (length(conditions_ordered) == 0) {
  warning("Could not match conditions to desired order. Using found conditions.")
  conditions_ordered <- conditions_found
  # Create default labels for unmatched conditions
  for (cond in conditions_ordered) {
    if (!cond %in% names(condition_labels)) {
      condition_labels[cond] <- cond
    }
  }
}

cat("Conditions ordered:", paste(conditions_ordered, collapse = ", "), "\n")

# Calculate average expression per condition
# Use AverageExpression from Seurat, explicitly using RNA assay with data slot
Idents(obj) <- obj$group

avg_expr <- AverageExpression(obj, features = genes_found, slot = "data", 
                               group.by = "group", assay = "RNA")

# Extract the expression matrix from RNA assay
# AverageExpression returns a list with assay name as key
if (is.list(avg_expr)) {
  if ("RNA" %in% names(avg_expr)) {
    expr_matrix <- avg_expr[["RNA"]]
  } else {
    # If RNA not found, try first element
    expr_matrix <- avg_expr[[1]]
    cat("Using first element from AverageExpression result\n")
  }
} else {
  # If it's not a list, it might be a matrix directly
  expr_matrix <- avg_expr
  cat("AverageExpression returned a matrix directly\n")
}

# Check if expr_matrix is valid
if (is.null(expr_matrix) || nrow(expr_matrix) == 0 || ncol(expr_matrix) == 0) {
  stop("Error: Expression matrix is empty or NULL. Check that genes and conditions are correct.")
}

cat("Expression matrix dimensions:", dim(expr_matrix), "\n")
cat("Matrix rownames (genes):", head(rownames(expr_matrix), 5), "...\n")
cat("Matrix colnames (conditions):", paste(colnames(expr_matrix), collapse = ", "), "\n")

# Transpose so conditions are rows and genes are columns
expr_matrix_t <- t(expr_matrix)

# Reorder genes to match exact original order: Otx2, Utf1, Rasgrp2, Dppa5a, Fgf5, Dnmt3a, Sox2, 
# Nodal, Nanog, Epha1, Eomes, Gsc, Wnt3, Mixl1, Wnt3a, Nkx1-2, T
# Map found genes back to original order, preserving the exact sequence
gene_order <- c()
genes_found_remaining <- genes_found  # Keep track of which genes haven't been added yet

for (gene in genes) {
  # Find matching gene in genes_found_remaining
  match_idx <- which(genes_found_remaining == gene | 
                     tolower(genes_found_remaining) == tolower(gene) |
                     gsub("\\.", "-", genes_found_remaining) == gene)
  if (length(match_idx) > 0) {
    matched_gene <- genes_found_remaining[match_idx[1]]
    gene_order <- c(gene_order, matched_gene)
    # Remove from remaining list to avoid duplicates
    genes_found_remaining <- genes_found_remaining[-match_idx[1]]
  }
}
# Add any remaining genes that weren't matched (shouldn't happen if all genes found)
if (length(genes_found_remaining) > 0) {
  gene_order <- c(gene_order, genes_found_remaining)
}

# Reorder columns to match gene order
expr_matrix_t <- expr_matrix_t[, gene_order, drop = FALSE]

# Reorder rows to match desired condition order and apply custom labels
# Only include conditions that exist in the matrix
conditions_in_matrix <- intersect(conditions_ordered, rownames(expr_matrix_t))
expr_matrix_t <- expr_matrix_t[conditions_in_matrix, , drop = FALSE]

# Apply custom labels to row names
row_labels <- condition_labels[rownames(expr_matrix_t)]
# Handle any missing labels (use original name if label not found)
row_labels[is.na(row_labels)] <- rownames(expr_matrix_t)[is.na(row_labels)]
rownames(expr_matrix_t) <- row_labels
expr_matrix_t <- as.matrix(expr_matrix_t)

# Create output directory
output_dir <- "MANUS_heatmap_48h_filtered"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create heatmap using pheatmap with viridis colors
# Generate viridis color palette
viridis_colors <- viridis(100)

# Create the heatmap with row labels on the left side
# Size: Width=8.5, Height=2 inches, compact with square cells 0.5x0.5 cm
# Convert 0.5 cm to inches: 0.5 cm ≈ 0.197 inches
cell_size_cm <- 0.5
cell_size_inches <- cell_size_cm / 2.54  # Convert cm to inches

png(file.path(output_dir, "MANUS_heatmap_48h_filtered.png"), 
    width = 8.5, height = 2, units = "in", res = 600)
pheatmap(expr_matrix_t,
         color = viridis_colors,
         cluster_rows = FALSE,  # Keep conditions in original order
         cluster_cols = FALSE,  # Keep genes in specified order
         scale = "none",  # Don't scale 
         main = NA,  # No title
         legend = FALSE,  # No legend
         fontsize = 6,
         fontsize_row = 7,
         fontsize_col = 6,
         angle_col = 45,
         display_numbers = FALSE,
         border_color = NA,
         show_rownames = TRUE,  # Ensure row names are shown on the left
         cellwidth = cell_size_inches * 72,  # Convert inches to points (72 points per inch)
         cellheight = cell_size_inches * 72,  # Square cells 0.5x0.5 cm
         gaps_row = 0,  # No gaps between rows for compact layout
         gaps_col = 0,  # No gaps between columns for compact layout
         treeheight_row = 0,  # No row dendrogram
         treeheight_col = 0)  # No column dendrogram
dev.off()

# Create another PNG version with legend
png(file.path(output_dir, "MANUS_heatmap_48h_filtered_with_legend.png"), 
    width = 9.5, height = 2, units = "in", res = 600)  # Slightly wider to accommodate legend
pheatmap(expr_matrix_t,
         color = viridis_colors,
         cluster_rows = FALSE,  # Keep conditions in original order
         cluster_cols = FALSE,  # Keep genes in specified order
         scale = "none",  # Don't scale - using scale.data which is already scaled
         main = NA,  # No title
         legend = TRUE,  # Include legend
         fontsize = 6,
         fontsize_row = 7,
         fontsize_col = 6,
         angle_col = 45,
         display_numbers = FALSE,
         border_color = NA,
         show_rownames = TRUE,  # Ensure row names are shown on the left
         cellwidth = cell_size_inches * 72,  # Convert inches to points (72 points per inch)
         cellheight = cell_size_inches * 72,  # Square cells 0.5x0.5 cm
         gaps_row = 0,  # No gaps between rows for compact layout
         gaps_col = 0,  # No gaps between columns for compact layout
         treeheight_row = 0,  # No row dendrogram
         treeheight_col = 0)  # No column dendrogram
dev.off()

# Also create PDF version with same compact dimensions
pdf(file.path(output_dir, "MANUS_heatmap_48h_filtered.pdf"), 
    width = 8.5, height = 2)
pheatmap(expr_matrix_t,
         color = viridis_colors,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         scale = "none",  # Don't scale 
         main = NA,  # No title
         legend = FALSE,  # No legend
         fontsize = 6,
         fontsize_row = 7,
         fontsize_col = 6,
         angle_col = 45,
         display_numbers = FALSE,
         border_color = NA,
         show_rownames = TRUE,
         cellwidth = cell_size_inches * 72,  # Square cells 0.5x0.5 cm
         cellheight = cell_size_inches * 72,
         gaps_row = 0,
         gaps_col = 0,
         treeheight_row = 0,
         treeheight_col = 0)
dev.off()

# Save the expression matrix as CSV
write.csv(expr_matrix_t, 
          file.path(output_dir, "MANUS_heatmap_48h_filtered_data.csv"),
          row.names = TRUE)

cat("Heatmap saved to", output_dir, "\n")
cat("Expression matrix saved to", file.path(output_dir, "MANUS_heatmap_48h_filtered_data.csv"), "\n")


