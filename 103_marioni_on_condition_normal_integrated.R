library(SingleCellExperiment)
library(Seurat)
library(scater)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

library(ComplexHeatmap)
library(cluster)
library(circlize)

# Create output directory
output_dir <- "103_marioni_on_condition_normal_integrated"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Load reference dataset
ref_seurat <- readRDS("Mariono_mouse_atlas/ExtendedMouseAtlas/marion_seurat_v2.rds")

# Check reference dataset normalization
cat("Reference dataset info:\n")
cat("Assays available:", names(ref_seurat@assays), "\n")
cat("Default assay:", DefaultAssay(ref_seurat), "\n")
cat("Reference dataset dimensions:", dim(ref_seurat), "\n")

# Ensure reference is properly normalized (if not already)
if (!"data" %in% names(ref_seurat@assays$RNA@layers)) {
  cat("Normalizing reference dataset...\n")
  ref_seurat <- NormalizeData(ref_seurat, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
}

# Define conditions to process
conditions <- c("ctrl", "nacl", "racl", "xal")

# Function to process each condition
process_condition <- function(condition) {
  cat("Processing condition:", condition, "\n")
  
  # Load the integrated dataset
  input_file <- paste0("102_normal_integrated_pipeline/102_", condition, "_integrated.rds")
  seurat_obj <- readRDS(input_file)
  
  cat("Query dataset info:\n")
  cat("Assays available:", names(seurat_obj@assays), "\n")
  cat("Default assay:", DefaultAssay(seurat_obj), "\n")
  cat("Query dataset dimensions:", dim(seurat_obj), "\n")
  
  # Ensure query dataset is using RNA assay for transfer
  if (DefaultAssay(seurat_obj) != "RNA") {
    cat("Switching to RNA assay for transfer...\n")
    DefaultAssay(seurat_obj) <- "RNA"
  }
  
  # Ensure query dataset is properly normalized
  if (!"data" %in% names(seurat_obj@assays$RNA@layers)) {
    cat("Normalizing query dataset...\n")
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
  }
  
  # Find overlapping features
  overlap_features <- intersect(rownames(ref_seurat), rownames(seurat_obj))
  cat("Number of overlapping features:", length(overlap_features), "\n")
  
  # Ensure we have enough overlapping features
  if (length(overlap_features) < 100) {
    cat("Warning: Very few overlapping features found!\n")
  }
  
  # Find transfer anchors
  cat("Finding transfer anchors...\n")
  anchors <- FindTransferAnchors(
    reference = ref_seurat,
    query = seurat_obj,
    reference.reduction = "pca",
    dims = 1:30,
    features = overlap_features
  )
  
  # Transfer data
  cat("Transferring cell type labels...\n")
  predictions <- TransferData(
    anchorset = anchors, 
    refdata = ref_seurat$celltype_extended_atlas, 
    dims = 1:30
  )
  
  # Add metadata
  seurat_obj <- AddMetaData(seurat_obj, metadata = predictions)
  
  # Create UMAP plot with predicted labels
  p <- DimPlot(
    seurat_obj, 
    reduction = "umap", 
    group.by = "predicted.id", 
    label = TRUE, 
    label.size = 3, 
    repel = TRUE
  ) + NoLegend()
  
  # Save UMAP plot
  ggsave(
    p, 
    filename = file.path(output_dir, paste0("103_marioni_on_", condition, "_normal_integrated.pdf")), 
    height = 20, 
    width = 24
  )
  
  # Create comprehensive analysis PDF
  pdf(file.path(output_dir, paste0("103_marioni_on_", condition, "_conf_score.pdf")), 
      height = 25, width = 30)
  
  # Confidence score plot
  FeaturePlot(seurat_obj, features = "prediction.score.max")
  
  # Filter for high confidence predictions
  query_high_conf <- subset(seurat_obj, subset = prediction.score.max > 0.7)
  
  # Plot 1: Original clusters
  p1 <- DimPlot(seurat_obj, group.by = "seurat_clusters", label = TRUE) +
    ggtitle("Original Seurat Clusters")
  
  # Plot 2: Predicted labels (all cells)
  p2 <- DimPlot(seurat_obj, group.by = "predicted.id", label = TRUE, repel = TRUE) +
    ggtitle("Predicted Cell Types (All)") + NoLegend()
  
  # Plot 3: Predicted labels (confidence > 0.7)
  p3 <- DimPlot(query_high_conf, group.by = "predicted.id", label = TRUE, repel = TRUE) +
    ggtitle("Predicted Cell Types (Confidence > 0.7)")
  
  # Combine plots
  p_combined <- p1 + p2 + p3
  print(p_combined)
  
  # Prepare prediction scores for heatmap
  pred_scores <- seurat_obj@meta.data %>%
    select(seurat_clusters, starts_with("prediction.score.")) %>%
    pivot_longer(-seurat_clusters, names_to = "label", values_to = "score") %>%
    mutate(label = gsub("prediction.score.", "", label)) %>%
    group_by(seurat_clusters, label) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = label, values_from = mean_score)
  
  pred_scores <- as.data.frame(pred_scores)
  rownames(pred_scores) <- pred_scores[,'seurat_clusters']
  
  # Filter non-zero columns
  non_zero_cols <- apply(pred_scores, 2, function(x) any(x != 0))
  pred_scores_filtered <- pred_scores[, non_zero_cols]
  pred_scores_filtered_mat <- data.matrix(pred_scores_filtered)
  rownames(pred_scores_filtered_mat) <- rownames(pred_scores_filtered)
  
  # Create heatmaps
  col_fun <- colorRamp2(c(0, 0.5, 1), c("white", "red", "darkred"))
  
  # Full prediction scores heatmap
  Heatmap(
    pred_scores[,-1], 
    cluster_rows = diana(pred_scores[,-1]), 
    cluster_columns = agnes(t(pred_scores[,-1])), 
    col = col_fun
  )
  
  # Filtered prediction scores heatmap
  Heatmap(
    pred_scores_filtered_mat[,-1], 
    cluster_rows = diana(pred_scores_filtered_mat[,-1]),
    cluster_columns = agnes(t(pred_scores_filtered_mat[,-1])), 
    col = col_fun
  )
  
  # Cell type assignment heatmap
  col_fun2 <- colorRamp2(c(0, 50, 100), c("white", "red", "darkred"))
  assigned.states <- data.frame(
    dcast(
      data.frame(table(seurat_obj@meta.data[,c("seurat_clusters","predicted.id")])), 
      seurat_clusters ~ predicted.id
    )
  )
  rownames(assigned.states) <- assigned.states$seurat_clusters
  assigned.states[,-1] <- 100 * assigned.states[,-1] / rowSums(assigned.states[,-1])
  
  Heatmap(
    assigned.states[,-1], 
    cluster_rows = diana(assigned.states[,-1]), 
    cluster_columns = agnes(t(assigned.states[,-1])), 
    col = col_fun2
  )
  
  dev.off()
  
  # Create barplots PDF
  pdf(file.path(output_dir, paste0("103_marioni_on_", condition, "_barplots.pdf")), 
      height = 15, width = 20)
  
  # Barplot 1: Total cell counts per predicted cell type
  cell_type_counts <- table(seurat_obj$predicted.id)
  cell_type_counts_df <- data.frame(
    CellType = names(cell_type_counts),
    Count = as.numeric(cell_type_counts)
  ) %>%
    arrange(desc(Count))
  
  p1 <- ggplot(cell_type_counts_df, aes(x = reorder(CellType, Count), y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    geom_text(aes(label = Count), vjust = -0.5, size = 3) +
    coord_flip() +
    labs(title = paste0(toupper(condition), " - Total Cell Counts by Predicted Cell Type"),
         x = "Predicted Cell Type",
         y = "Number of Cells") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 14, face = "bold"))
  
  print(p1)
  
  # Barplot 2: Cell type distribution by Seurat clusters
  cluster_celltype_counts <- table(seurat_obj$seurat_clusters, seurat_obj$predicted.id)
  cluster_celltype_df <- as.data.frame(cluster_celltype_counts) %>%
    rename(Cluster = Var1, CellType = Var2, Count = Freq) %>%
    filter(Count > 0) %>%
    arrange(Cluster, desc(Count))
  
  p2 <- ggplot(cluster_celltype_df, aes(x = Cluster, y = Count, fill = CellType)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = paste0(toupper(condition), " - Cell Type Distribution by Seurat Clusters"),
         x = "Seurat Clusters",
         y = "Number of Cells",
         fill = "Predicted Cell Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 14, face = "bold"),
          legend.position = "right")
  
  print(p2)
  
  # Barplot 3: Percentage of each cell type per cluster
  cluster_celltype_percent <- cluster_celltype_df %>%
    group_by(Cluster) %>%
    mutate(Percentage = Count / sum(Count) * 100) %>%
    ungroup()
  
  p3 <- ggplot(cluster_celltype_percent, aes(x = Cluster, y = Percentage, fill = CellType)) +
    geom_bar(stat = "identity", position = "stack") +
    labs(title = paste0(toupper(condition), " - Cell Type Percentage by Seurat Clusters"),
         x = "Seurat Clusters",
         y = "Percentage (%)",
         fill = "Predicted Cell Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 14, face = "bold"),
          legend.position = "right")
  
  print(p3)
  
  # Barplot 4: Top cell types per cluster (top 5)
  top_celltypes_per_cluster <- cluster_celltype_df %>%
    group_by(Cluster) %>%
    top_n(5, Count) %>%
    ungroup()
  
  p4 <- ggplot(top_celltypes_per_cluster, aes(x = Cluster, y = Count, fill = CellType)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = paste0(toupper(condition), " - Top 5 Cell Types per Cluster"),
         x = "Seurat Clusters",
         y = "Number of Cells",
         fill = "Predicted Cell Type") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 14, face = "bold"),
          legend.position = "right")
  
  print(p4)
  
  # Barplot 5: Confidence score distribution
  confidence_df <- data.frame(
    Confidence = seurat_obj$prediction.score.max
  )
  
  p5 <- ggplot(confidence_df, aes(x = Confidence)) +
    geom_histogram(bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
    geom_vline(xintercept = 0.7, color = "red", linetype = "dashed", size = 1) +
    labs(title = paste0(toupper(condition), " - Distribution of Prediction Confidence Scores"),
         x = "Prediction Confidence Score",
         y = "Number of Cells") +
    theme_minimal() +
    theme(plot.title = element_text(size = 14, face = "bold"))
  
  print(p5)
  
  dev.off()
  
  # Save the processed object
  saveRDS(
    seurat_obj, 
    file = file.path(output_dir, paste0("103_", condition, "_marioni_gast_transfer.rds"))
  )
  
  cat("Completed processing for condition:", condition, "\n")
  return(seurat_obj)
}

# Process all conditions
results <- list()
for (condition in conditions) {
  tryCatch({
    results[[condition]] <- process_condition(condition)
  }, error = function(e) {
    cat("Error processing condition", condition, ":", e$message, "\n")
  })
}

# Create a summary report
cat("Processing completed for all conditions.\n")
cat("Output files saved in:", output_dir, "\n")
cat("Processed conditions:", paste(conditions, collapse = ", "), "\n") 