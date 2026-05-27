library(Seurat)
library(pheatmap)
library(ggbiplot)
library(reshape2)
library(ggplot2)
library(dplyr)
library(tidyr)
library(cluster)
library(ComplexHeatmap)
library(circlize)

# Define conditions and create output directory
conditions <- c("ctrl", "nacl", "racl", "xal")
output_dir <- "104_EDA_marioni_on_condition_normal_integrated"

# Check if output directory exists and show warning if it does
if (dir.exists(output_dir)) {
  warning("Output directory '", output_dir, "' already exists. Files may be overwritten.")
} else {
  dir.create(output_dir, showWarnings = FALSE)
}

# Function to perform EDA analysis for a single condition
perform_eda_analysis <- function(condition) {
  cat("Processing condition:", condition, "\n")
  
  # Load the integrated seurat object
  input_file <- paste0("103_marioni_on_condition_normal_integrated/103_", condition, "_marioni_gast_transfer.rds")
  seurat_obj <- readRDS(input_file)
  
  # Arrange the time factor as 48h, 72h, 96h and 120h
  seurat_obj$time <- factor(seurat_obj$time, levels = c("48h", "72h", "96h", "120h"))
  
  # Create condition-specific output directory
  condition_output_dir <- file.path(output_dir, condition)
  
  # Check if condition directory exists and show warning if it does
  if (dir.exists(condition_output_dir)) {
    warning("Condition directory '", condition_output_dir, "' already exists. Files may be overwritten.")
  } else {
    dir.create(condition_output_dir, showWarnings = FALSE)
  }
  
  # 1. Prediction Score Analysis
  cat("  - Generating prediction score plots...\n")
  
  # Violin plot of prediction score
  vln_plot <- VlnPlot(seurat_obj, features = "prediction.score.max", 
                      group.by = "predicted.id", pt.size = 0.1) + 
    NoLegend() +
    ggtitle(paste0("Prediction Score Distribution - ", toupper(condition)))
  
  ggsave(file.path(condition_output_dir, paste0("104_", condition, "_prediction_score_violin.pdf")), 
         vln_plot, width = 12, height = 8)
  
  # Feature plot of prediction score
  feat_plot <- FeaturePlot(seurat_obj, features = "prediction.score.max", 
                          cols = c("lightgrey", "blue")) +
    ggtitle(paste0("Prediction Score Spatial Distribution - ", toupper(condition)))
  
  ggsave(file.path(condition_output_dir, paste0("104_", condition, "_prediction_score_feature.pdf")), 
         feat_plot, width = 10, height = 8)
  
  # 2. Chi-squared Analysis
  cat("  - Performing chi-squared analysis...\n")
  
  # Chi-squared test on predicted cell types vs time
  ct_table <- table(seurat_obj$predicted.id, seurat_obj$time)
  chi <- chisq.test(ct_table)
  
  # Save chi-squared results
  chi_results <- data.frame(
    statistic = chi$statistic,
    p_value = chi$p.value,
    df = chi$parameter
  )
  write.csv(chi_results, file.path(condition_output_dir, paste0("104_", condition, "_chi_squared_results.csv")), 
            row.names = FALSE)
  
  # Force residuals into proper matrix
  res_array <- chi$residuals
  res_mat <- matrix(as.numeric(res_array),
                    nrow = dim(res_array)[1],
                    dimnames = dimnames(res_array))
  
  # Identify top N most dynamic cell types (based on residual SD)
  res_sd <- apply(res_mat, 1, sd)
  top_n <- 12
  top_dynamic <- names(sort(res_sd, decreasing = TRUE)[1:top_n])
  
  # Subset and reshape to long format
  res_mat_top <- res_mat[top_dynamic, , drop = FALSE]
  res_df <- melt(res_mat_top)
  colnames(res_df) <- c("celltype", "time", "residual")
  
  # Ensure time is ordered correctly
  res_df$time <- factor(as.character(res_df$time), levels = c("48h", "72h", "96h", "120h"))
  
  # Plot residual trajectories
  residual_plot <- ggplot(res_df, aes(x = time, y = residual, group = celltype)) +
    geom_line(linewidth = 1.1, color = "#2C3E50") +
    geom_point(size = 2, color = "#E74C3C") +
    facet_wrap(~ celltype, scales = "free_y", ncol = 4) +
    theme_minimal(base_size = 14) +
    labs(title = paste0("Chi-squared Residual Trajectories - ", toupper(condition)),
         x = "Time", y = "Chi-squared Residual") +
    theme(strip.text = element_text(size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(condition_output_dir, paste0("104_", condition, "_residual_trajectories.pdf")), 
         residual_plot, width = 16, height = 12)
  
  # Heatmap of top dynamic cell types
  pheatmap_plot <- pheatmap(res_mat_top, cluster_rows = TRUE, cluster_cols = TRUE,
                           main = paste0("Top Dynamic Cell Types - ", toupper(condition)))
  
  pdf(file.path(condition_output_dir, paste0("104_", condition, "_dynamic_celltypes_heatmap.pdf")), 
      width = 10, height = 8)
  print(pheatmap_plot)
  dev.off()
  
  # PCA plot of top dynamic cell types
  pca <- prcomp(res_mat_top, scale. = TRUE, center = TRUE)
  
  pca_plot <- ggbiplot(pca, labels = rownames(res_mat_top), ellipse = TRUE, circle = TRUE) +
    ggtitle(paste0("PCA of Chi-squared Residual Profiles - ", toupper(condition))) +
    theme_minimal()
  
  ggsave(file.path(condition_output_dir, paste0("104_", condition, "_pca_residuals.pdf")), 
         pca_plot, width = 10, height = 8)
  
  # Chi-squared test on seurat clusters vs predicted cell types
  cat("  - Performing chi-squared analysis for seurat clusters vs predicted cell types...\n")
  
  # Chi-squared test on seurat clusters vs predicted cell types
  cluster_table <- table(seurat_obj$seurat_clusters, seurat_obj$predicted.id)
  chi_cluster <- chisq.test(cluster_table)
  
  # Save chi-squared results for clusters
  chi_cluster_results <- data.frame(
    statistic = chi_cluster$statistic,
    p_value = chi_cluster$p.value,
    df = chi_cluster$parameter
  )
  write.csv(chi_cluster_results, file.path(condition_output_dir, paste0("104_", condition, "_chi_squared_clusters_results.csv")), 
            row.names = FALSE)
  
  # Force residuals into proper matrix for clusters
  res_cluster_array <- chi_cluster$residuals
  res_cluster_mat <- matrix(as.numeric(res_cluster_array),
                           nrow = dim(res_cluster_array)[1],
                           dimnames = dimnames(res_cluster_array))
  
  # Identify top N most dynamic clusters (based on residual SD)
  res_cluster_sd <- apply(res_cluster_mat, 1, sd)
  top_n_cluster <- 12
  top_dynamic_clusters <- names(sort(res_cluster_sd, decreasing = TRUE)[1:top_n_cluster])
  
  # Subset and reshape to long format for clusters
  res_cluster_mat_top <- res_cluster_mat[top_dynamic_clusters, , drop = FALSE]
  res_cluster_df <- melt(res_cluster_mat_top)
  colnames(res_cluster_df) <- c("cluster", "celltype", "residual")
  
  # Plot residual trajectories for clusters
  residual_cluster_plot <- ggplot(res_cluster_df, aes(x = celltype, y = residual, group = cluster)) +
    geom_line(linewidth = 1.1, color = "#2C3E50") +
    geom_point(size = 2, color = "#E74C3C") +
    facet_wrap(~ cluster, scales = "free_y", ncol = 4) +
    theme_minimal(base_size = 14) +
    labs(title = paste0("Chi-squared Residual Profiles (Clusters vs Cell Types) - ", toupper(condition)),
         x = "Predicted Cell Type", y = "Chi-squared Residual") +
    theme(strip.text = element_text(size = 11),
          axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(file.path(condition_output_dir, paste0("104_", condition, "_residual_trajectories_clusters.pdf")), 
         residual_cluster_plot, width = 20, height = 12)
  
  # Heatmap of top dynamic clusters
  pheatmap_cluster_plot <- pheatmap(res_cluster_mat_top, cluster_rows = TRUE, cluster_cols = TRUE,
                                   main = paste0("Top Dynamic Clusters vs Cell Types - ", toupper(condition)))
  
  pdf(file.path(condition_output_dir, paste0("104_", condition, "_dynamic_clusters_heatmap.pdf")), 
      width = 12, height = 8)
  print(pheatmap_cluster_plot)
  dev.off()
  
  # PCA plot of top dynamic clusters
  pca_cluster <- prcomp(res_cluster_mat_top, scale. = TRUE, center = TRUE)
  
  pca_cluster_plot <- ggbiplot(pca_cluster, labels = rownames(res_cluster_mat_top), ellipse = TRUE, circle = TRUE) +
    ggtitle(paste0("PCA of Chi-squared Residual Profiles (Clusters vs Cell Types) - ", toupper(condition))) +
    theme_minimal()
  
  ggsave(file.path(condition_output_dir, paste0("104_", condition, "_pca_residuals_clusters.pdf")), 
         pca_cluster_plot, width = 10, height = 8)
  
  # 3. Cell Type Proportions Analysis
  cat("  - Analyzing cell type proportions...\n")
  
  # Compute proportions per cell type × time
  prop_df <- seurat_obj@meta.data %>%
    group_by(time, celltype = predicted.id) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(time) %>%
    mutate(total = sum(n), freq = n / total)
  
  # Save proportions to CSV
  write.csv(prop_df, file.path(condition_output_dir, paste0("104_", condition, "_celltype_proportions.csv")), 
            row.names = FALSE)
  
  # Filter for cell types with ≥1% frequency at any time point
  filtered_df <- prop_df %>%
    filter(freq >= 0.01)
  
  # Clean up time ordering
  filtered_df$time <- factor(filtered_df$time, levels = c("48h", "72h", "96h", "120h"))
  
  # Plot proportions
  prop_plot <- ggplot(filtered_df, aes(x = time, y = freq, group = 1)) +
    geom_line(linewidth = 1.1, color = "#2C3E50") +
    geom_point(size = 2, color = "#E74C3C") +
    facet_wrap(~ celltype, scales = "free_y", ncol = 4) +
    theme_minimal(base_size = 14) +
    labs(title = paste0("Cell Type Proportions Over Time - ", toupper(condition)),
         x = "Time", y = "Proportion of Cells") +
    theme(strip.text = element_text(size = 11))
  
  ggsave(file.path(condition_output_dir, paste0("104_", condition, "_celltype_proportions.pdf")), 
         prop_plot, width = 16, height = 12)
  
  # Sort proportions by time and frequency
  prop_df$time <- factor(prop_df$time, levels = c("48h", "72h", "96h", "120h"))
  prop_df_sorted <- prop_df %>%
    group_by(time) %>%
    arrange(time, desc(freq), .by_group = TRUE)
  
  # Save sorted proportions
  write.csv(prop_df_sorted, file.path(condition_output_dir, paste0("104_", condition, "_celltype_proportions_sorted.csv")), 
            row.names = FALSE)
  
  # 4. Prediction Score Analysis
  cat("  - Analyzing prediction scores...\n")
  
  # Prepare prediction score data
  pred_scores <- seurat_obj@meta.data %>%
    select(seurat_clusters, starts_with("prediction.score.")) %>%
    pivot_longer(-seurat_clusters, names_to = "label", values_to = "score") %>%
    mutate(label = gsub("prediction.score.", "", label)) %>%
    group_by(seurat_clusters, label) %>%
    summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = label, values_from = mean_score)
  
  pred_scores <- as.data.frame(pred_scores)
  rownames(pred_scores) <- pred_scores[,'seurat_clusters']
  
  # Save prediction scores
  write.csv(pred_scores, file.path(condition_output_dir, paste0("104_", condition, "_prediction_scores.csv")), 
            row.names = TRUE)
  
  # Filter non-zero columns
  non_zero_cols <- apply(pred_scores, 2, function(x) any(x != 0))
  pred_scores_filtered <- pred_scores[, non_zero_cols]
  pred_scores_filtered_mat <- data.matrix(pred_scores_filtered)
  rownames(pred_scores_filtered_mat) <- rownames(pred_scores_filtered)
  
  # Create heatmaps
  col_fun = colorRamp2(c(0, 0.5, 1), c("white", "red", "darkred"))
  
  # Full prediction scores heatmap
  pdf(file.path(condition_output_dir, paste0("104_", condition, "_prediction_scores_heatmap.pdf")), 
      width = 12, height = 10)
  print(Heatmap(pred_scores[,-1], 
                cluster_rows = diana(pred_scores[,-1]), 
                cluster_columns = agnes(t(pred_scores[,-1])), 
                col = col_fun,
                column_title = paste0("Prediction Scores - ", toupper(condition))))
  dev.off()
  
  # Filtered prediction scores heatmap
  pdf(file.path(condition_output_dir, paste0("104_", condition, "_prediction_scores_filtered_heatmap.pdf")), 
      width = 12, height = 10)
  print(Heatmap(pred_scores_filtered_mat[,-1], 
                cluster_rows = diana(pred_scores_filtered_mat[,-1]),
                cluster_columns = agnes(t(pred_scores_filtered_mat[,-1])), 
                col = col_fun,
                column_title = paste0("Filtered Prediction Scores - ", toupper(condition))))
  dev.off()
  
  # 5. Assigned States Analysis
  cat("  - Analyzing assigned states...\n")
  
  col_fun2 = colorRamp2(c(0, 50, 100), c("white", "red", "darkred"))
  assigned.states <- data.frame(dcast(data.frame(table(seurat_obj@meta.data[,c(10,19)])), seurat_clusters~predicted.id))
  rownames(assigned.states) <- assigned.states$seurat_clusters
  assigned.states[,-1] <- 100*assigned.states[,-1]/rowSums(assigned.states[,-1])
  
  # Save assigned states
  write.csv(assigned.states, file.path(condition_output_dir, paste0("104_", condition, "_assigned_states.csv")), 
            row.names = TRUE)
  
  # Assigned states heatmap
  pdf(file.path(condition_output_dir, paste0("104_", condition, "_assigned_states_heatmap.pdf")), 
      width = 12, height = 10)
  print(Heatmap(assigned.states[,-1], 
                cluster_rows = diana(assigned.states[,-1]), 
                cluster_columns = agnes(t(assigned.states[,-1])), 
                col = col_fun2,
                column_title = paste0("Assigned States - ", toupper(condition))))
  dev.off()
  
  cat("  - Completed analysis for", condition, "\n")
  
  # Return summary statistics
  return(list(
    condition = condition,
    total_cells = ncol(seurat_obj),
    unique_celltypes = length(unique(seurat_obj$predicted.id)),
    unique_clusters = length(unique(seurat_obj$seurat_clusters)),
    chi_squared_p = chi$p.value,
    chi_squared_clusters_p = chi_cluster$p.value,
    top_dynamic_celltypes = top_dynamic,
    top_dynamic_clusters = top_dynamic_clusters
  ))
}

# Process all conditions
cat("Starting EDA analysis for all conditions...\n")
results <- list()

for (condition in conditions) {
  tryCatch({
    results[[condition]] <- perform_eda_analysis(condition)
  }, error = function(e) {
    cat("Error processing condition", condition, ":", e$message, "\n")
  })
}

# Create summary report
cat("Creating summary report...\n")
summary_df <- do.call(rbind, lapply(results, function(x) {
  if (!is.null(x)) {
    data.frame(
      Condition = x$condition,
      Total_Cells = x$total_cells,
      Unique_Celltypes = x$unique_celltypes,
      Unique_Clusters = x$unique_clusters,
      Chi_Squared_Celltypes_vs_Time_P_Value = x$chi_squared_p,
      Chi_Squared_Clusters_vs_Celltypes_P_Value = x$chi_squared_clusters_p,
      Top_Dynamic_Celltypes = paste(x$top_dynamic_celltypes, collapse = "; "),
      Top_Dynamic_Clusters = paste(x$top_dynamic_clusters, collapse = "; ")
    )
  }
}))

if (!is.null(summary_df)) {
  write.csv(summary_df, file.path(output_dir, "104_analysis_summary.csv"), row.names = FALSE)
  
  # Create summary plot
  summary_plot <- ggplot(summary_df, aes(x = Condition, y = Total_Cells, fill = Condition)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Total Cells per Condition",
         x = "Condition", y = "Total Cells") +
    theme(legend.position = "none")
  
  ggsave(file.path(output_dir, "104_summary_total_cells.pdf"), summary_plot, width = 8, height = 6)
}

cat("EDA analysis completed! Results saved in:", output_dir, "\n")


