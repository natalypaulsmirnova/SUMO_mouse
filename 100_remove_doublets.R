library(sumoscrnaseq)
library(dplyr)
library(Seurat)
library(dittoSeq)
library(ggplot2)
library(gridExtra)
library(kableExtra)
library(tidyr)
library(patchwork)

# Function to remove doublets based on scDblFinder score
remove_doublets <- function(seurat_obj_path, output_dir = "100_doublet_removed", score_threshold = 0.7) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  # Load seurat object
  cat("Loading Seurat object from:", seurat_obj_path, "\n")
  seurat_obj <- readRDS(seurat_obj_path)
  
  # Extract basename for output naming
  seurat_basename <- paste0("100_", gsub("^[0-9]+_", "", tools::file_path_sans_ext(basename(seurat_obj_path))))
  
  # Print initial object info
  cat("Initial Seurat object dimensions:", dim(seurat_obj), "\n")
  cat("Initial number of cells:", ncol(seurat_obj), "\n")
  
  # Run scDblFinder if not already run
  if (!"scDblFinder_score" %in% colnames(seurat_obj@meta.data)) {
    cat("Running scDblFinder...\n")
    seurat_obj <- run_scdblfinder(seurat_obj = seurat_obj)
  } else {
    cat("scDblFinder already run, using existing scores.\n")
  }
  
  # Print doublet statistics before filtering
  doublet_stats_before <- seurat_obj@meta.data %>%
    group_by(scDblFinder_class) %>%
    summarise(
      count = n(),
      percentage = round(n() / nrow(.) * 100, 2),
      mean_score = round(mean(scDblFinder_score, na.rm = TRUE), 3),
      median_score = round(median(scDblFinder_score, na.rm = TRUE), 3)
    )
  
  cat("Doublet statistics before filtering:\n")
  print(doublet_stats_before)
  
  # Filter out doublets with score above threshold
  cat("Filtering out doublets with score >", score_threshold, "\n")
  seurat_obj_filtered <- subset(seurat_obj, scDblFinder_score <= score_threshold)
  
  # Print filtering results
  cat("Cells removed:", ncol(seurat_obj) - ncol(seurat_obj_filtered), "\n")
  cat("Cells remaining:", ncol(seurat_obj_filtered), "\n")
  cat("Percentage of cells kept:", round(ncol(seurat_obj_filtered) / ncol(seurat_obj) * 100, 2), "%\n")
  
  # Print doublet statistics after filtering
  doublet_stats_after <- seurat_obj_filtered@meta.data %>%
    group_by(scDblFinder_class) %>%
    summarise(
      count = n(),
      percentage = round(n() / nrow(.) * 100, 2),
      mean_score = round(mean(scDblFinder_score, na.rm = TRUE), 3),
      median_score = round(median(scDblFinder_score, na.rm = TRUE), 3)
    )
  
  cat("Doublet statistics after filtering:\n")
  print(doublet_stats_after)
  
  # Create comprehensive PDF with filtering results
  output_pdf <- file.path(output_dir, paste0(seurat_basename, "_doublet_removal_summary.pdf"))
  pdf(output_pdf, width = 12, height = 10)
  
  # Plot 1: UMAP before filtering
  p1 <- DimPlot(seurat_obj, group.by = "scDblFinder_class") + 
    ggtitle("Before Filtering: scDblFinder Classification")
  
  # Plot 2: UMAP after filtering
  p2 <- DimPlot(seurat_obj_filtered, group.by = "scDblFinder_class") + 
    ggtitle("After Filtering: scDblFinder Classification")
  
  # Plot 3: Feature plot of scDblFinder score before filtering
  p3 <- FeaturePlot(seurat_obj, features = "scDblFinder_score") + 
    ggtitle("Before Filtering: scDblFinder Score")
  
  # Plot 4: Feature plot of scDblFinder score after filtering
  p4 <- FeaturePlot(seurat_obj_filtered, features = "scDblFinder_score") + 
    ggtitle("After Filtering: scDblFinder Score")
  
  # Plot 5: Violin plot of scDblFinder score before filtering
  p5 <- VlnPlot(seurat_obj, features = "scDblFinder_score", group.by = "scDblFinder_class") + 
    ggtitle("Before Filtering: scDblFinder Score Distribution")
  
  # Plot 6: Violin plot of scDblFinder score after filtering
  p6 <- VlnPlot(seurat_obj_filtered, features = "scDblFinder_score", group.by = "scDblFinder_class") + 
    ggtitle("After Filtering: scDblFinder Score Distribution")
  
  # Plot 7: Histogram of scDblFinder scores
  score_data <- data.frame(
    score = seurat_obj@meta.data$scDblFinder_score,
    class = seurat_obj@meta.data$scDblFinder_class,
    filtered = "Before"
  )
  
  score_data_filtered <- data.frame(
    score = seurat_obj_filtered@meta.data$scDblFinder_score,
    class = seurat_obj_filtered@meta.data$scDblFinder_class,
    filtered = "After"
  )
  
  score_data_combined <- rbind(score_data, score_data_filtered)
  
  p7 <- ggplot(score_data_combined, aes(x = score, fill = class)) +
    geom_histogram(bins = 50, alpha = 0.7) +
    facet_wrap(~filtered, scales = "free_y") +
    geom_vline(xintercept = score_threshold, color = "red", linetype = "dashed") +
    labs(title = "Distribution of scDblFinder Scores",
         x = "scDblFinder Score",
         y = "Count") +
    theme_minimal()
  
  # Plot 8: QC metrics comparison
  qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  
  # Before filtering
  qc_before <- VlnPlot(
    seurat_obj,
    features = qc_features,
    group.by = "scDblFinder_class",
    pt.size = 0.1,
    ncol = 3
  ) + plot_annotation(title = "QC Metrics Before Filtering")
  
  # After filtering
  qc_after <- VlnPlot(
    seurat_obj_filtered,
    features = qc_features,
    group.by = "scDblFinder_class",
    pt.size = 0.1,
    ncol = 3
  ) + plot_annotation(title = "QC Metrics After Filtering")
  
  # Print all plots
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  print(p6)
  print(p7)
  print(qc_before)
  print(qc_after)
  
  dev.off()
  
  # Create summary tables
  # Table 1: Summary statistics
  summary_stats <- data.frame(
    Metric = c("Total cells before", "Total cells after", "Cells removed", "Percentage kept"),
    Value = c(
      ncol(seurat_obj),
      ncol(seurat_obj_filtered),
      ncol(seurat_obj) - ncol(seurat_obj_filtered),
      paste0(round(ncol(seurat_obj_filtered) / ncol(seurat_obj) * 100, 2), "%")
    )
  )
  
  # Create and save summary table
  summary_table_png <- file.path(output_dir, paste0(seurat_basename, "_summary_statistics.png"))
  summary_stats %>%
    kbl(caption = "Doublet Removal Summary Statistics") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) %>%
    column_spec(1, bold = TRUE) %>%
    save_kable(summary_table_png, zoom = 2)
  
  # Table 2: Doublet statistics by cluster before filtering
  doublet_by_cluster_before <- seurat_obj@meta.data %>%
    group_by(seurat_clusters, scDblFinder_class) %>%
    summarise(count = n(), .groups = 'drop') %>%
    pivot_wider(names_from = scDblFinder_class, values_from = count, values_fill = 0) %>%
    mutate(total_cells = singlet + doublet,
           doublet_percentage = round((doublet / total_cells) * 100, 2)) %>%
    arrange(seurat_clusters)
  
  # Create and save doublet by cluster table
  doublet_cluster_table_png <- file.path(output_dir, paste0(seurat_basename, "_doublet_by_cluster_before_filtering.png"))
  doublet_by_cluster_before %>%
    kbl(caption = "Doublet Distribution by Cluster Before Filtering") %>%
    kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE) %>%
    column_spec(1, bold = TRUE) %>%
    column_spec(5, color = "white", background = "#D7261E") %>%
    save_kable(doublet_cluster_table_png, zoom = 2)
  
  # Save the filtered Seurat object
  output_rds <- file.path(output_dir, paste0(seurat_basename, "_doublets_removed.rds"))
  saveRDS(seurat_obj_filtered, output_rds)
  
  # Save summary statistics as text file
  summary_txt <- file.path(output_dir, paste0(seurat_basename, "_doublet_removal_summary.txt"))
  sink(summary_txt)
  cat("Doublet Removal Summary\n")
  cat("======================\n\n")
  cat("Input file:", seurat_obj_path, "\n")
  cat("Output directory:", output_dir, "\n")
  cat("Score threshold:", score_threshold, "\n\n")
  
  cat("Initial Seurat object dimensions:", dim(seurat_obj), "\n")
  cat("Initial number of cells:", ncol(seurat_obj), "\n\n")
  
  cat("Doublet statistics before filtering:\n")
  print(doublet_stats_before)
  cat("\n")
  
  cat("Cells removed:", ncol(seurat_obj) - ncol(seurat_obj_filtered), "\n")
  cat("Cells remaining:", ncol(seurat_obj_filtered), "\n")
  cat("Percentage of cells kept:", round(ncol(seurat_obj_filtered) / ncol(seurat_obj) * 100, 2), "%\n\n")
  
  cat("Doublet statistics after filtering:\n")
  print(doublet_stats_after)
  cat("\n")
  
  cat("Output files:\n")
  cat("- Filtered Seurat object:", output_rds, "\n")
  cat("- Summary plots:", output_pdf, "\n")
  cat("- Summary statistics table:", summary_table_png, "\n")
  cat("- Doublet by cluster table:", doublet_cluster_table_png, "\n")
  cat("- Summary text file:", summary_txt, "\n")
  sink()
  
  cat("Analysis complete!\n")
  cat("Filtered Seurat object saved to:", output_rds, "\n")
  cat("Summary plots saved to:", output_pdf, "\n")
  cat("Summary statistics saved to:", summary_table_png, "\n")
  cat("Doublet by cluster table saved to:", doublet_cluster_table_png, "\n")
  cat("Summary text file saved to:", summary_txt, "\n")
  
  return(seurat_obj_filtered)
}

# Example usage:
# filtered_obj <- remove_doublets("05_sc_racl_reg_cc_clusters_filtered.rds", 
#                                 output_dir = "100_doublet_removed", 
#                                 score_threshold = 0.5) 