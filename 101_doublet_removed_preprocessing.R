library(sumoscrnaseq)
library(dplyr)
library(Seurat)
library(dittoSeq)
library(ggplot2)
library(gridExtra)
library(tidyr)
library(patchwork)

# Function to save table as PNG using base R
save_table_as_png <- function(data, filename, title = "Table") {
  # Create a temporary text file
  temp_file <- tempfile(fileext = ".txt")
  
  # Write table to text file
  sink(temp_file)
  cat(title, "\n")
  cat("=", paste(rep("=", nchar(title)), collapse = ""), "\n\n")
  print(data, row.names = FALSE)
  sink()
  
  # Read the text file and create a plot
  table_text <- readLines(temp_file)
  
  # Create PNG using base R graphics
  png(filename, width = 800, height = 600, res = 100)
  par(mar = c(1, 1, 1, 1))
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), 
       axes = FALSE, xlab = "", ylab = "", main = title)
  
  # Add text
  text(0.05, 0.95, paste(table_text, collapse = "\n"), 
       pos = 4, cex = 0.8, family = "mono", adj = c(0, 1))
  
  dev.off()
  
  # Clean up temp file
  unlink(temp_file)
}

# Function to perform preprocessing on doublet-removed datasets
preprocess_doublet_removed <- function(seurat_obj_path, output_dir = "101_doublet_removed_preprocessing") {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }
  
  # Load seurat object
  cat("Loading Seurat object from:", seurat_obj_path, "\n")
  seurat_obj <- readRDS(seurat_obj_path)
  
  # Extract basename for output naming
  seurat_basename <- paste0("101_", gsub("^100_", "", tools::file_path_sans_ext(basename(seurat_obj_path))))
  
  # Print initial object info
  cat("Initial Seurat object dimensions:", dim(seurat_obj), "\n")
  cat("Initial number of cells:", ncol(seurat_obj), "\n")
  
  # Step 1: Run process_seurat for preprocessing
  cat("Running process_seurat for preprocessing...\n")
  seurat_processed <- process_seurat(seurat_obj = seurat_obj)
  
  # Print processed object info
  cat("Processed Seurat object dimensions:", dim(seurat_processed), "\n")
  cat("Processed number of cells:", ncol(seurat_processed), "\n")
  
  # Step 2: Create QC plots
  cat("Creating QC plots...\n")
  qc_plots <- create_qc_plots(seurat_obj = seurat_processed)
  
  # Save QC plots
  qc_plots_pdf <- file.path(output_dir, paste0(seurat_basename, "_qc_plots.pdf"))
  pdf(qc_plots_pdf, width = 12, height = 10)
  
  # Print all QC plots
  for (i in seq_along(qc_plots)) {
    print(qc_plots[[i]])
  }
  
  dev.off()
  
  # Step 3: Create additional custom plots
  cat("Creating additional custom plots...\n")
  
  # Create comprehensive analysis PDF
  analysis_pdf <- file.path(output_dir, paste0(seurat_basename, "_comprehensive_analysis.pdf"))
  pdf(analysis_pdf, width = 12, height = 10)
  
  # Plot 1: UMAP with clusters
  p1 <- DimPlot(seurat_processed, reduction = "umap", group.by = "seurat_clusters") + 
    ggtitle("UMAP - Clusters")
  
  # Plot 2: UMAP with scDblFinder classification
  if ("scDblFinder_class" %in% colnames(seurat_processed@meta.data)) {
    p2 <- DimPlot(seurat_processed, reduction = "umap", group.by = "scDblFinder_class") + 
      ggtitle("UMAP - scDblFinder Classification")
  } else {
    p2 <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "scDblFinder_class not available") +
      theme_void() + ggtitle("UMAP - scDblFinder Classification")
  }
  
  # Plot 3: Feature plot of scDblFinder score
  if ("scDblFinder_score" %in% colnames(seurat_processed@meta.data)) {
    p3 <- FeaturePlot(seurat_processed, features = "scDblFinder_score") + 
      ggtitle("scDblFinder Score")
  } else {
    p3 <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "scDblFinder_score not available") +
      theme_void() + ggtitle("scDblFinder Score")
  }
  
  # Plot 4: QC metrics violin plots
  qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
  p4 <- VlnPlot(seurat_processed, features = qc_features, group.by = "seurat_clusters", 
                pt.size = 0.1, ncol = 3) + 
    plot_annotation(title = "QC Metrics by Cluster")
  
  # Plot 5: Cell cycle analysis if available
  if ("Phase" %in% colnames(seurat_processed@meta.data)) {
    p5 <- DimPlot(seurat_processed, reduction = "umap", group.by = "Phase") + 
      ggtitle("Cell Cycle Phase")
  } else {
    p5 <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Cell cycle data not available") +
      theme_void() + ggtitle("Cell Cycle Phase")
  }
  
  # Plot 6: Feature plot of top variable genes
  top_genes <- head(VariableFeatures(seurat_processed), 6)
  if (length(top_genes) > 0) {
    p6 <- FeaturePlot(seurat_processed, features = top_genes, ncol = 3) + 
      plot_annotation(title = "Top Variable Genes")
  } else {
    p6 <- ggplot() + annotate("text", x = 0.5, y = 0.5, label = "Variable genes not available") +
      theme_void() + ggtitle("Top Variable Genes")
  }
  
  # Print all plots
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p5)
  print(p6)
  
  dev.off()
  
  # Step 4: Create Phase barplot using dittoSeq
  cat("Creating Phase barplot using dittoSeq...\n")
  
  # Create Phase barplot PDF
  phase_barplot_pdf <- file.path(output_dir, paste0(seurat_basename, "_phase_barplot.pdf"))
  pdf(phase_barplot_pdf, width = 10, height = 8)
  
  if ("Phase" %in% colnames(seurat_processed@meta.data)) {
    # Create barplot using dittoSeq
    phase_barplot <- dittoBarPlot(
      seurat_processed, 
      var = "Phase", 
      group.by = "seurat_clusters",
      scale = "count",  # or "percent" for percentage
      color.panel = c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7"),
      main = "Cell Cycle Phase Distribution by Cluster",
      xlab = "Cluster",
      ylab = "Number of Cells"
    )
    print(phase_barplot)
    
    # Also create percentage version
    phase_barplot_percent <- dittoBarPlot(
      seurat_processed, 
      var = "Phase", 
      group.by = "seurat_clusters",
      scale = "percent",
      color.panel = c("#FF6B6B", "#4ECDC4", "#45B7D1", "#96CEB4", "#FFEAA7"),
      main = "Cell Cycle Phase Distribution by Cluster (%)",
      xlab = "Cluster",
      ylab = "Percentage of Cells"
    )
    print(phase_barplot_percent)
    
  } else {
    # Create placeholder plot if Phase data is not available
    placeholder_plot <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = "Cell cycle data not available", size = 6) +
      theme_void() + 
      ggtitle("Cell Cycle Phase Distribution by Cluster")
    print(placeholder_plot)
  }
  
  dev.off()
  
  # Step 5: Create summary statistics
  cat("Creating summary statistics...\n")
  
  # Summary statistics table
  summary_stats <- data.frame(
    Metric = c("Total cells", "Total genes", "Mean genes per cell", "Mean UMIs per cell", 
               "Median genes per cell", "Median UMIs per cell", "Number of clusters"),
    Value = c(
      ncol(seurat_processed),
      nrow(seurat_processed),
      round(mean(seurat_processed$nFeature_RNA), 1),
      round(mean(seurat_processed$nCount_RNA), 1),
      round(median(seurat_processed$nFeature_RNA), 1),
      round(median(seurat_processed$nCount_RNA), 1),
      length(unique(seurat_processed$seurat_clusters))
    )
  )
  
  # Save summary statistics table using base R
  summary_table_png <- file.path(output_dir, paste0(seurat_basename, "_summary_statistics.png"))
  save_table_as_png(summary_stats, summary_table_png, "Preprocessing Summary Statistics")
  
  # Step 6: Create cluster statistics
  cluster_stats <- seurat_processed@meta.data %>%
    group_by(seurat_clusters) %>%
    summarise(
      cell_count = n(),
      mean_genes = round(mean(nFeature_RNA), 1),
      mean_umis = round(mean(nCount_RNA), 1),
      mean_mt = round(mean(percent.mt), 2)
    ) %>%
    arrange(seurat_clusters)
  
  # Save cluster statistics table using base R
  cluster_table_png <- file.path(output_dir, paste0(seurat_basename, "_cluster_statistics.png"))
  save_table_as_png(cluster_stats, cluster_table_png, "Cluster Statistics")
  
  # Step 7: Create Phase statistics if available
  if ("Phase" %in% colnames(seurat_processed@meta.data)) {
    phase_stats <- seurat_processed@meta.data %>%
      group_by(seurat_clusters, Phase) %>%
      summarise(count = n(), .groups = 'drop') %>%
      pivot_wider(names_from = Phase, values_from = count, values_fill = 0) %>%
      mutate(total_cells = rowSums(across(-seurat_clusters))) %>%
      arrange(seurat_clusters)
    
    # Save Phase statistics table using base R
    phase_table_png <- file.path(output_dir, paste0(seurat_basename, "_phase_statistics.png"))
    save_table_as_png(phase_stats, phase_table_png, "Cell Cycle Phase Statistics by Cluster")
  }
  
  # Step 8: Save processed Seurat object
  output_rds <- file.path(output_dir, paste0(seurat_basename, "_processed.rds"))
  saveRDS(seurat_processed, output_rds)
  
  # Step 9: Create summary text file
  summary_txt <- file.path(output_dir, paste0(seurat_basename, "_preprocessing_summary.txt"))
  sink(summary_txt)
  cat("Preprocessing Summary\n")
  cat("====================\n\n")
  cat("Input file:", seurat_obj_path, "\n")
  cat("Output directory:", output_dir, "\n\n")
  
  cat("Initial Seurat object dimensions:", dim(seurat_obj), "\n")
  cat("Initial number of cells:", ncol(seurat_obj), "\n\n")
  
  cat("Processed Seurat object dimensions:", dim(seurat_processed), "\n")
  cat("Processed number of cells:", ncol(seurat_processed), "\n\n")
  
  cat("Summary Statistics:\n")
  print(summary_stats)
  cat("\n")
  
  cat("Cluster Statistics:\n")
  print(cluster_stats)
  cat("\n")
  
  if ("Phase" %in% colnames(seurat_processed@meta.data)) {
    cat("Phase Statistics:\n")
    print(phase_stats)
    cat("\n")
  }
  
  cat("Output files:\n")
  cat("- Processed Seurat object:", output_rds, "\n")
  cat("- QC plots:", qc_plots_pdf, "\n")
  cat("- Comprehensive analysis:", analysis_pdf, "\n")
  cat("- Phase barplot:", phase_barplot_pdf, "\n")
  cat("- Summary statistics table:", summary_table_png, "\n")
  cat("- Cluster statistics table:", cluster_table_png, "\n")
  if ("Phase" %in% colnames(seurat_processed@meta.data)) {
    cat("- Phase statistics table:", phase_table_png, "\n")
  }
  cat("- Summary text file:", summary_txt, "\n")
  sink()
  
  cat("Preprocessing complete!\n")
  cat("Processed Seurat object saved to:", output_rds, "\n")
  cat("QC plots saved to:", qc_plots_pdf, "\n")
  cat("Comprehensive analysis saved to:", analysis_pdf, "\n")
  cat("Phase barplot saved to:", phase_barplot_pdf, "\n")
  cat("Summary statistics saved to:", summary_table_png, "\n")
  cat("Cluster statistics saved to:", cluster_table_png, "\n")
  if ("Phase" %in% colnames(seurat_processed@meta.data)) {
    cat("Phase statistics saved to:", phase_table_png, "\n")
  }
  cat("Summary text file saved to:", summary_txt, "\n")
  
  return(list(
    seurat_obj = seurat_processed,
    summary_stats = summary_stats,
    cluster_stats = cluster_stats,
    phase_stats = if ("Phase" %in% colnames(seurat_processed@meta.data)) phase_stats else NULL,
    output_files = list(
      rds = output_rds,
      qc_plots = qc_plots_pdf,
      analysis = analysis_pdf,
      phase_barplot = phase_barplot_pdf,
      summary_table = summary_table_png,
      cluster_table = cluster_table_png,
      phase_table = if ("Phase" %in% colnames(seurat_processed@meta.data)) phase_table_png else NULL,
      summary_txt = summary_txt
    )
  ))
}

# List of doublet-removed files to process
doublet_removed_files <- c(
  "100_doublet_removed/100_sc_ctrl_reg_cc_clusters_filtered_doublets_removed.rds",
  "100_doublet_removed/100_sc_nacl_reg_cc_clusters_filtered_doublets_removed.rds",
  "100_doublet_removed/100_sc_racl_reg_cc_clusters_filtered_doublets_removed.rds",
  "100_doublet_removed/100_sc_xal_reg_cc_clusters_filtered_doublets_removed.rds"
)

# Process all files
results <- list()
for (file_path in doublet_removed_files) {
  if (file.exists(file_path)) {
    cat("Processing:", file_path, "\n")
    condition <- gsub("100_sc_|_reg_cc_clusters_filtered_doublets_removed.rds", "", basename(file_path))
    results[[condition]] <- preprocess_doublet_removed(file_path)
    cat("Completed processing for", condition, "\n\n")
  } else {
    cat("File not found:", file_path, "\n")
  }
}

# Create overall summary
cat("All preprocessing completed!\n")
cat("Total conditions processed:", length(results), "\n")
cat("Conditions:", paste(names(results), collapse = ", "), "\n") 