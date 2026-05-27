library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)

# Create output directory
output_dir <- "601_manual_labeling_marioni_filtering_time_based"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Define time points to process
time_points <- c("48h", "72h", "96h", "120h")

# Define cell types to remove
cell_types_to_remove <- c(
  "Amniotic ectoderm",
  "Midbrain progenitors",
  "Midbrain/Hindbrain boundary",
  "Midbrain/Hindbrain boundary",
  "Midbrain/Hindbrain border",
  "Neural crest",
  "Non-neural ectoderm",
  "Paraxial mesoderm",
  "Epidermis",
  "Lateral plate mesoderm"
)

# Define cell type merging rules
cell_type_merging <- list(
  "Cardiopharyngeal mesoderm" = c(
    "Anterior cardiopharyngeal progenitors",
    "Cardiopharyngeal progenitors",
    "Cardiopharyngeal progenitors FHF",
    "Cardiopharyngeal progenitors SHF",
    "Cardiomyocytes FHF 1",
    "Cardiomyocytes FHF",
    "Cardiomyocytes SHF 1",
    "Cardiomyocytes SHF",
    "Pharyngeal mesoderm",
    "Cardiomyocytes"
  ),
  "Foregut" = c(
    "Foregut",
    "Pharyngeal endoderm",
    "Thyroid primordium"
  )
)

# Define genes for PGC to Naive pluripotency conversion
naive_pluripotency_genes <- c("Nanog", "Zfp42", "Klf2", "Klf4")

# Function to check if a cell expresses naive pluripotency genes
check_naive_pluripotency <- function(seurat_obj, cell_indices) {
  # Get expression data for the genes
  gene_expr <- FetchData(seurat_obj, vars = naive_pluripotency_genes, cells = cell_indices)
  
  # Check if all 4 genes are expressed (count > 0)
  expressed_genes <- rowSums(gene_expr > 0)
  return(expressed_genes == 4)
}


# Function to process each time point
process_time_point <- function(time_point) {
  cat("Processing time point:", time_point, "\n")
  
  # Load the Marioni-labeled dataset
  input_file <- paste0("401_integrate_by_time/", time_point, "/401_integrated_", time_point, ".rds")
  seurat_obj <- readRDS(input_file)
  
  cat("Original dataset dimensions:", dim(seurat_obj), "\n")
  cat("Original cell types:", length(unique(seurat_obj$sc_merged_gast_transfer_label)), "\n")
  
  # Create a copy of the original sc_merged_gast_transfer_label for reference
  seurat_obj$original_predicted.id <- seurat_obj$sc_merged_gast_transfer_label
  
  # Create manual labels based on the rules
  manual_labels <- seurat_obj$sc_merged_gast_transfer_label
  
  # Apply merging rules
  for (new_label in names(cell_type_merging)) {
    old_labels <- cell_type_merging[[new_label]]
    manual_labels[manual_labels %in% old_labels] <- new_label
  }
  
  # Special handling for PGC -> Naive pluripotency conversion based on gene expression
  pgc_cells <- which(seurat_obj$sc_merged_gast_transfer_label == "PGC")
  if (length(pgc_cells) > 0) {
    cat("Checking", length(pgc_cells), "PGC cells for naive pluripotency gene expression...\n")
    
    # Check gene expression for PGC cells
    naive_expr <- check_naive_pluripotency(seurat_obj, pgc_cells)
    
    # Convert PGC to Naive pluripotency if genes are expressed
    pgc_to_convert <- pgc_cells[naive_expr]
    manual_labels[pgc_to_convert] <- "Naive pluripotency"
    
    cat("Converted", length(pgc_to_convert), "PGC cells to Naive pluripotency based on gene expression\n")
    cat("Kept", length(pgc_cells) - length(pgc_to_convert), "PGC cells as PGC\n")
  }
  
  # Add manual labels to metadata
  seurat_obj$manual_marioni_label <- manual_labels
  
  # Create filtering vector (TRUE = keep cell, FALSE = remove cell)
  keep_cells <- !(manual_labels %in% cell_types_to_remove)
  
  # Add filtering info to metadata
  seurat_obj$keep_cell <- keep_cells
  seurat_obj$removed_reason <- ifelse(keep_cells, "Kept", "Removed")
  
  # Count cells before filtering
  cells_before <- ncol(seurat_obj)
  removed_cells <- sum(!keep_cells)
  
  cat("Cells before filtering:", cells_before, "\n")
  cat("Cells to be removed:", removed_cells, "\n")
  cat("Cells after filtering:", cells_before - removed_cells, "\n")
  
  # Create summary of removed cell types
  removed_summary <- table(manual_labels[!keep_cells])
  cat("Removed cell types:\n")
  print(removed_summary)
  
  # Create summary of merged cell types
  merged_summary <- list()
  for (new_label in names(cell_type_merging)) {
    old_labels <- cell_type_merging[[new_label]]
    merged_cells <- sum(seurat_obj$sc_merged_gast_transfer_label %in% old_labels)
    merged_summary[[new_label]] <- merged_cells
  }
  
  # Add PGC to Naive pluripotency conversion summary
  if (length(pgc_cells) > 0) {
    pgc_converted <- sum(manual_labels[pgc_cells] == "Naive pluripotency")
    pgc_kept <- sum(manual_labels[pgc_cells] == "PGC")
    merged_summary[["PGC -> Naive pluripotency (gene-based)"]] <- pgc_converted
    merged_summary[["PGC (kept as PGC)"]] <- pgc_kept
  }
  
  cat("Merged cell types:\n")
  print(merged_summary)
  
  # Create visualization before filtering
  pdf(file.path(output_dir, paste0("601_", time_point, "_before_filtering.pdf")), 
      height = 20, width = 25)
  
  # Plot 1: Original Marioni predictions
  p1 <- DimPlot(seurat_obj, group.by = "original_predicted.id", label = TRUE, repel = TRUE) +
    ggtitle(paste0(toupper(time_point), " - Original Marioni Predictions")) +
    theme(legend.position = "none")
  
  # Plot 2: Manual labels (after merging)
  p2 <- DimPlot(seurat_obj, group.by = "manual_marioni_label", label = TRUE, repel = TRUE) +
    ggtitle(paste0(toupper(time_point), " - Manual Labels (After Merging)")) +
    theme(legend.position = "none")
  
  # Plot 3: Cells to keep vs remove
  p3 <- DimPlot(seurat_obj, group.by = "removed_reason", cols = c("Kept" = "blue", "Removed" = "red")) +
    ggtitle(paste0(toupper(time_point), " - Cells to Keep (Blue) vs Remove (Red)"))
  
  # Plot 4: Confidence scores
  p4 <- FeaturePlot(seurat_obj, features = "prediction.score.max") +
    ggtitle(paste0(toupper(time_point), " - Prediction Confidence Scores"))
  
  # Plot 5: Naive pluripotency gene expression (if PGC cells exist)
  if (length(pgc_cells) > 0) {
    p5 <- FeaturePlot(seurat_obj, features = naive_pluripotency_genes, ncol = 2) +
      ggtitle(paste0(toupper(time_point), " - Naive Pluripotency Gene Expression"))
    print(p5)
  }
  
  # Plot 6: Cell type composition before filtering
  cell_type_counts_before <- table(manual_labels)
  cell_type_counts_df_before <- data.frame(
    CellType = names(cell_type_counts_before),
    Count = as.numeric(cell_type_counts_before),
    Percentage = round(as.numeric(cell_type_counts_before) / sum(cell_type_counts_before) * 100, 2)
  ) %>%
    arrange(desc(Count))
  
  p6 <- ggplot(cell_type_counts_df_before, aes(x = reorder(CellType, Count), y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    geom_text(aes(label = paste0(Count, " (", Percentage, "%)")), 
              hjust = -0.1, size = 3) +
    coord_flip() +
    labs(title = paste0(toupper(time_point), " - Cell Type Composition Before Filtering"),
         x = "Cell Type",
         y = "Number of Cells") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 14, face = "bold"))
  
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p6)
  
  dev.off()
  
  # Create PNG versions for presentation
  png(file.path(output_dir, paste0("601_", time_point, "_before_filtering.png")), 
      height = 20*300, width = 25*300, res = 300)
  
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  print(p6)
  
  dev.off()
  
  # Filter the dataset
  seurat_obj_filtered <- subset(seurat_obj, subset = keep_cell == TRUE)
  
  cat("Filtered dataset dimensions:", dim(seurat_obj_filtered), "\n")
  cat("Filtered cell types:", length(unique(seurat_obj_filtered$manual_marioni_label)), "\n")
  
  # Create visualization after filtering
  pdf(file.path(output_dir, paste0("601_", time_point, "_after_filtering.pdf")), 
      height = 18, width = 25)
  
  # Plot 1: Filtered manual labels
  p1 <- DimPlot(seurat_obj_filtered, group.by = "manual_marioni_label", label = TRUE, repel = TRUE) +
    ggtitle(paste0(toupper(time_point), " - Filtered Manual Labels")) +
    theme(legend.position = "none")
  
  # Plot 2: Original Seurat clusters
  p2 <- DimPlot(seurat_obj_filtered, group.by = "seurat_clusters", label = TRUE) +
    ggtitle(paste0(toupper(time_point), " - Original Seurat Clusters"))
  
  # Plot 3: Confidence scores (filtered)
  p3 <- FeaturePlot(seurat_obj_filtered, features = "prediction.score.max") +
    ggtitle(paste0(toupper(time_point), " - Prediction Confidence Scores (Filtered)"))
  
  # Plot 4: Cell type composition after filtering
  cell_type_counts_after <- table(seurat_obj_filtered$manual_marioni_label)
  cell_type_counts_df_after <- data.frame(
    CellType = names(cell_type_counts_after),
    Count = as.numeric(cell_type_counts_after),
    Percentage = round(as.numeric(cell_type_counts_after) / sum(cell_type_counts_after) * 100, 2)
  ) %>%
    arrange(desc(Count))
  
  p4 <- ggplot(cell_type_counts_df_after, aes(x = reorder(CellType, Count), y = Count)) +
    geom_bar(stat = "identity", fill = "darkgreen", alpha = 0.8) +
    geom_text(aes(label = paste0(Count, " (", Percentage, "%)")), 
              hjust = -0.1, size = 3) +
    coord_flip() +
    labs(title = paste0(toupper(time_point), " - Cell Type Composition After Filtering"),
         x = "Cell Type",
         y = "Number of Cells") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 14, face = "bold"))
  
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  
  dev.off()
  
  # Create PNG versions for presentation
  png(file.path(output_dir, paste0("601_", time_point, "_after_filtering.png")), 
      height = 18*300, width = 25*300, res = 300)
  
  print(p1)
  print(p2)
  print(p3)
  print(p4)
  
  dev.off()
  
  # Create summary statistics
  summary_stats <- data.frame(
    TimePoint = time_point,
    Cells_Before = cells_before,
    Cells_Removed = removed_cells,
    Cells_After = cells_before - removed_cells,
    Cell_Types_Before = length(unique(seurat_obj$sc_merged_gast_transfer_label)),
    Cell_Types_After = length(unique(seurat_obj_filtered$manual_marioni_label))
  )
  
  # Create cell type composition table
  cell_type_composition <- table(seurat_obj_filtered$manual_marioni_label)
  cell_type_composition_df <- data.frame(
    CellType = names(cell_type_composition),
    Count = as.numeric(cell_type_composition),
    Percentage = round(as.numeric(cell_type_composition) / sum(cell_type_composition) * 100, 2)
  ) %>%
    arrange(desc(Count))
  
  # Create barplot of cell type composition
  pdf(file.path(output_dir, paste0("601_", time_point, "_cell_type_composition.pdf")), 
      height = 10, width = 15)
  
  p <- ggplot(cell_type_composition_df, aes(x = reorder(CellType, Count), y = Count)) +
    geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
    geom_text(aes(label = paste0(Count, " (", Percentage, "%)")), 
              hjust = -0.1, size = 3) +
    coord_flip() +
    labs(title = paste0(toupper(time_point), " - Cell Type Composition After Filtering"),
         x = "Cell Type",
         y = "Number of Cells") +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10),
          plot.title = element_text(size = 14, face = "bold"))
  
  print(p)
  dev.off()
  
  # Create PNG version for presentation
  png(file.path(output_dir, paste0("601_", time_point, "_cell_type_composition.png")), 
      height = 10*300, width = 15*300, res = 300)
  print(p)
  dev.off()
  
  # Save summary statistics
  write.csv(summary_stats, file.path(output_dir, paste0("601_", time_point, "_summary_stats.csv")), 
            row.names = FALSE)
  write.csv(cell_type_composition_df, file.path(output_dir, paste0("601_", time_point, "_cell_type_composition.csv")), 
            row.names = FALSE)
  
  # Save the filtered dataset
  saveRDS(seurat_obj_filtered, 
          file.path(output_dir, paste0("601_", time_point, "_filtered.rds")))
  
  
  cat("Completed processing for time point:", time_point, "\n")
  return(list(
    original = seurat_obj,
    filtered = seurat_obj_filtered,
    summary = summary_stats,
    composition = cell_type_composition_df
  ))
}

# Process all time points
results <- list()
for (time_point in time_points) {
  tryCatch({
    results[[time_point]] <- process_time_point(time_point)
  }, error = function(e) {
    cat("Error processing time point", time_point, ":", e$message, "\n")
  })
}

# Create overall summary
overall_summary <- do.call(rbind, lapply(results, function(x) x$summary))
write.csv(overall_summary, file.path(output_dir, "601_overall_summary_time_based.csv"), row.names = FALSE)

# Create combined cell type composition
all_compositions <- do.call(rbind, lapply(names(results), function(time_point) {
  comp <- results[[time_point]]$composition
  comp$TimePoint <- time_point
  return(comp)
}))

write.csv(all_compositions, file.path(output_dir, "601_all_timepoints_cell_type_composition.csv"), 
          row.names = FALSE)

# Create detailed merging analysis
cat("\n=== DETAILED MERGING ANALYSIS ===\n")

# Analyze Foregut merging
foregut_analysis <- all_compositions %>%
  filter(CellType %in% c("Foregut", "Pharingeal endoderm", "Thyroid primordium")) %>%
  group_by(TimePoint) %>%
  summarise(
    Foregut_merged = sum(Count[CellType == "Foregut"]),
    Pharingeal_endoderm = sum(Count[CellType == "Pharingeal endoderm"]),
    Thyroid_primordium = sum(Count[CellType == "Thyroid_primordium"]),
    Total_foregut_cells = sum(Count),
    .groups = "drop"
  )

# Analyze Cardiopharingeal merging (look for individual components)
cardiopharingeal_analysis <- all_compositions %>%
  filter(grepl("cardiophar|cardiomyo|pharyngeal", CellType, ignore.case = TRUE)) %>%
  group_by(TimePoint) %>%
  summarise(
    Total_cardiopharingeal_related = sum(Count),
    Cell_types = paste(unique(CellType), collapse = ", "),
    .groups = "drop"
  )

# Analyze PGC to Naive pluripotency conversion
pgc_analysis <- all_compositions %>%
  filter(CellType %in% c("PGC", "Naive pluripotency")) %>%
  group_by(TimePoint) %>%
  summarise(
    PGC_count = sum(Count[CellType == "PGC"]),
    Naive_pluripotency_count = sum(Count[CellType == "Naive pluripotency"]),
    Total_pgc_related = sum(Count),
    .groups = "drop"
  )

# Print detailed results
cat("\n=== FOREGUT MERGING ANALYSIS ===\n")
print(foregut_analysis)

cat("\n=== CARDIOPHARINGEAL MERGING ANALYSIS ===\n")
print(cardiopharingeal_analysis)

cat("\n=== PGC TO NAIVE PLURIPOTENCY CONVERSION ===\n")
print(pgc_analysis)

# Create summary
cat("\n=== MERGING SUMMARY ===\n")
cat("Foregut merging:\n")
for(i in 1:nrow(foregut_analysis)) {
  time_pt <- foregut_analysis$TimePoint[i]
  total <- foregut_analysis$Total_foregut_cells[i]
  cat(sprintf("- %s: %d cells merged into Foregut\n", toupper(time_pt), total))
}

cat("\nPGC to Naive pluripotency conversion:\n")
for(i in 1:nrow(pgc_analysis)) {
  time_pt <- pgc_analysis$TimePoint[i]
  converted <- pgc_analysis$Naive_pluripotency_count[i]
  pgc_kept <- pgc_analysis$PGC_count[i]
  conversion_rate <- round(converted / (converted + pgc_kept) * 100, 1)
  cat(sprintf("- %s: %d converted (%.1f%%), %d kept as PGC\n", 
              toupper(time_pt), converted, conversion_rate, pgc_kept))
}

# Save detailed merging analysis to CSV files
write.csv(foregut_analysis, file.path(output_dir, "601_foregut_merging_analysis_time_based.csv"), row.names = FALSE)
write.csv(cardiopharingeal_analysis, file.path(output_dir, "601_cardiopharingeal_analysis_time_based.csv"), row.names = FALSE)
write.csv(pgc_analysis, file.path(output_dir, "601_pgc_conversion_analysis_time_based.csv"), row.names = FALSE)

# Create combined visualization
pdf(file.path(output_dir, "601_all_timepoints_comparison.pdf"), height = 20, width = 25)

# Plot 1: Cell counts comparison
p1 <- ggplot(overall_summary, aes(x = TimePoint, y = Cells_After, fill = TimePoint)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_text(aes(label = Cells_After), vjust = -0.5, size = 4) +
  labs(title = "Cell Counts After Filtering by Time Point",
       x = "Time Point",
       y = "Number of Cells") +
  theme_minimal() +
  theme(legend.position = "none")

# Plot 2: Cell type composition by time point
p2 <- ggplot(all_compositions, aes(x = TimePoint, y = Count, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Cell Type Composition by Time Point",
       x = "Time Point",
       y = "Number of Cells",
       fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot 3: Cell type percentage by time point
p3 <- ggplot(all_compositions, aes(x = TimePoint, y = Percentage, fill = CellType)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Cell Type Percentage by Time Point",
       x = "Time Point",
       y = "Percentage (%)",
       fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p1)
print(p2)
print(p3)

dev.off()

# Create PNG version for presentation
png(file.path(output_dir, "601_all_timepoints_comparison.png"), 
    height = 20*300, width = 25*300, res = 300)

print(p1)
print(p2)
print(p3)

dev.off()

cat("Processing completed for all time points.\n")
cat("Output files saved in:", output_dir, "\n")
cat("Processed time points:", paste(time_points, collapse = ", "), "\n")
cat("Summary statistics saved in: 601_overall_summary_time_based.csv\n")
cat("Cell type compositions saved in: 601_all_timepoints_cell_type_composition.csv\n")
cat("Detailed merging analysis saved in:\n")
cat("  - 601_foregut_merging_analysis_time_based.csv\n")
cat("  - 601_cardiopharingeal_analysis_time_based.csv\n")
cat("  - 601_pgc_conversion_analysis_time_based.csv\n")
cat("\nPNG files for presentation created:\n")
cat("  - 601_{timepoint}_before_filtering.png\n")
cat("  - 601_{timepoint}_after_filtering.png\n")
cat("  - 601_{timepoint}_cell_type_composition.png\n")
cat("  - 601_all_timepoints_comparison.png\n")
cat("\nAll PNG files are high-resolution (300 DPI) for presentation quality.\n")
