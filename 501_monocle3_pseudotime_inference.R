library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(Polychrome)
library(tidyverse)
library(Seurat)
library(viridis)

# Memory management function
get_memory_usage <- function() {
  mem <- gc()
  cat("Memory usage - Vcells:", round(mem[1,2]/1024^3, 2), "GB, Ncells:", round(mem[2,2]/1024^3, 2), "GB\n")
}

# Clean up function
cleanup_memory <- function() {
  gc()
  invisible()
}

# Define all conditions to process
conditions <- c("ctrl", "nacl", "racl", "xal")

# Create output directory for plots
output_dir <- "501_monocle3_pseudotime_output"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

# Function to process each condition
process_condition <- function(condition) {
  cat("Processing condition:", condition, "\n")
  get_memory_usage()
  
  # Define file paths
  marioni_file <- paste0("103_marioni_on_condition_normal_integrated/103_", condition, "_marioni_gast_transfer.rds")
  aizarani_file <- paste0("301_aizarani_on_condition_normal_integrated/301_", condition, "_aizarani_gast_transfer.rds")
  
  # Check if files exist
  if (!file.exists(marioni_file)) {
    cat("Warning: Marioni file not found:", marioni_file, "\n")
    return(NULL)
  }
  if (!file.exists(aizarani_file)) {
    cat("Warning: Aizarani file not found:", aizarani_file, "\n")
    return(NULL)
  }
  
  tryCatch({
    # Load Seurat objects one at a time to manage memory
    cat("Loading Marioni object...\n")
    seurat_obj_tpm <- readRDS(marioni_file)
    get_memory_usage()
    
    # Extract only the metadata we need to save memory
    cat("Extracting Marioni metadata...\n")
    marioni_meta <- seurat_obj_tpm@meta.data %>% select(predicted.id)
    colnames(marioni_meta) <- c("predicted_id_marioni")
    
    # Remove the large Marioni object immediately
    rm(seurat_obj_tpm)
    cleanup_memory()
    cat("Marioni object removed from memory\n")
    get_memory_usage()
    
    cat("Loading Aizarani object...\n")
    seurat_obj <- readRDS(aizarani_file)
    get_memory_usage()
    
    # Set active assay
    seurat_obj@active.assay <- 'RNA'
    
    # Add Marioni predictions to Aizarani object
    cat("Adding Marioni predictions...\n")
    seurat_obj <- AddMetaData(seurat_obj, marioni_meta)
    
    # Clean up metadata object
    rm(marioni_meta)
    cleanup_memory()
    
    # Select root cell for pseudotime
    cat("Selecting root cell...\n")
    # Genes that are expressed in Naive pluripotency cells
    # Genes that should not be expressed in Naive pluripotency cells
    # Subset the temporary Seurat object to select candidate root cells based on gene expression.
    # For non-control conditions, select cells expressing Nanog, Zfp42, Klf2, Klf4,
    # and not expressing Dazl, Tfap2c, Fgf5, or Gata2.
    if (condition != "ctrl") {
      tmp_seurat_obj <- subset(
        seurat_obj,
        subset = Nanog > 0 & Zfp42 > 0 & Klf2 > 0 & Klf4 > 0 &
                 Dazl == 0 & Tfap2c == 0 & Fgf5 == 0 & Gata2 == 0
      )
    }
    # For the control condition, select cells expressing Nanog, Zfp42, Klf2, Klf4,
    # and not expressing Fgf5 or Gata2 (manual workaround for control).
    else if (condition == "ctrl") {
      tmp_seurat_obj <- subset(
        seurat_obj,
        subset = Nanog > 0 & Zfp42 > 0 & Klf2 > 0 & Klf4 > 0 & Gata2 == 0
      )
    }
    # Filter the metadata of the temporary Seurat object to select cells at 48h timepoint
    # and with a high prediction score (> 0.7) for the root cell selection.
    selected_cell_tmp <- tmp_seurat_obj@meta.data %>%
      filter(time == "48h") %>%
      filter(prediction.score.max > 0.7)
    
    # Select the first cell that matches the criteria as the root cell.
    selected_cell <- rownames(selected_cell_tmp)[1]
    
    # Remove the temporary Seurat object from memory to free up resources.
    rm(tmp_seurat_obj)
    cleanup_memory()
    
    if (is.na(selected_cell) || length(selected_cell) == 0) {
      cat("Warning: No suitable root cell found for condition", condition, "\n")
      cleanup_memory()
      return(NULL)
    }
    
    cat("Selected root cell:", selected_cell, "\n")
    
    # Set active assay
    seurat_obj@active.assay <- 'integrated'
    
    # Convert to Cell Data Set
    cat("Converting to Cell Data Set...\n")
    cds <- as.cell_data_set(seurat_obj)
    get_memory_usage()
    
    # Remove the large Seurat object after conversion
    rm(seurat_obj)
    cleanup_memory()
    cat("Seurat object removed from memory\n")
    get_memory_usage()
    
    # Process CDS step by step with memory cleanup
    cat("Clustering cells...\n")
    cds <- cluster_cells(cds = cds)
    cleanup_memory()
    
    cat("Learning graph...\n")
    cds <- learn_graph(cds, use_partition = FALSE)
    cleanup_memory()
    
    cat("Ordering cells...\n")
    cds <- order_cells(cds, root_cells = selected_cell)
    cleanup_memory()
    
    # Generate plots with memory management
    cat("Generating plots...\n")
    
    # Plot 1: Pseudotime
    p1 <- plot_cells(cds, color_cells_by = "pseudotime", 
                     label_branch_points=FALSE, label_leaves=FALSE)
    ggsave(paste0(output_dir, "/", condition, "_pseudotime.pdf"), p1, width = 10, height = 8)
    rm(p1)
    cleanup_memory()
    
    # Plot 2: Time
    p2 <- plot_cells(cds, color_cells_by = "time", 
                     label_branch_points=FALSE, label_leaves=FALSE)
    ggsave(paste0(output_dir, "/", condition, "_time.pdf"), p2, width = 10, height = 8)
    rm(p2)
    cleanup_memory()
    
    # Plot 3: Marioni predictions
    p3 <- plot_cells(cds, color_cells_by = "predicted_id_marioni")
    ggsave(paste0(output_dir, "/", condition, "_marioni_predictions.pdf"), p3, width = 10, height = 8)
    rm(p3)
    cleanup_memory()
    
    # Plot 4: Marioni predictions with labels
    p4 <- plot_cells(cds, group_label_size = 4, label_cell_groups = TRUE, 
                     color_cells_by = "predicted_id_marioni")
    ggsave(paste0(output_dir, "/", condition, "_marioni_predictions_labeled.pdf"), p4, width = 12, height = 10)
    rm(p4)
    cleanup_memory()
    
    # Add pseudotime to metadata
    cat("Calculating pseudotime...\n")
    colData(cds)$pseudotime <- pseudotime(cds)
    cleanup_memory()
    
    # Plot 5: Pseudotime with trajectory graph
    p5 <- plot_cells(
      cds,
      color_cells_by = "pseudotime",
      show_trajectory_graph = TRUE,
      label_cell_groups = FALSE,
      label_leaves = FALSE,
      label_branch_points = FALSE
    )
    ggsave(paste0(output_dir, "/", condition, "_pseudotime_trajectory.pdf"), p5, width = 10, height = 8)
    rm(p5)
    cleanup_memory()
    
    # Plot 6: Faceted by time with viridis color
    p6 <- plot_cells(
      cds,
      color_cells_by = "pseudotime",
      show_trajectory_graph = TRUE,
      label_cell_groups = FALSE,
      label_leaves = FALSE,
      label_branch_points = FALSE
    ) + 
      facet_wrap(~time) +
      scale_color_gradientn(colors = viridis::viridis(100)) +
      theme_minimal() +
      theme(strip.text = element_text(size = 12))
    ggsave(paste0(output_dir, "/", condition, "_pseudotime_faceted_time.pdf"), p6, width = 15, height = 10)
    rm(p6)
    cleanup_memory()
    
    # Plot 7: Boxplot of pseudotime by time and cluster
    p7 <- as.data.frame(colData(cds)) %>%
      ggplot(aes(x = time, y = pseudotime, fill = seurat_clusters)) +
      geom_boxplot(outlier.size = 0.3) +
      theme_minimal() +
      labs(y = "Pseudotime", fill = "Cluster", title = paste("Condition:", condition))
    ggsave(paste0(output_dir, "/", condition, "_pseudotime_boxplot.pdf"), p7, width = 12, height = 8)
    rm(p7)
    cleanup_memory()
    
    # Save the processed CDS object
    cat("Saving CDS object...\n")
    saveRDS(cds, paste0(output_dir, "/", condition, "_monocle3_cds.rds"))

    # Note: seurat_obj was removed from memory earlier to save space
    # The pseudotime information is available in the CDS object

    
    cat("Successfully processed condition:", condition, "\n")
    get_memory_usage()
    return(cds)
    
  }, error = function(e) {
    cat("Error processing condition", condition, ":", e$message, "\n")
    # Clean up any remaining objects on error
    if (exists("seurat_obj_tpm")) rm(seurat_obj_tpm)
    if (exists("seurat_obj")) rm(seurat_obj)
    if (exists("marioni_meta")) rm(marioni_meta)
    if (exists("cds")) rm(cds)
    cleanup_memory()
    return(NULL)
  })
}

# Process all conditions
cat("Starting Monocle3 pseudotime analysis for all conditions...\n")
cat("Conditions to process:", paste(conditions, collapse = ", "), "\n\n")
get_memory_usage()

# Process each condition
results <- list()
for (condition in conditions) {
  cat(paste(rep("=", 50), collapse = ""), "\n")
  cat("Processing condition:", condition, "\n")
  get_memory_usage()
  
  result <- process_condition(condition)
  if (!is.null(result)) {
    results[[condition]] <- result
  }
  
  # Force cleanup between conditions
  cleanup_memory()
  cat("Memory cleanup completed for condition:", condition, "\n")
  get_memory_usage()
  cat(paste(rep("=", 50), collapse = ""), "\n\n")
}

# Summary
cat("Analysis complete!\n")
cat("Successfully processed conditions:", names(results), "\n")
cat("Failed conditions:", setdiff(conditions, names(results)), "\n")
cat("Output saved to:", output_dir, "\n")

# Optional: Create a combined summary plot if multiple conditions succeeded
if (length(results) > 1) {
  cat("Creating combined summary...\n")
  get_memory_usage()
  
  # Combine pseudotime data from all conditions with memory management
  combined_data <- do.call(rbind, lapply(names(results), function(cond) {
    cds <- results[[cond]]
    data <- data.frame(
      condition = cond,
      pseudotime = pseudotime(cds),
      time = colData(cds)$time,
      cluster = colData(cds)$seurat_clusters
    )
    # Clean up intermediate data
    cleanup_memory()
    return(data)
  }))
  
  # Combined boxplot
  p_combined <- combined_data %>%
    ggplot(aes(x = time, y = pseudotime, fill = condition)) +
    geom_boxplot(outlier.size = 0.3) +
    theme_minimal() +
    labs(y = "Pseudotime", fill = "Condition", title = "Pseudotime across all conditions") +
    facet_wrap(~condition, scales = "free_x")
  
  ggsave(paste0(output_dir, "/combined_pseudotime_all_conditions.pdf"), p_combined, width = 16, height = 10)
  rm(p_combined, combined_data)
  cleanup_memory()
  
  cat("Combined summary plot saved!\n")
}

# Final cleanup
cat("Final memory cleanup...\n")
if (length(results) > 0) {
  rm(results)
}
cleanup_memory()
get_memory_usage()
cat("Script completed successfully!\n")

# Post run comments:
# Root cells selected based on gene expression:
# ctrl: Ctrl48h_CTTCAATAGGCTGGTAACTTTAGG-1
# nacl: NACL48h_TGAGGTCCAGCACTAAAGTAGGCT-1
# racl: RACL48h_CTTCGAATCGTTGCACAACGGGAA-1
# xal: XAL48h_GGTTTGACAGCCAGCAATGTTGAC-1