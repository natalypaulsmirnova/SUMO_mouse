# =============================================================================
# CellChat Analysis Script - All Conditions (Generalized)
# =============================================================================
# This script performs cell-cell communication analysis using CellChat
# on all conditions (Ctrl, NACL, RACL, XAL) with user-specified cell type subsets
#
# Usage: Rscript 701_cellchat_all_condition_cellchat_generalized.R --celltypes "Epiblast,Primitive Streak,Anterior Primitive Streak" --suffix "subset1"
# Arguments:
#   --celltypes: Comma-separated list of cell types to analyze (required)
#   --suffix: Suffix to add to output directory name (optional, default: "generalized")

# Load required libraries
library(Seurat)
library(CellChat)
library(dplyr)
library(ggplot2)
library(optparse)

# =============================================================================
# Command Line Argument Parsing
# =============================================================================

# Define command line options
option_list <- list(
  make_option(c("--celltypes"), type="character", default=NULL,
              help="Comma-separated list of cell types to analyze (required)"),
  make_option(c("--suffix"), type="character", default="generalized",
              help="Suffix to add to output directory name (default: generalized)")
)

# Parse command line arguments
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$celltypes)) {
  print_help(opt_parser)
  stop("--celltypes argument is required. Please specify cell types as comma-separated values.")
}

# Parse cell types from command line argument
celltypes_to_keep <- trimws(unlist(strsplit(opt$celltypes, ",")))

# =============================================================================
# Configuration and Setup
# =============================================================================

# Define all conditions to process
conditions <- c("ctrl", "nacl", "racl", "xal")

# Create output directory based on script basename and suffix
script_name <- tools::file_path_sans_ext(basename("701_cellchat_all_condition_cellchat_generalized.R"))
output_dir <- paste0(script_name, "_", opt$suffix)

# Check if output directory exists and show warning if it does
if (dir.exists(output_dir)) {
  warning("Output directory '", output_dir, "' already exists. Files may be overwritten.")
} else {
  dir.create(output_dir, showWarnings = FALSE)
}

# Analysis parameters
celltype_col <- "manual_marioni_label"  # Column name containing cell type annotations
assay_to_use <- "RNA"                   # Assay to use for expression data
min_cells_per_group <- 10               # Minimum cells per group to stabilize probabilities

# Log the configuration
cat("Configuration:\n")
cat("  Cell types to analyze:", paste(celltypes_to_keep, collapse = ", "), "\n")
cat("  Output directory:", output_dir, "\n")
cat("  Conditions to process:", paste(conditions, collapse = ", "), "\n")
cat("\n")

# =============================================================================
# Helper Functions
# =============================================================================

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

# Function to perform CellChat analysis for a single condition
perform_cellchat_analysis <- function(condition) {
  cat("Processing condition:", condition, "\n")
  
  # Define input file path
  input_file <- paste0("601_manual_labeling_marioni_filtering_time_based/601_", condition, "_filtered.rds")
  
  # Check if input file exists
  if (!file.exists(input_file)) {
    cat("  ERROR: Input file not found:", input_file, "\n")
    return(NULL)
  }
  
  tryCatch({
    # Load the Seurat object
    cat("  Loading Seurat object...\n")
    obj <- readRDS(input_file)
    
    # Join layers if the object has multiple layers (e.g., from SCTransform)
    obj <- JoinLayers(obj)
    
    # Set cell type labels as the active identity in the Seurat object
    if (!is.null(celltype_col)) {
      stopifnot(celltype_col %in% colnames(obj@meta.data))
      Idents(obj) <- obj[[celltype_col, drop = TRUE]]
    }
    
    # Check cell type availability and get present types
    present_celltypes <- check_cell_types(obj, condition, celltypes_to_keep)
    
    if (length(present_celltypes) == 0) {
      cat("  ERROR: No target cell types found in", condition, "\n")
      return(NULL)
    }
    
    # Subset to only include the present cell types
    cat("  Subsetting to", length(present_celltypes), "cell types...\n")
    obj_sub <- subset(obj, idents = present_celltypes)
    
    # Remove the original object to free up memory
    rm(obj)
    gc()
    
    # Remove unused factor levels from the cell type identities
    Idents(obj_sub) <- droplevels(Idents(obj_sub))
    
    # Further filter to remove cell types with fewer than min_cells_per_group cells
    obj_sub <- subset(obj_sub, idents = names(which(table(Idents(obj_sub)) >= min_cells_per_group)))
    
    # Check if we still have enough cell types after filtering
    if (length(levels(Idents(obj_sub))) < 2) {
      cat("  ERROR: Not enough cell types remaining after filtering in", condition, "\n")
      return(NULL)
    }
    
    cat("  Final cell types:", length(levels(Idents(obj_sub))), "\n")
    cat("  Final cell types:", paste(levels(Idents(obj_sub)), collapse = ", "), "\n")
    
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
    
    cat("  Creating CellChat object...\n")
    # Create CellChat object with expression data and metadata
    cellchat <- createCellChat(object = expr, meta = meta, group.by = "labels")
    
    # Set the database to mouse (CellChatDB.mouse contains mouse ligand-receptor pairs)
    cellchat@DB <- CellChatDB.mouse
    
    cat("  Running CellChat analysis pipeline...\n")
    # Subset the database to only include genes expressed in the dataset
    cellchat <- subsetData(cellchat)
    
    # Identify over-expressed genes in each cell group
    cellchat <- identifyOverExpressedGenes(cellchat)
    
    # Identify over-expressed ligand-receptor interactions
    cellchat <- identifyOverExpressedInteractions(cellchat)
    
    # Compute communication probabilities using triMean method
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
    # Create Output Directory and Save Results
    # =============================================================================
    
    # Create condition-specific output directory
    condition_output_dir <- file.path(output_dir, condition)
    
    if (dir.exists(condition_output_dir)) {
      warning("Condition directory '", condition_output_dir, "' already exists. Files may be overwritten.")
    } else {
      dir.create(condition_output_dir, showWarnings = FALSE)
    }
    
    # Save CellChat object
    saveRDS(cellchat, file.path(condition_output_dir, paste0("701_", condition, "_", opt$suffix, "_cellchat.rds")))
    
    # =============================================================================
    # Generate Visualizations
    # =============================================================================
    
    cat("  Generating visualizations...\n")
    
    # Calculate group sizes for visualization
    group_size <- as.numeric(table(cellchat@idents))
    
    # Create circular network plot showing interaction strength
    pdf(file.path(condition_output_dir, paste0("701_", condition, "_", opt$suffix, "_interaction_strength.pdf")), 
        width = 10, height = 10)
    netVisual_circle(cellchat@net$weight, vertex.weight = group_size, 
                    title.name = paste0("Interaction Strength - ", toupper(condition), " (", opt$suffix, ")"))
    dev.off()
    
    # =============================================================================
    # Generate Summary Statistics
    # =============================================================================
    
    # Extract communication information
    pair_pathways <- subsetCommunication(cellchat)
    
    # Save communication data
    write.csv(pair_pathways, 
              file.path(condition_output_dir, paste0("701_", condition, "_", opt$suffix, "_communication_pairs.csv")),
              row.names = FALSE)
    
    # Create summary statistics
    summary_stats <- data.frame(
      condition = condition,
      suffix = opt$suffix,
      n_cell_types = length(levels(cellchat@idents)),
      cell_types = paste(levels(cellchat@idents), collapse = ";"),
      n_cells = sum(group_size),
      n_interactions = nrow(pair_pathways),
      n_pathways = length(unique(pair_pathways$pathway_name)),
      stringsAsFactors = FALSE
    )
    
    # Save summary statistics
    write.csv(summary_stats, 
              file.path(condition_output_dir, paste0("701_", condition, "_", opt$suffix, "_summary_stats.csv")),
              row.names = FALSE)
    
    cat("  Analysis completed successfully for", condition, "\n")
    cat("  - Cell types:", length(levels(cellchat@idents)), "\n")
    cat("  - Total cells:", sum(group_size), "\n")
    cat("  - Interactions found:", nrow(pair_pathways), "\n")
    cat("  - Pathways identified:", length(unique(pair_pathways$pathway_name)), "\n")
    
    return(summary_stats)
    
  }, error = function(e) {
    cat("  ERROR in", condition, ":", e$message, "\n")
    return(NULL)
  })
}

# =============================================================================
# Main Analysis Loop
# =============================================================================

# Initialize summary data storage
all_summary_stats <- list()

# Process each condition
for (condition in conditions) {
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("Processing condition:", toupper(condition), "\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  summary_stats <- perform_cellchat_analysis(condition)
  
  if (!is.null(summary_stats)) {
    all_summary_stats[[condition]] <- summary_stats
  }
  
  # Clean up memory
  gc()
}

# =============================================================================
# Generate Overall Summary
# =============================================================================

if (length(all_summary_stats) > 0) {
  # Combine all summary statistics
  overall_summary <- do.call(rbind, all_summary_stats)
  
  # Save overall summary
  write.csv(overall_summary, 
            file.path(output_dir, paste0("701_all_conditions_", opt$suffix, "_summary.csv")),
            row.names = FALSE)
  
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("ANALYSIS COMPLETE\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  cat("Successfully processed", nrow(overall_summary), "conditions\n")
  cat("Output directory:", output_dir, "\n")
  cat("\nSummary by condition:\n")
  print(overall_summary)
  
} else {
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("ANALYSIS FAILED\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  cat("No conditions were successfully processed.\n")
}

cat("\nScript completed.\n")
