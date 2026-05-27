library(Seurat)
library(ggplot2)
library(dplyr)

####
# Stacked violin plot comparing NMPs and NMPs/Mesoderm-biased at 96h
# Data: 601_96h_filtered.rds
####

# Function to plot a stacked violin plot
plot_stacked_violin <- function(seurat_obj, genes_to_plot, output_dir, filename = "stacked_violin.png", 
                               cols = c("#0070ff", "#cd3278"), lw = 0, width = 4, height = 6, dpi = 600) {
  p <- VlnPlot(object = seurat_obj, features = genes_to_plot, 
               fill.by = 'ident', cols = cols, 
               stack = TRUE, flip = TRUE) + NoLegend()
  p <- p +
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = 6),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y.left = element_line(linewidth = 0.1),
      line = element_line(linewidth = 0.1)
    )
  p <- p & theme(
    strip.text = element_blank(),
    strip.text.y = element_blank(),
    strip.text.y.right = element_blank(),
    strip.background = element_blank(),
  )
  p <- thin_lines(p, lw = lw)
  ggsave(
    file.path(output_dir, filename),
    p, 
    width = width, 
    height = height, 
    units = "cm", 
    dpi = dpi
  )
  return(p)
}

# Accessory function to thin the lines of the plot
thin_lines <- function(p, lw = 0.1) {

  if (inherits(p, "ggplot")) {
    for (k in seq_along(p$layers)) {
      geom <- p$layers[[k]]$geom

      if (inherits(geom, c("GeomViolin","GeomPolygon","GeomSegment","GeomHline","GeomVline","GeomBoxplot"))) {
        p$layers[[k]]$aes_params$linewidth <- lw  # ggplot2 >= 3.4
        p$layers[[k]]$aes_params$size <- lw       # fallback for older ggplot2
      }
    }
    return(p)
  }
  
  stop("Object is not a ggplot")
}

# Create output directory
output_dir <- "MANUS_96h_nmps_stacked_violin_252312"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Load the Seurat object
obj <- readRDS("601_manual_labeling_marioni_filtering_time_based/601_96h_filtered.rds")

# Subset based on group
obj_ctrl <- subset(obj, group == "Ctrl")
obj_racl <- subset(obj, group == "RACL")
obj_naclb <- subset(obj, group == "NACL")
obj_xal <- subset(obj, group == "XAL")

# Check available cell types in the object
cat("Available cell types:\n")
print(unique(obj[["manual_marioni_label", drop = TRUE]]))

# Set identity to manual_marioni_label
Idents(obj) <- obj[["manual_marioni_label", drop = TRUE]]

# Define cell types to compare
cell_types_to_compare <- c("NMPs", "NMPs/Mesoderm-biased")
#obj_subset <- subset(obj, idents = cell_types_to_compare)

# Check if cell types exist in the object
available_types <- levels(Idents(obj))
cat("\nCell types to compare:\n")
for (ct in cell_types_to_compare) {
  if (ct %in% available_types) {
    n_cells <- sum(Idents(obj) == ct)
    cat(sprintf("  %s: %d cells\n", ct, n_cells))
  } else {
    cat(sprintf("  %s: NOT FOUND\n", ct))
  }
}

# Set identity to manual_marioni_label for each condition
Idents(obj_ctrl) <- obj_ctrl[["manual_marioni_label", drop = TRUE]]
Idents(obj_racl) <- obj_racl[["manual_marioni_label", drop = TRUE]]
Idents(obj_naclb) <- obj_naclb[["manual_marioni_label", drop = TRUE]]
Idents(obj_xal) <- obj_xal[["manual_marioni_label", drop = TRUE]]

# Subset object to only include the cell types of interest for each condition
obj_ctrl_subset <- subset(obj_ctrl, idents = cell_types_to_compare)
obj_racl_subset <- subset(obj_racl, idents = cell_types_to_compare)
obj_naclb_subset <- subset(obj_naclb, idents = cell_types_to_compare)
obj_xal_subset <- subset(obj_xal, idents = cell_types_to_compare)

# Check how many cells are in each condition and how many cells are in the cell types of interest
n_cells_ctrl <- ncol(obj_ctrl_subset)
n_cells_racl <- ncol(obj_racl_subset)
n_cells_naclb <- ncol(obj_naclb_subset)
n_cells_xal <- ncol(obj_xal_subset)

# Print the number of cells in each condition and the cell types of interest
cat("Number of cells in Control:", n_cells_ctrl, "\n")
cat("Number of cells in RACL:", n_cells_racl, "\n")
cat("Number of cells in NACLB:", n_cells_naclb, "\n")
cat("Number of cells in XAL:", n_cells_xal, "\n")

# Define genes to plot (TO BE FILLED IN)
# Add your genes of interest here
genes_to_plot <- c(
  # Example genes - replace with your genes of interest
  "T", "Sox2", "Sox1", "Prr11", "Dll1", "Tpx2", "Rspo3", "Tbx6", "Msgn1", "Cdx1", "Cdx4"
)

# Plot the stacked violin plot
plot_stacked_violin(obj_naclb_subset, genes_to_plot, output_dir, "MANUS_96h_nmps_stacked_violin_naclb.png")
plot_stacked_violin(obj_racl_subset, genes_to_plot, output_dir, "MANUS_96h_nmps_stacked_violin_racl.png")
plot_stacked_violin(obj_ctrl_subset, genes_to_plot, output_dir, "MANUS_96h_nmps_stacked_violin_ctrl.png")
plot_stacked_violin(obj_xal_subset, genes_to_plot, output_dir, "MANUS_96h_nmps_stacked_violin_xal.png")
