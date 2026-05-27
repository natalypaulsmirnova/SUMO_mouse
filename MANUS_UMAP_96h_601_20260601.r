library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
####
#UMAP projection of clusters at 96h
#Data 601_96h_filtered

#Separate UMAPs for each condition: Control, RACL, NACLB, XAL
#Size 5x5 cm
#Only color code + version with number code (as assigned in the cluster list)

obj <- readRDS("601_manual_labeling_marioni_filtering_time_based/601_96h_filtered.rds")

# Remove ExE ectoderm and ExE endoderm from the object
# First, set idents to the manual_marioni_label to access cell types
Idents(obj) <- obj[["manual_marioni_label", drop = TRUE]]
# Remove ExE ectoderm and ExE endoderm
exclude_types <- c("ExE ectoderm", "ExE endoderm")
obj <- subset(obj, idents = setdiff(levels(Idents(obj)), exclude_types))

# RACL only
obj_racl <- subset(obj, group == "RACL")

# NACLB only
obj_naclb <- subset(obj, group == "NACL")

# XAL only
obj_xal <- subset(obj, group == "XAL")

# Control only
obj_control <- subset(obj, group == "Ctrl")

# Select cell types to plot
#1.Anterior primitive streak/AD0314
#2.Nascent mesoderm/0666E3
#3.Cardiopharingeal Mesoderm/6F032E
#4.Primitive Streak/5D1B85
#5.Caudal Epiblast/1F968B
#6.NMPs/0070FF
#7.NMPs/Mesoderm-biased/CD3278
#8.Caudal Mesoderm/A2FC3C
#9.Presomitic mesoderm/9D6CD2
#10.Somitic mesoderm/A5C3EB
#11.Posterior somitic tissue/9D6CD2
#12.Dermamyotome/3E3993
#13.Embryo proper endothelium/04B311
#14.Epiblast/55C667
#15.Ectoderm/F8E71C
#16.Neural tube/F1B405
#17.Spinal cord progenitors/FFB420
#18.Dorsal spinal cord progenitors/96DDB5
#19.Dorsal hindbrain progenitors/F6BB97
#20.Ventral hindbrain progenitors/8EFE49
#21.Dorsal midbrain neurones/359EAA
#22.Hindbrain neural progenitors/A92395
#23.Hindbrain floor plate/CC489A
#24.Intermediate Mesoderm/6BB2F1
#25.Kidney primordium/4E81AD
#26.Limb mesoderm/0D0D88
#27.Node/EE91A9
#28.Notochord/D45172
#29.Foregut/F89441
#30.Gut tube/E56B5D
#31.Midgut/4C02A1
#32.Hindgut/7E03A8
#33.PGC/B8E986
#34.Naive pluripotency/7ED321
#35.Parietal Endoderm/9013FE
#36.Visceral Endoderm/BD10E0

cell_types <- c("Anterior Primitive Streak", "Nascent mesoderm", "Cardiopharingeal mesoderm", 
                "Primitive Streak", "Caudal epiblast", "NMPs", "NMPs/Mesoderm-biased", 
                "Caudal mesoderm", "Presomitic mesoderm", "Somitic mesoderm", "Posterior somitic tissue", 
                "Dermamyotome", "Embryo proper endothelium", "Epiblast", "Ectoderm", "Neural tube", 
                "Spinal cord progenitors", "Dorsal spinal cord progenitors", "Dorsal hindbrain progenitors", 
                "Ventral hindbrain progenitors", "Dorsal midbrain neurones", "Hindbrain neural progenitors", 
                "Hindbrain floor plate", "Intermediate mesoderm", "Kidney primordium", "Limb mesoderm", 
                "Node", "Notochord", "Foregut", "Gut tube", "Midgut", "Hindgut", "PGC", 
                "Naive pluripotency", "Parietal endoderm", "Visceral endoderm")
# Color code:
# Anterior Primitive Streak: #AD0314
# Nascent mesoderm: #0666E3
# Cardiopharingeal mesoderm: #6F032E
# Primitive Streak: #5D1B85
# Caudal epiblast: #1F968B
# NMPs: #0070FF
# NMPs/Mesoderm-biased: #CD3278
# Caudal mesoderm: #A2FC3C
# Presomitic mesoderm: #9D6CD2
# Somitic mesoderm: #A5C3EB
# Posterior somitic tissue: #9D6CD2
# Dermamyotome: #3E3993
# Embryo proper endothelium: #04B311
# Epiblast: #55C667
# Ectoderm: #F8E71C
# Neural tube: #F1B405
# Spinal cord progenitors: #FFB420
# Dorsal spinal cord progenitors: #96DDB5
# Dorsal hindbrain progenitors: #F6BB97
# Ventral hindbrain progenitors: #8EFE49
# Dorsal midbrain neurones: #359EAA
# Hindbrain neural progenitors: #A92395
# Hindbrain floor plate: #CC489A
# Intermediate mesoderm: #6BB2F1
# Kidney primordium: #4E81AD
# Limb mesoderm: #0D0D88
# Node: #EE91A9
# Notochord: #D45172
# Foregut: #F89441
# Gut tube: #E56B5D
# Midgut: #4C02A1
# Hindgut: #7E03A8
# PGC: #B8E986
# Naive pluripotency: #7ED321
# Parietal endoderm: #9013FE
# Visceral endoderm: #BD10E0

global_palette <- c("Anterior Primitive Streak" = "#AD0314", "Nascent mesoderm" = "#0666E3", 
                    "Cardiopharingeal mesoderm" = "#6F032E", "Primitive Streak" = "#5D1B85", 
                    "Caudal epiblast" = "#1F968B", "NMPs" = "#0070FF", "NMPs/Mesoderm-biased" = "#CD3278", 
                    "Caudal mesoderm" = "#A2FC3C", "Presomitic mesoderm" = "#9D6CD2", 
                    "Somitic mesoderm" = "#A5C3EB", "Posterior somitic tissue" = "#9D6CD2", 
                    "Dermamyotome" = "#3E3993", "Embryo proper endothelium" = "#04B311", 
                    "Epiblast" = "#55C667", "Ectoderm" = "#F8E71C", "Neural tube" = "#F1B405", 
                    "Spinal cord progenitors" = "#FFB420", "Dorsal spinal cord progenitors" = "#96DDB5", 
                    "Dorsal hindbrain progenitors" = "#F6BB97", "Ventral hindbrain progenitors" = "#8EFE49", 
                    "Dorsal midbrain neurones" = "#359EAA", "Hindbrain neural progenitors" = "#A92395", 
                    "Hindbrain floor plate" = "#CC489A", "Intermediate mesoderm" = "#6BB2F1", 
                    "Kidney primordium" = "#4E81AD", "Limb mesoderm" = "#0D0D88", "Node" = "#EE91A9", 
                    "Notochord" = "#D45172", "Foregut" = "#F89441", "Gut tube" = "#E56B5D", 
                    "Midgut" = "#4C02A1", "Hindgut" = "#7E03A8", "PGC" = "#B8E986", 
                    "Naive pluripotency" = "#7ED321", "Parietal endoderm" = "#9013FE", 
                    "Visceral endoderm" = "#BD10E0")

# Change the idents to the cell_types
Idents(obj_racl) <- obj_racl[["manual_marioni_label", drop = TRUE]]
Idents(obj_naclb) <- obj_naclb[["manual_marioni_label", drop = TRUE]]
Idents(obj_xal) <- obj_xal[["manual_marioni_label", drop = TRUE]]
Idents(obj_control) <- obj_control[["manual_marioni_label", drop = TRUE]]

# Subset the objects to only include the desired cell types
# Only use cell types that actually exist in each object to avoid errors
obj_racl_filtered <- subset(obj_racl, idents = intersect(cell_types, levels(Idents(obj_racl))))
obj_naclb_filtered <- subset(obj_naclb, idents = intersect(cell_types, levels(Idents(obj_naclb))))
obj_xal_filtered <- subset(obj_xal, idents = intersect(cell_types, levels(Idents(obj_xal))))
obj_control_filtered <- subset(obj_control, idents = intersect(cell_types, levels(Idents(obj_control))))


# Reorder identity factor levels to match desired order
Idents(obj_racl_filtered) <- factor(Idents(obj_racl_filtered), levels = cell_types)
Idents(obj_naclb_filtered) <- factor(Idents(obj_naclb_filtered), levels = cell_types)
Idents(obj_xal_filtered) <- factor(Idents(obj_xal_filtered), levels = cell_types)
Idents(obj_control_filtered) <- factor(Idents(obj_control_filtered), levels = cell_types)

# Create mapping of cell types to numbers (based on order in lines 27-62)
cell_type_numbers <- 1:36
names(cell_type_numbers) <- cell_types

# Function to create UMAP plot with numbered labels (similar to MANUS_umap.r style)
make_umap_plot <- function(seurat_obj, title = NULL, point_size = 0.3, alpha = 0.6, show_numbers = TRUE, xlim = NULL, ylim = NULL) {
  
  # Extract UMAP coordinates
  umap_df <- Embeddings(seurat_obj, "umap") %>% as.data.frame()
  
  # Detect coordinate column names
  coord_names <- colnames(umap_df)
  if (length(coord_names) < 2) stop("UMAP embedding must have at least 2 dimensions.")
  colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
  
  # Get cell labels
  umap_df$label <- Idents(seurat_obj)
  
  # Enforce consistent factor order
  umap_df$label <- factor(umap_df$label, levels = cell_types)
  
  # Safe summarization (handle missing values)
  centers <- umap_df %>%
    group_by(label) %>%
    summarize(
      UMAP_1 = median(UMAP_1, na.rm = TRUE),
      UMAP_2 = median(UMAP_2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(!is.na(label))
  
  # Map each label to a stable numeric index matching cell_types
  label_to_index <- setNames(seq_along(cell_types), cell_types)
  centers$label_num <- label_to_index[as.character(centers$label)]
  
  # Plot
  p <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = label)) +
    geom_point(size = point_size, alpha = alpha, na.rm = TRUE) +
    scale_color_manual(
      values = global_palette,
      limits = cell_types,
      drop = TRUE
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14)
    ) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2")
  
  # Apply axis limits if provided
  if (!is.null(xlim)) {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  } else if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }
  
  # Add numeric labels at cluster centers if requested (repelled, non-overlapping with white halo)
  if (show_numbers && nrow(centers)) {
    set.seed(123)
    p <- p +
      ggrepel::geom_text_repel(
        data = centers,
        aes(x = UMAP_1, y = UMAP_2, label = label_num),
        color = "white",
        size = 5.6,
        fontface = "bold",
        box.padding = 0.35,
        point.padding = 0.35,
        max.overlaps = Inf,
        segment.size = 0,
        min.segment.length = Inf,
        show.legend = FALSE,
        inherit.aes = FALSE
      ) +
      ggrepel::geom_text_repel(
        data = centers,
        aes(x = UMAP_1, y = UMAP_2, label = label_num),
        color = "black",
        size = 5.2,
        fontface = "bold",
        box.padding = 0.35,
        point.padding = 0.35,
        max.overlaps = Inf,
        segment.size = 0,
        min.segment.length = Inf,
        show.legend = FALSE,
        inherit.aes = FALSE
      )
  }
  
  # Always hide legend on plots
  p <- p + NoLegend()
  
  return(p)
}

# Function to create UMAP plot highlighting only NMPs and NMPs/Mesoderm-biased
make_nmps_highlight_plot <- function(seurat_obj, title = NULL, point_size = 0.3, alpha = 0.6, xlim = NULL, ylim = NULL) {
  
  # Extract UMAP coordinates
  umap_df <- Embeddings(seurat_obj, "umap") %>% as.data.frame()
  
  # Detect coordinate column names
  coord_names <- colnames(umap_df)
  if (length(coord_names) < 2) stop("UMAP embedding must have at least 2 dimensions.")
  colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
  
  # Get cell labels
  umap_df$label <- Idents(seurat_obj)
  
  # Create color mapping: NMPs and NMPs/Mesoderm-biased get their colors, others get grey
  nmps_colors <- c("NMPs" = "#0070FF", "NMPs/Mesoderm-biased" = "#CD3278")
  grey_color <- "#dfe3ee"
  
  # Create color vector for each cell
  umap_df$color <- ifelse(umap_df$label %in% names(nmps_colors), 
                          nmps_colors[as.character(umap_df$label)], 
                          grey_color)
  
  # Create alpha vector: lower alpha for grey points, full alpha for NMPs
  grey_alpha <- 0.2  # Reduced alpha for grey points
  umap_df$point_alpha <- ifelse(umap_df$label %in% names(nmps_colors), 
                                alpha, 
                                grey_alpha)
  
  # Plot
  p <- ggplot(umap_df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = color, alpha = point_alpha), size = point_size, na.rm = TRUE) +
    scale_color_identity() +
    scale_alpha_identity() +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14)
    ) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2") +
    NoLegend()
  
  # Apply axis limits if provided
  if (!is.null(xlim)) {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  } else if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }
  
  return(p)
}

# Function to plot gene expression on UMAP using Seurat's FeaturePlot
make_gene_expression_plot <- function(seurat_obj, gene_name, title = NULL, point_size = 0.3, alpha = 0.6, show_legend = TRUE, xlim = NULL, ylim = NULL) {
  # Check if gene exists
  if (!gene_name %in% rownames(seurat_obj)) {
    warning(paste("Gene", gene_name, "not found in the object. Skipping..."))
    return(NULL)
  }
  
  # Use FeaturePlot with custom settings
  # FeaturePlot returns a list when combine=TRUE, so we extract the first (and only) plot
  p_list <- FeaturePlot(
    seurat_obj,
    features = gene_name,
    reduction = "umap",
    pt.size = point_size,
    alpha = alpha,
    cols = c("#dfe3ee", "#011f4b"),  # Min to Max colors
    combine = FALSE  # Get list of plots
  )
  
  # Extract the plot (FeaturePlot returns a list)
  p <- p_list[[1]]
  
  # Customize the plot
  p <- p +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14)
    ) +
    labs(
      title = ifelse(is.null(title), gene_name, title),
      x = "UMAP 1",
      y = "UMAP 2"
    )
  
  # Apply axis limits if provided
  if (!is.null(xlim)) {
    p <- p + coord_cartesian(xlim = xlim, ylim = ylim)
  } else if (!is.null(ylim)) {
    p <- p + coord_cartesian(ylim = ylim)
  }
  
  # Control legend visibility
  if (!show_legend) {
    p <- p + NoLegend()
  } else {
    p <- p + theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(size = 5),
      legend.key.size = unit(0.4, "cm"),
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(0.3, "cm")
    ) +
    guides(color = guide_colorbar(
      title = "",
      title.position = "top",
      barwidth = 0.4
    ))
  }
  
  return(p)
}

# Function to create legend plot with or without numbers on top of color swatches
# available_types: vector of cell types present in this condition (resets numbering to 1, 2, 3...)
# show_numbers: if TRUE, displays numbers on swatches; if FALSE, just colors
make_legend_plot <- function(available_types, circle_size = 8, font_size = 14, line_gap = 0.3, show_numbers = TRUE) {
  # Filter to only include available types, order by original numbering but reset to 1, 2, 3...
  legend_df <- data.frame(
    label = available_types,
    original_num = cell_type_numbers[available_types],  # Keep original for ordering
    color = global_palette[available_types],
    stringsAsFactors = FALSE
  )
  # Order by the original number
  legend_df <- legend_df[order(legend_df$original_num), ]
  
  # Reset numbering to 1, 2, 3...
  legend_df$label_num <- seq_along(legend_df$label)
  
  # vertical layout: one item per row
  legend_df$x <- 1
  legend_df$y <- rev(seq_along(legend_df$label)) * line_gap
  
  # Create the plot
  p <- ggplot(legend_df, aes(x = x, y = y)) +
    # Draw colored circles
    geom_point(aes(color = color), size = circle_size, show.legend = FALSE) +
    scale_color_identity()
  
  # Add numbers on swatches only if requested
  if (show_numbers) {
    p <- p +
      # Add white halo for numbers (for better visibility)
      geom_text(aes(label = label_num), color = "white", fontface = "bold", 
                size = font_size * 0.4, show.legend = FALSE) +
      # Add black numbers on top
      geom_text(aes(label = label_num), color = "black", fontface = "bold", 
                size = font_size * 0.35, show.legend = FALSE)
  }
  
  # Add text labels to the right
  p <- p +
    geom_text(aes(x = x + 0.8, y = y, label = label), 
              hjust = 0, size = font_size * 0.35, show.legend = FALSE) +
    theme_void() +
    theme(
      plot.margin = margin(20, 20, 20, 20)
    ) +
    coord_cartesian(xlim = c(0.5, 7))
  
  return(p)
}

# Get available cell types for each condition (preserving original numbering)
get_available_types <- function(seurat_obj) {
  available <- levels(Idents(seurat_obj))
  available <- available[available %in% names(cell_type_numbers)]
  return(available)
}

available_racl <- get_available_types(obj_racl_filtered)
available_nacl <- get_available_types(obj_naclb_filtered)
available_control <- get_available_types(obj_control_filtered)
available_xal <- get_available_types(obj_xal_filtered)

# Create output directory with script name
output_dir <- "MANUS_UMAP_96h_601_20260601"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Calculate global UMAP axis limits across all conditions (must be done before creating plots)
all_umap_coords <- list()
all_umap_coords[[1]] <- Embeddings(obj_racl_filtered, "umap")[, 1:2]
all_umap_coords[[2]] <- Embeddings(obj_naclb_filtered, "umap")[, 1:2]
all_umap_coords[[3]] <- Embeddings(obj_control_filtered, "umap")[, 1:2]
all_umap_coords[[4]] <- Embeddings(obj_xal_filtered, "umap")[, 1:2]

# Combine all coordinates
all_coords <- do.call(rbind, all_umap_coords)
colnames(all_coords) <- c("UMAP_1", "UMAP_2")

# Calculate global min and max for both axes
global_xlim <- c(min(all_coords[, "UMAP_1"], na.rm = TRUE), max(all_coords[, "UMAP_1"], na.rm = TRUE))
global_ylim <- c(min(all_coords[, "UMAP_2"], na.rm = TRUE), max(all_coords[, "UMAP_2"], na.rm = TRUE))

# Generate UMAP plots with numbers
p_racl_numbers <- make_umap_plot(seurat_obj = obj_racl_filtered, title = "RACL", show_numbers = TRUE, xlim = global_xlim, ylim = global_ylim)
p_nacl_numbers <- make_umap_plot(seurat_obj = obj_naclb_filtered, title = "NACL", show_numbers = TRUE, xlim = global_xlim, ylim = global_ylim)
p_control_numbers <- make_umap_plot(seurat_obj = obj_control_filtered, title = "Control", show_numbers = TRUE, xlim = global_xlim, ylim = global_ylim)
p_xal_numbers <- make_umap_plot(seurat_obj = obj_xal_filtered, title = "XAL", show_numbers = TRUE, xlim = global_xlim, ylim = global_ylim)

# Generate and save a single legend plot with all 36 cell types
legend_line_gap <- 1.1  # Minimal gap for very compact legend
legend_font_size <- 16  # Increased font size

# Create legend with all 36 cell types (numbered 1-36)
all_cell_types <- cell_types  # Use all 36 cell types
legend_df <- data.frame(
  label = all_cell_types,
  label_num = 1:36,  # Use original numbering 1-36
  color = global_palette[all_cell_types],
  stringsAsFactors = FALSE
)

# Vertical layout: one item per row
legend_df$x <- 1
legend_df$y <- rev(seq_along(legend_df$label)) * legend_line_gap

# Create legend plot with numbers
legend_all_numbers <- ggplot(legend_df, aes(x = x, y = y)) +
  # Draw colored circles
  geom_point(aes(color = color), size = 8, show.legend = FALSE) +
  scale_color_identity() +
  # Add white halo for numbers (for better visibility)
  #geom_text(aes(label = label_num), color = "white", fontface = "bold", 
  #          size = legend_font_size * 0.35, show.legend = FALSE) +
  # Add black numbers on top
  geom_text(aes(label = label_num), color = "black", fontface = "bold", 
            size = legend_font_size * 0.3, show.legend = FALSE) +
  # Add text labels to the right
  geom_text(aes(x = x + 0.8, y = y, label = label), 
            hjust = 0, size = legend_font_size * 0.35, show.legend = FALSE) +
  theme_void() +
  theme(
    plot.margin = margin(20, 20, 20, 20)
  ) +
  coord_cartesian(xlim = c(0.5, 7))

# Create legend plot without numbers
legend_all_no_numbers <- ggplot(legend_df, aes(x = x, y = y)) +
  # Draw colored circles
  geom_point(aes(color = color), size = 8, show.legend = FALSE) +
  scale_color_identity() +
  # Add text labels to the right
  geom_text(aes(x = x + 0.8, y = y, label = label), 
            hjust = 0, size = legend_font_size * 0.35, show.legend = FALSE) +
  theme_void() +
  theme(
    plot.margin = margin(20, 20, 20, 20)
  ) +
  coord_cartesian(xlim = c(0.5, 7))

# Calculate height for all 36 cell types
legend_height_all <- max(8, length(all_cell_types) * (0.28 * legend_line_gap))

# Save the legend plots
ggsave(file.path(output_dir, "MANUS_umap_96h_legend_all.png"), legend_all_numbers, width = 7.5, height = legend_height_all, units = "in", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_96h_legend_all_no_numbers.png"), legend_all_no_numbers, width = 7.5, height = legend_height_all, units = "in", dpi = 600)

# High-quality PNG outputs (numbered)
ggsave(file.path(output_dir, "MANUS_umap_racl_96h.png"), p_racl_numbers, width = 5, height = 5, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_nacl_96h.png"), p_nacl_numbers, width = 5, height = 5, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_control_96h.png"), p_control_numbers, width = 5, height = 5, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_xal_96h.png"), p_xal_numbers, width = 5, height = 5, units = "cm", dpi = 600)

# Generate UMAP plots without numbers
p_racl <- make_umap_plot(seurat_obj = obj_racl_filtered, title = "RACL", show_numbers = FALSE, xlim = global_xlim, ylim = global_ylim)
p_nacl <- make_umap_plot(seurat_obj = obj_naclb_filtered, title = "NACL", show_numbers = FALSE, xlim = global_xlim, ylim = global_ylim)
p_control <- make_umap_plot(seurat_obj = obj_control_filtered, title = "Control", show_numbers = FALSE, xlim = global_xlim, ylim = global_ylim)
p_xal <- make_umap_plot(seurat_obj = obj_xal_filtered, title = "XAL", show_numbers = FALSE, xlim = global_xlim, ylim = global_ylim)

ggsave(file.path(output_dir, "MANUS_umap_racl_96h_no_numbers.png"), p_racl, width = 5, height = 5, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_nacl_96h_no_numbers.png"), p_nacl, width = 5, height = 5, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_control_96h_no_numbers.png"), p_control, width = 5, height = 5, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_xal_96h_no_numbers.png"), p_xal, width = 5, height = 5, units = "cm", dpi = 600)

# Generate versions without axes and titles (no numbers)
p_racl_no_axes <- p_racl + theme_void() + theme(plot.title = element_blank()) + NoLegend()
p_nacl_no_axes <- p_nacl + theme_void() + theme(plot.title = element_blank()) + NoLegend()
p_control_no_axes <- p_control + theme_void() + theme(plot.title = element_blank()) + NoLegend()
p_xal_no_axes <- p_xal + theme_void() + theme(plot.title = element_blank()) + NoLegend()

ggsave(file.path(output_dir, "MANUS_umap_racl_96h_no_axes.png"), p_racl_no_axes, width = 5, height = 5, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_nacl_96h_no_axes.png"), p_nacl_no_axes, width = 5, height = 5, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_control_96h_no_axes.png"), p_control_no_axes, width = 5, height = 5, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_xal_96h_no_axes.png"), p_xal_no_axes, width = 5, height = 5, units = "cm", dpi = 600)

# Generate NMPs highlight plots (only NMPs and NMPs/Mesoderm-biased colored, rest grey)
p_racl_nmps <- make_nmps_highlight_plot(seurat_obj = obj_racl_filtered, title = "RACL", xlim = global_xlim, ylim = global_ylim)
p_nacl_nmps <- make_nmps_highlight_plot(seurat_obj = obj_naclb_filtered, title = "NACL", xlim = global_xlim, ylim = global_ylim)
p_control_nmps <- make_nmps_highlight_plot(seurat_obj = obj_control_filtered, title = "Control", xlim = global_xlim, ylim = global_ylim)
p_xal_nmps <- make_nmps_highlight_plot(seurat_obj = obj_xal_filtered, title = "XAL", xlim = global_xlim, ylim = global_ylim)

# Save NMPs highlight plots (4x4 cm)
ggsave(file.path(output_dir, "MANUS_umap_racl_96h_nmps.png"), p_racl_nmps, width = 4, height = 4, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_nacl_96h_nmps.png"), p_nacl_nmps, width = 4, height = 4, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_control_96h_nmps.png"), p_control_nmps, width = 4, height = 4, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_xal_96h_nmps.png"), p_xal_nmps, width = 4, height = 4, units = "cm", dpi = 600)

# Generate versions without axes and titles
p_racl_nmps_no_axes <- p_racl_nmps + theme_void() + theme(plot.title = element_blank())
p_nacl_nmps_no_axes <- p_nacl_nmps + theme_void() + theme(plot.title = element_blank())
p_control_nmps_no_axes <- p_control_nmps + theme_void() + theme(plot.title = element_blank())
p_xal_nmps_no_axes <- p_xal_nmps + theme_void() + theme(plot.title = element_blank())

ggsave(file.path(output_dir, "MANUS_umap_racl_96h_nmps_no_axes.png"), p_racl_nmps_no_axes, width = 4, height = 4, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_nacl_96h_nmps_no_axes.png"), p_nacl_nmps_no_axes, width = 4, height = 4, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_control_96h_nmps_no_axes.png"), p_control_nmps_no_axes, width = 4, height = 4, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_xal_96h_nmps_no_axes.png"), p_xal_nmps_no_axes, width = 4, height = 4, units = "cm", dpi = 600)


