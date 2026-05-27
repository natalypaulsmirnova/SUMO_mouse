library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
####
#UMAP projection of clusters at 72h
#Data 601_72h_filtered

#Separate UMAPs for each condition: Control, RACL, NACLB, XAL
#Size 5x5 cm
#Only color code + version with number code (as assigned in the cluster list)

obj <- readRDS("601_manual_labeling_marioni_filtering_time_based/601_72h_filtered.rds")

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
#3.Primitive Streak/5D1B85
#4.Caudal Epiblast/1F968B
#5.Epiblast/55C667
#6.Ectoderm/F8E71C
#7.Neural tube/F1B405
#8.Intermediate Mesoderm/6BB2F1
#9.PGC/B8E986
#10.Naive pluripotency/7ED321
#11.Parietal Endoderm/9013FE
#12.Visceral Endoderm/BD10E0

cell_types <- c("Anterior Primitive Streak", "Nascent mesoderm", "Primitive Streak", 
                "Caudal epiblast", "Epiblast", "Ectoderm", "Neural tube", "Intermediate mesoderm", 
                "PGC", "Naive pluripotency", "Parietal endoderm", "Visceral endoderm")
# Color code:
# Anterior Primitive Streak: #AD0314
# Nascent mesoderm: #0666E3
# Primitive Streak: #5D1B85
# Caudal epiblast: #1F968B
# Epiblast: #55C667
# Ectoderm: #F8E71C
# Neural tube: #F1B405
# Intermediate mesoderm: #6BB2F1
# PGC: #B8E986
# Naive pluripotency: #7ED321
# Parietal endoderm: #9013FE
# Visceral endoderm: #BD10E0

global_palette <- c("Anterior Primitive Streak" = "#AD0314", "Nascent mesoderm" = "#0666E3","Primitive Streak" = "#5D1B85", "Caudal epiblast" = "#1F968B", "Epiblast" = "#55C667", "Ectoderm" = "#F8E71C", "Neural tube" = "#F1B405","Intermediate mesoderm" = "#6BB2F1", 
"PGC" = "#B8E986", "Naive pluripotency" = "#7ED321", "Parietal endoderm" = "#9013FE", "Visceral endoderm" = "#BD10E0")

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

# Remove cell types that are less than 10 cells
obj_racl_filtered <- subset(obj_racl_filtered, idents = names(which(table(Idents(obj_racl_filtered)) >= 10)))
obj_naclb_filtered <- subset(obj_naclb_filtered, idents = names(which(table(Idents(obj_naclb_filtered)) >= 10)))
obj_xal_filtered <- subset(obj_xal_filtered, idents = names(which(table(Idents(obj_xal_filtered)) >= 10)))
obj_control_filtered <- subset(obj_control_filtered, idents = names(which(table(Idents(obj_control_filtered)) >= 10)))

# Reorder identity factor levels to match desired order
Idents(obj_racl_filtered) <- factor(Idents(obj_racl_filtered), levels = cell_types)
Idents(obj_naclb_filtered) <- factor(Idents(obj_naclb_filtered), levels = cell_types)
Idents(obj_xal_filtered) <- factor(Idents(obj_xal_filtered), levels = cell_types)
Idents(obj_control_filtered) <- factor(Idents(obj_control_filtered), levels = cell_types)

# Create mapping of cell types to numbers (based on order in lines 26-37)
cell_type_numbers <- 1:12
names(cell_type_numbers) <- cell_types

# Function to create UMAP plot with numbered labels (similar to MANUS_umap.r style)
make_umap_plot <- function(seurat_obj, title = NULL, point_size = 0.3, alpha = 0.6, show_numbers = TRUE) {
  
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

# Function to plot gene expression on UMAP using Seurat's FeaturePlot
make_gene_expression_plot <- function(seurat_obj, gene_name, title = NULL, point_size = 0.3, alpha = 0.6, show_legend = TRUE) {
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
output_dir <- "MANUS_UMAP_72h_601_20251203"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Generate UMAP plots with numbers
p_racl_numbers <- make_umap_plot(seurat_obj = obj_racl_filtered, title = "RACL", show_numbers = TRUE)
p_nacl_numbers <- make_umap_plot(seurat_obj = obj_naclb_filtered, title = "NACL", show_numbers = TRUE)
p_control_numbers <- make_umap_plot(seurat_obj = obj_control_filtered, title = "Control", show_numbers = TRUE)
p_xal_numbers <- make_umap_plot(seurat_obj = obj_xal_filtered, title = "XAL", show_numbers = TRUE)

# Generate and save legend plots for each condition
legend_line_gap <- 0.3  # Minimal gap for very compact legend
legend_font_size <- 14  # Increased font size

# RACL legends (with and without numbers)
legend_racl_numbers <- make_legend_plot(available_racl, circle_size = 8, font_size = legend_font_size, line_gap = legend_line_gap, show_numbers = TRUE)
legend_racl_no_numbers <- make_legend_plot(available_racl, circle_size = 8, font_size = legend_font_size, line_gap = legend_line_gap, show_numbers = FALSE)
legend_height_racl <- max(8, length(available_racl) * (0.28 * legend_line_gap))
ggsave(file.path(output_dir, "MANUS_umap_racl_72h_legend.png"), legend_racl_numbers, width = 7.5, height = legend_height_racl, units = "in", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_racl_72h_legend_no_numbers.png"), legend_racl_no_numbers, width = 7.5, height = legend_height_racl, units = "in", dpi = 600)

# NACL legends (with and without numbers)
legend_nacl_numbers <- make_legend_plot(available_nacl, circle_size = 8, font_size = legend_font_size, line_gap = legend_line_gap, show_numbers = TRUE)
legend_nacl_no_numbers <- make_legend_plot(available_nacl, circle_size = 8, font_size = legend_font_size, line_gap = legend_line_gap, show_numbers = FALSE)
legend_height_nacl <- max(8, length(available_nacl) * (0.28 * legend_line_gap))
ggsave(file.path(output_dir, "MANUS_umap_nacl_72h_legend.png"), legend_nacl_numbers, width = 7.5, height = legend_height_nacl, units = "in", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_nacl_72h_legend_no_numbers.png"), legend_nacl_no_numbers, width = 7.5, height = legend_height_nacl, units = "in", dpi = 600)

# Control legends (with and without numbers)
legend_control_numbers <- make_legend_plot(available_control, circle_size = 8, font_size = legend_font_size, line_gap = legend_line_gap, show_numbers = TRUE)
legend_control_no_numbers <- make_legend_plot(available_control, circle_size = 8, font_size = legend_font_size, line_gap = legend_line_gap, show_numbers = FALSE)
legend_height_control <- max(8, length(available_control) * (0.28 * legend_line_gap))
ggsave(file.path(output_dir, "MANUS_umap_control_72h_legend.png"), legend_control_numbers, width = 7.5, height = legend_height_control, units = "in", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_control_72h_legend_no_numbers.png"), legend_control_no_numbers, width = 7.5, height = legend_height_control, units = "in", dpi = 600)

# XAL legends (with and without numbers)
legend_xal_numbers <- make_legend_plot(available_xal, circle_size = 8, font_size = legend_font_size, line_gap = legend_line_gap, show_numbers = TRUE)
legend_xal_no_numbers <- make_legend_plot(available_xal, circle_size = 8, font_size = legend_font_size, line_gap = legend_line_gap, show_numbers = FALSE)
legend_height_xal <- max(8, length(available_xal) * (0.28 * legend_line_gap))
ggsave(file.path(output_dir, "MANUS_umap_xal_72h_legend.png"), legend_xal_numbers, width = 7.5, height = legend_height_xal, units = "in", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_xal_72h_legend_no_numbers.png"), legend_xal_no_numbers, width = 7.5, height = legend_height_xal, units = "in", dpi = 600)

# High-quality PNG outputs (numbered)
ggsave(file.path(output_dir, "MANUS_umap_racl_72h.png"), p_racl_numbers, width = 5, height = 5, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_nacl_72h.png"), p_nacl_numbers, width = 5, height = 5, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_control_72h.png"), p_control_numbers, width = 5, height = 5, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_xal_72h.png"), p_xal_numbers, width = 5, height = 5, units = "cm", dpi = 600)

# Generate UMAP plots without numbers
p_racl <- make_umap_plot(seurat_obj = obj_racl_filtered, title = "RACL", show_numbers = FALSE)
p_nacl <- make_umap_plot(seurat_obj = obj_naclb_filtered, title = "NACL", show_numbers = FALSE)
p_control <- make_umap_plot(seurat_obj = obj_control_filtered, title = "Control", show_numbers = FALSE)
p_xal <- make_umap_plot(seurat_obj = obj_xal_filtered, title = "XAL", show_numbers = FALSE)

ggsave(file.path(output_dir, "MANUS_umap_racl_72h_no_numbers.png"), p_racl, width = 5, height = 5, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_nacl_72h_no_numbers.png"), p_nacl, width = 5, height = 5, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_control_72h_no_numbers.png"), p_control, width = 5, height = 5, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_xal_72h_no_numbers.png"), p_xal, width = 5, height = 5, units = "cm", dpi = 600)

# Generate versions without axes and titles (no numbers)
p_racl_no_axes <- p_racl + theme_void() + theme(plot.title = element_blank()) + NoLegend()
p_nacl_no_axes <- p_nacl + theme_void() + theme(plot.title = element_blank()) + NoLegend()
p_control_no_axes <- p_control + theme_void() + theme(plot.title = element_blank()) + NoLegend()
p_xal_no_axes <- p_xal + theme_void() + theme(plot.title = element_blank()) + NoLegend()

ggsave(file.path(output_dir, "MANUS_umap_racl_72h_no_axes.png"), p_racl_no_axes, width = 5, height = 5, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_nacl_72h_no_axes.png"), p_nacl_no_axes, width = 5, height = 5, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_control_72h_no_axes.png"), p_control_no_axes, width = 5, height = 5, units = "cm", dpi = 600)
ggsave(file.path(output_dir, "MANUS_umap_xal_72h_no_axes.png"), p_xal_no_axes, width = 5, height = 5, units = "cm", dpi = 600)

# Gene expression plots on UMAP
genes_to_plot <- c("T", "Eomes", "Nodal", "Lhx1", "Mixl1", "Mesp1", "Dkk1", "Lefty1", 
                    "Cer1", "Msgn1", "Otx2", "Rgma", "Sox1", "Sox2", "Wnt3", "Wnt3a", 
                    "Tdgf1", "Cdx2")

# Create a list of conditions and their corresponding objects
conditions <- list(
  "racl" = obj_racl_filtered,
  "nacl" = obj_naclb_filtered,
  "control" = obj_control_filtered,
  "xal" = obj_xal_filtered
)

# Generate gene expression plots for each gene and each condition
for (gene in genes_to_plot) {
  for (cond_name in names(conditions)) {
    cond_obj <- conditions[[cond_name]]
    cond_title <- toupper(cond_name)
    if (cond_title == "NACL") cond_title <- "NACL"
    if (cond_title == "CONTROL") cond_title <- "Control"
    
    # Create the plot with legend
    p_gene <- make_gene_expression_plot(seurat_obj = cond_obj, gene_name = gene, 
                                        title = paste(cond_title, "-", gene), show_legend = TRUE)
    
    if (!is.null(p_gene)) {
      # Save with legend
      ggsave(file.path(output_dir, paste0("MANUS_umap_", cond_name, "_72h_", gene, ".png")), 
             p_gene, width = 5, height = 5, units = "cm", dpi = 600)
      
      # Save without legend
      p_gene_no_legend <- p_gene + theme(legend.position = "none")
      ggsave(file.path(output_dir, paste0("MANUS_umap_", cond_name, "_72h_", gene, "_no_legend.png")), 
             p_gene_no_legend, width = 5, height = 5, units = "cm", dpi = 600)
      
      # Save without axes and title (for cleaner presentation)
      p_gene_no_axes <- p_gene + theme_void() + theme(plot.title = element_blank())
      ggsave(file.path(output_dir, paste0("MANUS_umap_", cond_name, "_72h_", gene, "_no_axes.png")), 
             p_gene_no_axes, width = 5, height = 5, units = "cm", dpi = 600)
      
      # Save without axes, title, and legend (cleanest version)
      p_gene_no_axes_no_legend <- p_gene + theme_void() + theme(plot.title = element_blank()) + NoLegend()
      ggsave(file.path(output_dir, paste0("MANUS_umap_", cond_name, "_72h_", gene, "_no_axes_no_legend.png")), 
             p_gene_no_axes_no_legend, width = 3, height = 3, units = "cm", dpi = 600)
    }
  }
}

