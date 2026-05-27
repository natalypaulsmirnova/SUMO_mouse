library(Seurat)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(RColorBrewer)
library(grid)


# seurat_obj: a Seurat object with UMAP and your annotation column
# anno_col:   name of your annotation column (e.g., "annotation" or "celltype")
# keep_union_legend: if TRUE, shows legend with all union labels; if FALSE, only those present
# show_numbers: if TRUE, displays numeric labels at cluster centers
make_umap_plot <- function(seurat_obj, anno_col, title = NULL,
                           keep_union_legend = FALSE, point_size = 0.3, alpha = 0.8,
                           show_numbers = FALSE) {
  
  Idents(seurat_obj) <- anno_col
  
  # Enforce consistent factor order
  seurat_obj[[anno_col]] <- factor(seurat_obj[[anno_col]][,1], levels = global_labels)
  
  # Extract UMAP coordinates
  umap_df <- Embeddings(seurat_obj, "umap") %>% as.data.frame()
  
  # Detect coordinate column names
  coord_names <- colnames(umap_df)
  if (length(coord_names) < 2) stop("UMAP embedding must have at least 2 dimensions.")
  colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
  
  umap_df$label <- seurat_obj[[anno_col]][,1]
  
  # Safe summarization (handle missing values)
  centers <- umap_df %>%
    group_by(label) %>%
    summarize(
      UMAP_1 = median(UMAP_1, na.rm = TRUE),
      UMAP_2 = median(UMAP_2, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(!is.na(label))

  # Map each label to a stable numeric index matching global_labels
  label_to_index <- setNames(seq_along(global_labels), global_labels)
  centers$label_num <- label_to_index[as.character(centers$label)]
  
  # Plot
  p <- ggplot(umap_df, aes(UMAP_1, UMAP_2, color = label)) +
    geom_point(size = point_size, alpha = alpha, na.rm = TRUE) +
    scale_color_manual(
      values = global_palette,
      limits = global_labels,
      drop = !keep_union_legend,
      guide = "none"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 18),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14)
    ) +
    labs(title = title, x = "UMAP 1", y = "UMAP 2")
  
  # Add numeric labels at cluster centers if requested (repelled, non-overlapping)
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
        show.legend = FALSE
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
        show.legend = FALSE
      )
  }
  
  # Always hide legend on plots
  p <- p + NoLegend()
  
  p
}


racl <- readRDS("601_manual_labeling_marioni_filtering/601_racl_filtered.rds")
ctrl <- readRDS("601_manual_labeling_marioni_filtering/601_ctrl_filtered.rds")
nacl <- readRDS("601_manual_labeling_marioni_filtering/601_nacl_filtered.rds")
xal <- readRDS("601_manual_labeling_marioni_filtering/601_xal_filtered.rds")

set1 <- union(racl$manual_marioni_label, ctrl$manual_marioni_label)
set2 <- union(set1, nacl$manual_marioni_label)
set3 <- union(set2, xal$manual_marioni_label)

global_labels <- c(
  "Posterior somitic tissues","Embryo proper endothelium","Spinal cord progenitors",
  "Dorsal spinal cord progenitors","Cardiopharyngeal mesoderm","Neural tube","NMPs",
  "Dermomyotome","Anterior Primitive Streak","Hindbrain floor plate",
  "Ventral hindbrain progenitors","Parietal endoderm","Sclerotome","PGC","Ectoderm",
  "Hindgut","Visceral endoderm","Notochord","Kidney primordium","Gut tube",
  "Somitic mesoderm","Endotome","Foregut","NMPs/Mesoderm-biased","Presomitic mesoderm",
  "ExE ectoderm","Midgut","Naive pluripotency","Anterior somitic tissues",
  "Dorsal midbrain neurons","Epiblast","Hindbrain neural progenitors",
  "Venous endothelium","Limb mesoderm","Dorsal hindbrain progenitors",
  "Caudal epiblast","Caudal mesoderm","Migratory neural crest",
  "Intermediate mesoderm","Primitive Streak","Node","Nascent mesoderm","ExE endoderm"
)

groups <- list(
  Pluripotent = c("Naive pluripotency","Epiblast","Caudal epiblast"),
  Streak_Node_Nascent = c("Primitive Streak","Anterior Primitive Streak","Node","Nascent mesoderm"),
  NMP_axis = c("NMPs","NMPs/Mesoderm-biased","Caudal mesoderm"),
  Paraxial_Somitic = c("Presomitic mesoderm","Somitic mesoderm","Anterior somitic tissues",
                       "Posterior somitic tissues","Dermomyotome","Sclerotome",
                       "Limb mesoderm","Cardiopharyngeal mesoderm"),   # ← added here
  Intermed_Kidney = c("Intermediate mesoderm","Kidney primordium"),
  Endothelium = c("Embryo proper endothelium","Venous endothelium"),
  Endoderm = c("Parietal endoderm","Visceral endoderm","Foregut","Midgut",
               "Hindgut","Gut tube","Endotome","ExE endoderm"),
  Neural = c("Neural tube","Spinal cord progenitors","Dorsal spinal cord progenitors",
             "Hindbrain floor plate","Ventral hindbrain progenitors","Hindbrain neural progenitors",
             "Dorsal hindbrain progenitors","Dorsal midbrain neurons","Migratory neural crest"),
  Ecto_exE = c("Ectoderm","ExE ectoderm"),
  Specials = c("PGC","Notochord")
)


grad <- function(n, pal) colorRampPalette(brewer.pal(max(3, min(9, length(brewer.pal.info[pal,"maxcolors"]))), pal))(n)

cols <- c(
  setNames(grad(length(groups$Pluripotent), "Reds"), groups$Pluripotent),
  setNames(grad(length(groups$Streak_Node_Nascent), "YlOrRd"), groups$Streak_Node_Nascent),
  setNames(grad(length(groups$NMP_axis), "YlGn"), groups$NMP_axis),
  setNames(grad(length(groups$Paraxial_Somitic), "Greens"), groups$Paraxial_Somitic),
  setNames(grad(length(groups$Intermed_Kidney), "BuGn"), groups$Intermed_Kidney),
  setNames(grad(length(groups$Endothelium), "Blues"), groups$Endothelium),
  setNames(grad(length(groups$Endoderm), "PuBuGn"), groups$Endoderm),
  setNames(grad(length(groups$Neural), "Purples"), groups$Neural),
  setNames(grad(length(groups$Ecto_exE), "Greys"), groups$Ecto_exE),
  c("PGC" = "#D946EF",        # vivid magenta (stands out)
    "Notochord" = "#006D5B")   # deep teal
)



missing <- setdiff(global_labels, names(cols))
if (length(missing)) stop("Missing colors for: ", paste(missing, collapse=", "))

# This is your **global, fixed** palette (named vector)
global_palette <- cols[global_labels]


# Optional: save for reuse across scripts/projects
saveRDS(global_palette, file = "MANUS_umap_global_palette.rds")

# Function to create legend plot with numbers on top of color swatches (one per line)
make_legend_plot <- function(circle_size = 8, font_size = 10, line_gap = 1.25) {
  # Create data frame for legend
  legend_df <- data.frame(
    label = global_labels,
    label_num = seq_along(global_labels),
    color = global_palette,
    stringsAsFactors = FALSE
  )
  # vertical layout: one item per row
  legend_df$x <- 1
  legend_df$y <- rev(seq_along(global_labels)) * line_gap

  # Create the plot
  p <- ggplot(legend_df, aes(x = x, y = y)) +
    # Draw colored circles
    geom_point(aes(color = color), size = circle_size, show.legend = FALSE) +
    scale_color_identity() +
    # Add white halo for numbers (for better visibility)
    geom_text(aes(label = label_num), color = "white", fontface = "bold", 
              size = font_size * 0.4, show.legend = FALSE) +
    # Add black numbers on top
    geom_text(aes(label = label_num), color = "black", fontface = "bold", 
              size = font_size * 0.35, show.legend = FALSE) +
    # Add text labels to the right
    geom_text(aes(x = x + 0.8, y = y, label = label), 
              hjust = 0, size = font_size * 0.35, show.legend = FALSE) +
    theme_void() +
    theme(
      plot.margin = margin(20, 20, 20, 20)
    ) +
    coord_cartesian(xlim = c(0.5, 7))
  
  return(p)
}

# Generate and save legend plot with automatic height
legend_line_gap <- 1.25
legend_plot <- make_legend_plot(circle_size = 8, font_size = 10, line_gap = legend_line_gap)
legend_height <- max(8, length(global_labels) * (0.28 * legend_line_gap))  # inches
ggsave("MANUS_umap_legend.pdf", legend_plot, width = 7.5, height = legend_height, units = "in")

# Generate UMAP plots with numbers
p_racl_numbers <- make_umap_plot(seurat_obj = racl, anno_col = "manual_marioni_label",
                                  title = "RACL", keep_union_legend = TRUE, show_numbers = TRUE)
p_nacl_numbers <- make_umap_plot(seurat_obj = nacl, anno_col = "manual_marioni_label",
                                  title = "NACL", keep_union_legend = TRUE, show_numbers = TRUE)
p_ctrl_numbers <- make_umap_plot(seurat_obj = ctrl, anno_col = "manual_marioni_label",
                                  title = "CTRL", keep_union_legend = TRUE, show_numbers = TRUE)
p_xal_numbers <- make_umap_plot(seurat_obj = xal, anno_col = "manual_marioni_label",
                                 title = "XAL", keep_union_legend = TRUE, show_numbers = TRUE)

ggsave("MANUS_umap_racl_numbers.pdf", p_racl_numbers, width = 8, height = 8, units = "in")
ggsave("MANUS_umap_nacl_numbers.pdf", p_nacl_numbers, width = 8, height = 8, units = "in")
ggsave("MANUS_umap_ctrl_numbers.pdf", p_ctrl_numbers, width = 8, height = 8, units = "in")
ggsave("MANUS_umap_xal_numbers.pdf", p_xal_numbers, width = 8, height = 8, units = "in")

# High-quality PNG outputs (numbered)
ggsave("MANUS_umap_racl_numbers.png", p_racl_numbers, width = 8, height = 8, units = "in", dpi = 600)
ggsave("MANUS_umap_nacl_numbers.png", p_nacl_numbers, width = 8, height = 8, units = "in", dpi = 600)
ggsave("MANUS_umap_ctrl_numbers.png", p_ctrl_numbers, width = 8, height = 8, units = "in", dpi = 600)
ggsave("MANUS_umap_xal_numbers.png", p_xal_numbers, width = 8, height = 8, units = "in", dpi = 600)

# Generate UMAP plots without numbers
p_racl <- make_umap_plot(seurat_obj = racl, anno_col = "manual_marioni_label",
                         title = "RACL", keep_union_legend = TRUE, show_numbers = FALSE)
p_nacl <- make_umap_plot(seurat_obj = nacl, anno_col = "manual_marioni_label",
                         title = "NACL", keep_union_legend = TRUE, show_numbers = FALSE)
p_ctrl <- make_umap_plot(seurat_obj = ctrl, anno_col = "manual_marioni_label",
                         title = "CTRL", keep_union_legend = TRUE, show_numbers = FALSE)
p_xal <- make_umap_plot(seurat_obj = xal, anno_col = "manual_marioni_label",
                        title = "XAL", keep_union_legend = TRUE, show_numbers = FALSE)

ggsave("MANUS_umap_racl.pdf", p_racl, width = 8, height = 8, units = "in")
ggsave("MANUS_umap_nacl.pdf", p_nacl, width = 8, height = 8, units = "in")
ggsave("MANUS_umap_ctrl.pdf", p_ctrl, width = 8, height = 8, units = "in")
ggsave("MANUS_umap_xal.pdf", p_xal, width = 8, height = 8, units = "in")

# High-quality PNG outputs (no numbers)
ggsave("MANUS_umap_racl.png", p_racl, width = 8, height = 8, units = "in", dpi = 600)
ggsave("MANUS_umap_nacl.png", p_nacl, width = 8, height = 8, units = "in", dpi = 600)
ggsave("MANUS_umap_ctrl.png", p_ctrl, width = 8, height = 8, units = "in", dpi = 600)
ggsave("MANUS_umap_xal.png", p_xal, width = 8, height = 8, units = "in", dpi = 600)

# Generate versions without axes and titles (numbered)
p_racl_numbers_no_axes <- p_racl_numbers + theme_void() + theme(plot.title = element_blank())
p_nacl_numbers_no_axes <- p_nacl_numbers + theme_void() + theme(plot.title = element_blank())
p_ctrl_numbers_no_axes <- p_ctrl_numbers + theme_void() + theme(plot.title = element_blank())
p_xal_numbers_no_axes <- p_xal_numbers + theme_void() + theme(plot.title = element_blank())

ggsave("MANUS_umap_racl_numbers_no_axes.pdf", p_racl_numbers_no_axes, width = 8, height = 8, units = "in")
ggsave("MANUS_umap_nacl_numbers_no_axes.pdf", p_nacl_numbers_no_axes, width = 8, height = 8, units = "in")
ggsave("MANUS_umap_ctrl_numbers_no_axes.pdf", p_ctrl_numbers_no_axes, width = 8, height = 8, units = "in")
ggsave("MANUS_umap_xal_numbers_no_axes.pdf", p_xal_numbers_no_axes, width = 8, height = 8, units = "in")

ggsave("MANUS_umap_racl_numbers_no_axes.png", p_racl_numbers_no_axes, width = 8, height = 8, units = "in", dpi = 600)
ggsave("MANUS_umap_nacl_numbers_no_axes.png", p_nacl_numbers_no_axes, width = 8, height = 8, units = "in", dpi = 600)
ggsave("MANUS_umap_ctrl_numbers_no_axes.png", p_ctrl_numbers_no_axes, width = 8, height = 8, units = "in", dpi = 600)
ggsave("MANUS_umap_xal_numbers_no_axes.png", p_xal_numbers_no_axes, width = 8, height = 8, units = "in", dpi = 600)

# Generate versions without axes and titles (no numbers)
p_racl_no_axes <- p_racl + theme_void() + theme(plot.title = element_blank())
p_nacl_no_axes <- p_nacl + theme_void() + theme(plot.title = element_blank())
p_ctrl_no_axes <- p_ctrl + theme_void() + theme(plot.title = element_blank())
p_xal_no_axes <- p_xal + theme_void() + theme(plot.title = element_blank())

ggsave("MANUS_umap_racl_no_axes.pdf", p_racl_no_axes, width = 8, height = 8, units = "in")
ggsave("MANUS_umap_nacl_no_axes.pdf", p_nacl_no_axes, width = 8, height = 8, units = "in")
ggsave("MANUS_umap_ctrl_no_axes.pdf", p_ctrl_no_axes, width = 8, height = 8, units = "in")
ggsave("MANUS_umap_xal_no_axes.pdf", p_xal_no_axes, width = 8, height = 8, units = "in")

ggsave("MANUS_umap_racl_no_axes.png", p_racl_no_axes, width = 8, height = 8, units = "in", dpi = 600)
ggsave("MANUS_umap_nacl_no_axes.png", p_nacl_no_axes, width = 8, height = 8, units = "in", dpi = 600)
ggsave("MANUS_umap_ctrl_no_axes.png", p_ctrl_no_axes, width = 8, height = 8, units = "in", dpi = 600)
ggsave("MANUS_umap_xal_no_axes.png", p_xal_no_axes, width = 8, height = 8, units = "in", dpi = 600)
