library(Seurat)
library(viridis)
library(ggplot2)


obj <- readRDS("601_manual_labeling_marioni_filtering_time_based/601_72h_filtered.rds")

genes <- c("T","Eomes", "Mixl1", "Gsc", "Foxa2", "Lhx1", "Nodal", "Tdgf1", "Wnt3", 
           "Fgf5", "Otx2", "Trh", "Flt1", "Sox17", "Kdr", "Mesp1", "Tbx6", "Msgn1", 
           "Cdh1", "Cdh2", "Wnt3a", "Wnt8a", "Epha1", "Nkx1-2", "Cdx1", "Cdx2", "Sox1",
           "Rgma", "Sox2", "Utf1", "Dppa5a", "Dnmt3a", "Zfp42", "Dppa3", "Klf2", "Afp", 
           "Ttr", "Gata6", "Epcam", "Dab2", "Snai1", "Vim", "Pth1r", "Pdgfra")
Idents(obj) <- obj[["manual_marioni_label", drop = TRUE]]
DefaultAssay(obj) <- "RNA"

# Filter to only show specific cell types and maintain order
cell_types <- c("Anterior Primitive Streak", "Nascent mesoderm", "Primitive Streak", 
                "Caudal epiblast", "Neural tube", "Ectoderm", "Epiblast", "PGC", 
                "Visceral endoderm", "Parietal endoderm")

# Subset object to only include desired cell types
obj_filtered <- subset(obj, idents = cell_types)

# Reorder identity factor levels to match desired order
Idents(obj_filtered) <- factor(Idents(obj_filtered), levels = cell_types)

# Create DotPlot with viridis color scheme (yellow to purple gradient)
dp <- DotPlot(object = obj_filtered, features = genes,
              assay = "RNA", cluster.idents = FALSE, cols = viridis(2))

# Plot 1: Standard viridis (yellow = high, purple = low)
dp1 <- dp + scale_y_discrete(limits = rev(cell_types)) + 
           scale_colour_viridis_c(option = "viridis") + 
            ylab("Cell types") + xlab("Genes") +
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12))

# Plot 2: Reversed viridis (purple = high, yellow = low)
dp2 <- dp + scale_y_discrete(limits = rev(cell_types)) + 
           scale_colour_viridis_c(option = "viridis", direction = -1) + 
            ylab("Cell types") + xlab("Genes") +
           theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 12))

# Get base filename from script name
script_name <- "MANUS_dotplot_req_20251203"  # Base name for output files

# Save plots as high-resolution PNG files
ggsave(filename = paste0(script_name, "_viridis.png"), 
       plot = dp1, 
       width = 10.3, 
       height = 3.7, 
       dpi = 300, 
       units = "in")

ggsave(filename = paste0(script_name, "_viridis_reversed.png"), 
       plot = dp2, 
       width = 10.3, 
       height = 3.7, 
       dpi = 300, 
       units = "in")

# Display plots
print(dp1)
print(dp2)






