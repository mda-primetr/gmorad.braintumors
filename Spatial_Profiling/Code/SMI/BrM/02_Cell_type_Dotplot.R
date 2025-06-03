#----------------------------------------------------------------------------------------------------
# Source files
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")

# Set reproducibility seed
set.seed(42)

# Load slim seurat object
brm_slim <- readRDS("Input_files/SMI/BrM/brm_slim_seurat.rds") 

#----------------------------------------------------------------------------------------------------
# Pre-processing
#----------------------------------------------------------------------------------------------------

# Clean cell type names for UMAP
brm_slim@meta.data <- brm_slim@meta.data %>%
  mutate(
    cell_type_umap =
      case_when(
        manual_annotation_v6 == "Tumor" & Patient == "P35" ~ "T35",
        manual_annotation_v6 == "Tumor" & Patient == "P43" ~ "T43",
        manual_annotation_v6 == "Tumor" & Patient == "P46" ~ "T46",
        manual_annotation_v6 == "Tumor" & Patient == "P41" ~ "T41",
        manual_annotation_v6 == "Tumor" & Patient == "P40" ~ "T40",
        manual_annotation_v6 == "Tumor" & Patient == "P32" ~ "T32",
        manual_annotation_v6 == "Tumor" & Patient == "P47" ~ "T47",
        manual_annotation_v6 == "Tumor" & Patient == "P39" ~ "T39",
        manual_annotation_v6 == "Tumor" & Patient == "P33" ~ "T33",
        manual_annotation_v6 == "Mac" ~ "Mac/Mic",
        TRUE ~ manual_annotation_v6
      )
  )

# Cell type relevant markers
marker_list <- list(
  EC = c("CLDN5", "VWF", "PLVAP"),
  Mural = c("RGS5", "TAGLN", "ACTA2"),
  MSC = c("COL1A1", "FN1", "POSTN"),
  Astroglial = c("CHI3L1", "PLP1", "GFAP"),
  Erythrocyte = c("HBB", "HBA1/2"),
  Mac = c("CD68", "CD14"),
  DC = c("ITGAX"),
  Neutrophil = c("CXCL8", "S100A8", "S100A9"),
  Mast = c("TPSAB1/2", "CPA3"),
  NK = c("GNLY", "PRF1", "KLRK1"),
  B = c("MS4A1"),
  Plasma = c("IGHA1", "IGKC", "JCHAIN"),
  CD4 = c("IL7R", "CD4", "CD3D"),
  CD8 = c("CCL5", "CD8A", "CD8B"),
  Treg = c("IL2RA", "CTLA4", "FOXP3"),
  Tumor = c("ERBB2", "CCND1", "TFF3", "AGR2", "MSLN", "MUC1", "LCN2", "KRT6A", "KRT7",
    "MDK", "SCG2", "NKX2-2", "CRP", "KRT8", "SPINK1", "MKI67")
)

# Set default assay to RNA
DefaultAssay(brm_slim) <- "RNA"

# Identity mapping for correct order
identity_mapping <- c(
  "EC",
  "Mural",
  "MSC-like",
  "Astroglial",
  "Erythrocyte",
  "Mac/Mic",
  "DC",
  "Neutrophil",
  "Mast",
  "NK",
  "B cell",
  "Plasma",
  "CD4 T",
  "CD8 T",
  "Treg",
  "T32",
  "T33",
  "T35",
  "T39",
  "T40",
  "T41",
  "T43",
  "T46",
  "T47"
)

# Make sure the column is a factor with correct levels
brm_slim@meta.data$cell_type_umap <- factor(brm_slim@meta.data$cell_type_umap, levels = identity_mapping)

# Apply it to identities
Idents(brm_slim) <- brm_slim@meta.data$cell_type_umap

# Flatten gene list
gene_order <- unique(unlist(marker_list, use.names = FALSE))

#----------------------------------------------------------------------------------------------------
# Dotplot
#----------------------------------------------------------------------------------------------------

p <- DotPlot(brm_slim, features = gene_order) +
  scale_fill_gradientn(colors = c("#FFFDFD", "#99123C"), name = "Average Expression") +
  scale_size(range = c(0, 5)) +
  coord_flip() +
  labs(x = "", y = "", fill = "Average Expression") +
  theme_bw(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.ticks.x = element_blank())

# Modify layer to use fill + border
p$layers[[1]]$mapping$aes_params <- NULL
p$layers[[1]]$aes_params <- list(shape = 21,
                                 stroke = 0.3,
                                 colour = "grey30") # Border color

# Adjust mapping to use fill instead of color
p$layers[[1]]$mapping$fill <- p$layers[[1]]$mapping$colour
p$layers[[1]]$mapping$colour <- NULL

# Save plot
ggsave(
  file = "Output_files/SMI/Figures.FigureS2B.pdf",
  plot = p,
  width = 7,
  height = 9
)       
