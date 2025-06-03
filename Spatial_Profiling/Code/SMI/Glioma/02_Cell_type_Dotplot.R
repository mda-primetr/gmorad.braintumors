#----------------------------------------------------------------------------------------------------
# Source files
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")

# Set reproducibility seed
set.seed(42)

# Load slim seurat object
glioma_slim <- readRDS("Input_files/SMI/Glioma/glioma_slim_seurat.rds") 

#----------------------------------------------------------------------------------------------------
# Pre-processing
#----------------------------------------------------------------------------------------------------

# Clean cell type names for UMAP

# Clean cell type names for UMAP
glioma_slim@meta.data <- glioma_slim@meta.data %>%
  mutate(
    cell_type_umap =
      case_when(
        cell_type == "Tumor" & Patient == "P2" ~ "T2",
        cell_type == "Tumor" & Patient == "P3" ~ "T3",
        cell_type == "Tumor" & Patient == "P4" ~ "T4",
        cell_type == "Tumor" & Patient == "P6" ~ "T6",
        cell_type == "Tumor" & Patient == "P7" ~ "T7",
        cell_type == "Tumor" & Patient == "P8" ~ "T8",
        cell_type == "Tumor" & Patient == "P10" ~ "T10",
        cell_type == "Tumor" & Patient == "P11" ~ "T11",
        cell_type == "Tumor" & Patient == "P12" ~ "T12",
        cell_type == "Tumor" & Patient == "P15" ~ "T15",
        cell_type == "Tumor" & Patient == "P16" ~ "T16",
        cell_type == "Tumor" & Patient == "P17" ~ "T17",
        cell_type == "Tumor" & Patient == "P20" ~ "T20",
        cell_type == "Tumor" & Patient == "P22" ~ "T22",
        cell_type == "Tumor" & Patient == "P24" ~ "T24",
        cell_type == "Tumor" & Patient == "P25" ~ "T25",
        cell_type == "Tumor" & Patient == "P26" ~ "T26",
        cell_type == "Tumor" & Patient == "P29" ~ "T29",
        TRUE ~ cell_type  
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
  DC = c("ITGAX","HA"),
  Neutrophil = c("CXCL8", "S100A8", "S100A9"),
  Mast = c("TPSAB1/2", "CPA3"),
  NK = c("GNLY", "PRF1", "KLRK1"),
  B = c("MS4A1"),
  Plasma = c("IGHA1", "IGKC", "JCHAIN"),
  CD4 = c("IL7R", "CD4", "CD3D"),
  CD8 = c("CCL5", "CD8A", "CD8B"),
  Treg = c("IL2RA", "CTLA4", "FOXP3"),
  Tumor = c("VEGFA", "EGFR", "CRYAB")
)

# Set default assay to RNA
DefaultAssay(glioma_slim) <- "RNA"

# Identity mapping for correct order
identity_mapping <- c(
  "Vasculature",
  "MSC-like",
  "Astrocyte",
  "Neuron",
  "Oligodendrocyte",
  "OPC",
  "Erythrocytes",
  "Mac/Mic",
  "Mono",
  "DC",
  "Neutrophil",
  "NK",
  "Plasma",
  "CD8 T",
  "Treg",
  "T2",
  "T3",
  "T4",
  "T6",
  "T7",
  "T8",
  "T10",
  "T11",
  "T12",
  "T15",
  "T16",
  "T17",
  "T20",
  "T22",
  "T24",
  "T25",
  "T26",
  "T29"
)

# Make sure the column is a factor with correct levels
glioma_slim@meta.data$cell_type_umap <- factor(glioma_slim@meta.data$cell_type_umap, levels = identity_mapping)

# Apply to identities
Idents(glioma_slim) <- glioma_slim@meta.data$cell_type_umap

# Flatten gene list
gene_order <- unique(unlist(marker_list, use.names = FALSE))

#----------------------------------------------------------------------------------------------------
# Dotplot
#----------------------------------------------------------------------------------------------------

p <- DotPlot(glioma_slim, features = gene_order) +
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
  file = "Output_files/SMI/Figures.FigureS2A.pdf",
  plot = p,
  width = 7,
  height = 9
)        
