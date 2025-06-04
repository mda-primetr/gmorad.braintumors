#----------------------------------------------------------------------------------------------------
# Source files
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")

# Set reproducibility seed
set.seed(42)

# Load slim seurat object
glioma_slim <- readRDS("Input_files/SMI/Glioma/glioma_slim_seurat.rds") 

#----------------------------------------------------------------------------------------------------
# UMAP with cell types
#----------------------------------------------------------------------------------------------------

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

# Assign colors to cell types
custom_colors <- c(
  "Astrocyte"      = "#88CCEE",
  "CD8 T"          = "#764E9F",   
  "DC"             = "#6B7FD7",   
  "Erythrocytes"   = "#DDCC77",
  "Mac/Mic"         = "#6A2013",   
  "Mono"           = "#AA4499",
  "MSC-like"       = "#F15167",   
  "Neuron"         = "#0072B2",   
  "Neutrophil"     = "darkblue",
  "NK"             = "#CC61B0",   
  "Oligodendrocyte"= "#D2691E",   
  "OPC"            = "#FC8D62",   
  "Plasma"         = "#3D2842",   
  "T10"            = "#E69F00",
  "T11"            = "#56B4E9",
  "T12"            = "#009E73",
  "T15"            = "#F0E442",
  "T16"            = "#708090",   
  "T17"            = "#BDB76B",   
  "T20"            = "#D55E00",
  "T22"            = "#CC79A7",
  "T24"            = "#66C2A5",
  "T25"            = "#8DA0CB",
  "T26"            = "#8FBC8F",
  "T29"            = "#E78AC3",
  "T3"             = "#A6D854",
  "T4"             = "#FFD92F",
  "T6"             = "#E5C494",
  "T7"             = "#B3B3B3",
  "T8"             = "#006400",   
  "Treg"           = "#E58606",   
  "Vasculature"    = "#99C945"    
)

# Get embeddings
umap_df <- glioma_slim@reductions$umap@cell.embeddings %>%
  as.data.frame() %>% tibble::rownames_to_column("cell_ID")

# Get relevant metadata
meta_df <- glioma_slim@meta.data

# Join embeddings and metadata
plot_df <- left_join(meta_df, umap_df, by = "cell_ID") %>%
  mutate(cell_type = factor(cell_type_umap))  

# UMAP plot
ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  ggrastr::geom_point_rast(size = 0.01,
                           alpha = 0.5,
                           shape = 16) +
  scale_color_manual(values = custom_colors) +
  coord_fixed() +
  theme_classic(base_size = 9) +
  theme(legend.position = "none", panel.grid = element_blank()) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  ggrepel::geom_text_repel(
    data = plot_df %>%
      group_by(cell_type) %>%
      summarize(
        x = median(UMAP_1, na.rm = TRUE),
        y = median(UMAP_2, na.rm = TRUE),
        .groups = "drop"
      ),
    aes(x = x, y = y, label = cell_type),
    color = "black",
    # <-- sets label text to black
    size = 2.5,
    segment.color = "grey40",
    segment.size = 0.1,
    box.padding = 0.3,
    max.overlaps = Inf
  )

# Save plot
ggsave(
  filename = "Output_files/SMI/Figures/Figure2E1.pdf",
  width = 4,
  height = 4,
  dpi = 600
)

#----------------------------------------------------------------------------------------------------
# UMAP with 16S Signal
#----------------------------------------------------------------------------------------------------

# Add 16S norm signal and intracellular column
plot_df <- plot_df %>%
  mutate(
    signal = glioma_slim@meta.data$Microbial16S_bg_norm,
    is_highlighted = Cell_location_Microbial16S_bg == "Intracellular_Microbial16S" # High confidence 16S+ cells
  )

# Plot UMAP
ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = signal)) +
  ggrastr::geom_point_rast(
    data = filter(plot_df, !is_highlighted),
    color = "lightgrey",
    alpha = 0.5,
    size = 0.01,
    dpi = 600,
    shape = 16
  ) +
  ggrastr::geom_point_rast(
    data = filter(plot_df, is_highlighted),
    aes(color = signal),
    alpha = 0.8,
    size = 0.8,
    dpi = 600,
    shape = 16
  ) +
  scale_color_viridis_c(
    name = "16S signal",
    option = "D",
    breaks = c(11, 50, 100, 500),
    labels = c("11", "50", "100", "500"),
    guide = guide_legend(
      override.aes = list(size = 4, shape = 15),
      title.position = "top",
      label.theme = element_text(size = 11),
      title.theme = element_text( size = 14),
      direction = "vertical"
    )
  ) +
  coord_fixed() +
  theme_classic(base_size = 9) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")

# Save plot
ggsave(
  filename = "Output_files/SMI/Figures/Figure2F1.pdf",
  width = 4,
  height = 4,
  dpi = 600
)
