#----------------------------------------------------------------------------------------------------
# Source files
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")

# Set reproducibility seed
set.seed(42)

# Load slim seurat object
brm_slim <- readRDS("Input_files/SMI/BrM/brm_slim_seurat.rds") 

#----------------------------------------------------------------------------------------------------
# UMAP with cell types
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

# Assign colors to cell types
custom_colors <- c(
  "T46" = "#BDB76B",
  "T43" = "#A5AA99",
  "T41" = "#708090",
  "T47" = "#CD5C5C",
  "T32" = "#8FBC8F",
  "T39" = "#D2691E",
  "T40" = "#4682B4",
  "Mac/Mic" = "#f683a2",
  "T35" = "#006400",
  "MSC-like" = "#3b267e",
  "T33" = "#5D69B1",
  "EC" = "#6a2013",
  "Mural" = "#f15167",
  "Plasma" = "#2F8AC4",
  "Astroglial" = "#DAA51B",
  "DC" = "#6B7FD7",
  "Neutrophil" = "#24796C",
  "NK" = "#CC61B0",
  "Erythrocyte" = "#99C945",
  "CD4 T" = "#40E0D0",
  "Treg" = "#E58606",
  "CD8 T" = "#764E9F"
)

# Get embeddings
umap_df <- brm_slim@reductions$umap@cell.embeddings %>%
  as.data.frame() %>% tibble::rownames_to_column("cell_ID")

# Get relevant metadata
meta_df <- brm_slim@meta.data

# Join embeddings and metadata
plot_df <- left_join(meta_df, umap_df, by = "cell_ID") %>%
  mutate(cell_type = factor(cell_type_umap))

# UMAP plot
ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2, color = cell_type)) +
  ggrastr::geom_point_rast(size = 0.01, 
                           # To reduce image size
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
    size = 2.5,
    segment.color = "grey40",
    segment.size = 0.1,
    box.padding = 0.3,
    max.overlaps = Inf
  )

# Save plot
ggsave(
  filename = "Output_files/SMI/Figures/Figure2E2.pdf",
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
    signal = brm_slim@meta.data$Microbial16S_bg_norm,
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
      title.theme = element_text(size = 14),
      direction = "vertical"
    )
  ) +
  coord_fixed() +
  theme_classic(base_size = 9) +
  theme(legend.position = "right") +
  labs(x = "UMAP 1", y = "UMAP 2")

# Save plot
ggsave(
  filename = "Output_files/SMI/Figures/Figure2F2.pdf",
  width = 4,
  height = 4,
  dpi = 600
)
