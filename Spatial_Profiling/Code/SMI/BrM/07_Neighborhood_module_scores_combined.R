#----------------------------------------------------------------------------------------------------
# Source files
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")

# Set reproducibility seed
set.seed(42)

# Load slimmed Seurat object with metadata and AUCell scores
brm_slim <- readRDS("Input_files/SMI/BrM/brm_slim_seurat.rds") 

#----------------------------------------------------------------------------------------------------
# Pre-processing
#----------------------------------------------------------------------------------------------------
# Patients to exclude due to low neighborhood cell counts 
remove_patient_nh <- c("P33", "P39", "P47", "P40", "P32", "P35")

# Define gene sets
gene_sets <- c(
  "MicrobeHost_Interaction_AUC",
  "TumorImmuneRecruitment_AUC",
  "TumorLipid_Metabolism_AUC",
  "stress_apoptosis_genes_AUC",
  "TumorAngiogenesis_AUC"
)

# Define clean labels
gene_labels <- c(
  "Microbe–Host interactions",
  "Immune recruitment",
  "Lipid metabolism",
  "Stress response",
  "Angiogenesis"
)

# Set names
names(gene_labels) <- gene_sets

# Neighborhood group colors
neigh_colors <- c("16S- NH" = "#218C24", "16S+ NH" = "#FF00FD")

# Build combined metadata for all gene sets
plot_data <- brm_slim@meta.data %>%
  dplyr::select(Patient, neighborhood_intra_label, all_of(gene_sets)) %>%
  filter(
    neighborhood_intra_label %in% c("16S+ Neighborhood", "16S- Neighborhood"),
    !Patient %in% remove_patient_nh
  ) %>%
  pivot_longer(cols = all_of(gene_sets),
               names_to = "gene_set",
               values_to = "auc_score") %>%
  filter(!is.na(auc_score)) %>%
  mutate(
    gene_label = gene_labels[gene_set],
    neighborhood_intra_label = factor(
      neighborhood_intra_label,
      levels = c("16S- Neighborhood", "16S+ Neighborhood"),
      labels = c("16S- NH", "16S+ NH")
    )
  )

# Set plotting order of facet
plot_data$gene_label <- factor(
  plot_data$gene_label,
  levels = c(
    "Microbe–Host interactions",
    "Immune recruitment",
    "Lipid metabolism",
    "Stress response",
    "Angiogenesis"
  )
)

#----------------------------------------------------------------------------------------------------
# Violin plot comparing 16S– vs. 16S+ neighborhoods
#----------------------------------------------------------------------------------------------------

p <- ggplot(
  plot_data,
  aes(x = neighborhood_intra_label, y = auc_score, fill = neighborhood_intra_label)
) +
  geom_violin(
    scale = "width",
    adjust = 1,
    trim = TRUE,
    alpha = 0.6,
    width = 0.8
  ) +
  geom_boxplot(
    width = 0.2,
    outlier.shape = NA,
    color = "black",
    size = 0.2,
    fill = "white"
  ) +
  facet_wrap(
    ~ gene_label,
    scales = "free_y",
    nrow = 1,
    labeller = labeller(gene_label = label_wrap_gen(width = 15))
  ) +
  scale_fill_manual(values = neigh_colors) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 10, face = "bold"),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5, ),
    axis.ticks.x = element_blank()
  ) +
  scale_x_discrete(
    labels = function(x)
      str_wrap(x, width = 5)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

# Save plot
ggsave(
  filename = "Output_files/SMI/Figures/Figure5F.pdf",
  plot = p,
  width = 8,
  height = 2.25
)
