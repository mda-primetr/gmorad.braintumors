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

# Define gene sets
gene_sets <- c(
  "MicrobeHost_Interaction_AUC",
  "TumorImmuneRecruitment_AUC",
  "TumorLipid_Metabolism_AUC",
  "stress_apoptosis_genes_AUC",
  "TumorAngiogenesis_AUC"
)

#Defince clean labels
gene_labels <- c(
  "Microbe–Host interactions",
  "Immune recruitment",
  "Lipid metabolism",
  "Stress response",
  "Angiogenesis"
)

# Set names
names(gene_labels) <- gene_sets

# Colors for neighborhoods
neigh_colors <- c("16S- NH" = "#218C24", "16S+ NH" = "#FF00FD")

# Build data for plotting
plot_data <- brm_slim@meta.data %>%
  dplyr::select(Patient, neighborhood_intra_label, all_of(gene_sets)) %>%
  filter(neighborhood_intra_label %in% c("16S- Neighborhood", "16S+ Neighborhood"), Patient %in% c("P41", "P43", "P46")) %>%
  pivot_longer(
    cols = all_of(gene_sets),
    names_to = "gene_set",
    values_to = "auc_score"
  ) %>%
  filter(!is.na(auc_score)) %>%
  mutate(
    gene_label = gene_labels[gene_set],
    Patient = factor(Patient),
    neighborhood_intra_label = factor(neighborhood_intra_label, 
                                      levels = c("16S- Neighborhood", "16S+ Neighborhood"),
                                      labels = c("16S- NH", "16S+ NH"))
  )

# Force order for facets
plot_data$gene_label <- factor(
  plot_data$gene_label,
  levels = c("Microbe–Host interactions", "Immune recruitment", "Lipid metabolism", "Stress response", "Angiogenesis")
)

#----------------------------------------------------------------------------------------------------
# Violin plot with per-patient comparison
#----------------------------------------------------------------------------------------------------

p <- ggplot(plot_data, aes(x = neighborhood_intra_label, y = auc_score, fill = neighborhood_intra_label)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE, alpha = 0.6, width = 0.8) +
  geom_boxplot(width = 0.2, outlier.shape = NA, color = "black", size = 0.2, fill = "white") +
  facet_grid(rows = vars(Patient), cols = vars(gene_label), scales = "free",
             labeller = labeller(gene_label = label_wrap_gen(10))) +
  stat_compare_means(
    aes(group = neighborhood_intra_label),
    method = "wilcox.test",
    label = "p.format",
    label.y.npc = "top",
    label.x.npc = "center",
    size = 2.8,
    digits = 3
  ) +
  scale_fill_manual(values = neigh_colors) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
      strip.background = element_blank(),  
    panel.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 9.5),
    axis.title = element_blank(),
    axis.text.x = element_text(size = 8),
    axis.ticks.x = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 5))

# Save plot
ggsave(
  filename = "Output_files/SMI/Figures/FigureS5C.pdf",
  plot = p,
   width = 5.5,
  height = 3.5
)
