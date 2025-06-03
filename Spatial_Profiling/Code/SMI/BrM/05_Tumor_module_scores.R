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

# Patients to exclude (those with < 20 16S+ tumor cells)
remove_patient <- c("P33", "P39", "P47", "P40", "P32")

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

# Tumor group colors
tumor_colors <- c("16S- Tumor" = "darkblue",
                  "16S+ Tumor" = "darkred")

# Build combined metadata for all gene sets
plot_data <- brm_slim@meta.data %>%
  dplyr::select(Patient, DE_tumor_bg, all_of(gene_sets)) %>%
  filter(DE_tumor_bg %in% c("16Sneg_tumor", "16Spos_tumor"),
         !Patient %in% remove_patient) %>% # Filter patients to those with > 20 16S+ cells
  pivot_longer(cols = all_of(gene_sets),
               names_to = "gene_set",
               values_to = "auc_score") %>%
  filter(!is.na(auc_score)) %>%
  mutate(gene_label = gene_labels[gene_set],
         DE_tumor_bg = factor(DE_tumor_bg, levels = c("16Sneg_tumor", "16Spos_tumor")))

# Convert DE_tumor_bg to readable factors
plot_data$DE_tumor_bg <- factor(
  plot_data$DE_tumor_bg,
  levels = c("16Sneg_tumor", "16Spos_tumor"),
  labels = c("16S- Tumor", "16S+ Tumor")
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
# Violin plot
#---------------------------------------------------------------------------------------------------- 

p <- ggplot(plot_data, aes(x = DE_tumor_bg, y = auc_score, fill = DE_tumor_bg)) +
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
    labeller = labeller(gene_label = label_wrap_gen(width = 10))
  ) +
  scale_fill_manual(values = tumor_colors) +
  theme_bw(base_size = 10) +
  theme(
    panel.grid = element_blank(),
    legend.position = "none",
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 10, face = "bold"),
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    axis.ticks.x = element_blank()
  ) +
  scale_x_discrete(
    labels = function(x)
      str_wrap(x, width = 5)
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs()

# Save plot
ggsave(
  filename = "Output_files/SMI/Figures/Figure5C.pdf"),
  plot = p,
  width = 8,
  height = 2.25)
