#----------------------------------------------------------------------------------------------------
# Source files
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")

# Set reproducibility seed
set.seed(42)

#----------------------------------------------------------------------------------------------------
# Pre-processing
#----------------------------------------------------------------------------------------------------

# Load DGE results across all patients (P11, P15, P26) and compute additional stats
glioma_dge_neighborhood_all <- read_csv("Input_files/SMI/Glioma/glioma_dge_neighborhood.csv") %>%
  dplyr::select(contrast, term, target, ncells_1, ncells_2, fold_change, p.value) %>%
  mutate(
    Patient = "ALL",                          # Label all rows with a unified patient tag
    log2FC = log2(fold_change),               # Compute log2 fold change
    padj = p.adjust(p.value, method = "BH")   # Apply Benjamini–Hochberg correction
  )

# Load curated signature gene list
signature_genes_df <- read_csv("Input_files/SMI/Glioma/glioma_neighborhood_gene_list.csv")

# Filter DGE to include only signature genes with significant upregulation
glioma_dge_neighborhood_all_sig <- glioma_dge_neighborhood_all %>%
  inner_join(signature_genes_df, by = c("target")) %>%
  filter(p.value < 0.05 & log2FC > 0.58) %>%
  arrange(Category, Patient) %>%
  mutate(Order = rep(1:29))

# Set order
glioma_dge_neighborhood_all_sig$Category <- factor(
  glioma_dge_neighborhood_all_sig$Category,
  levels = c(
    "Anti-microbial and immune signatures",
    "Metabolism and lipid homeostasis",
    "Cell stress",
    "Angiogenesis"
  )
)

#----------------------------------------------------------------------------------------------------
# Dot plot
#----------------------------------------------------------------------------------------------------

# Create dot plot of signature genes

ggplot(
  glioma_dge_neighborhood_all_sig,
  aes(
    x = Patient,
    y = reorder(target, Order),
    fill = p.value,
    size = log2FC
  )
) +
  geom_point(shape = 21,
             # Shape 21 allows both fill and color
             color = "grey30",
             # This is the border color
             stroke = 0.3) +
  # Color scale (for fill)
  scale_fill_distiller(palette = "PuBu",
                       direction = 1,
                       trans = "reverse") +
  
  scale_size_continuous(limits = c(0.58, 4)) +
  facet_grid(~ Category, scales = "free", space = "free_x") +
  coord_flip() +
  theme_bw(base_size = 9) +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.6, "cm"),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.placement = "outside",
    strip.background = element_rect(colour = "white", fill = "white"),
    strip.text = element_text(face = "bold")
  ) +
  labs(x = NULL,
       y = NULL,
       size = "log2FC",
       fill = "p-value") +
  guides(fill = guide_colorbar(order = 1), size = guide_legend(order = 2))

# Save plot
ggsave(
  filename = "Output_files/SMI/Figures/Figure5E1.pdf",
  width = 10,
  height = 1.8
)
