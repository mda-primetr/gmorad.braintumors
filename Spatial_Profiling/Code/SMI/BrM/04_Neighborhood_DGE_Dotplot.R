#----------------------------------------------------------------------------------------------------
# Source files
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")

# Set reproducibility seed
set.seed(42)

#----------------------------------------------------------------------------------------------------
# Pre-processing
#----------------------------------------------------------------------------------------------------

# Load DGE results across all patients (P41, P43, P46) and compute additional stats
brm_dge_neighborhood_all <- read_csv("Input_files/SMI/BrM/brm_dge_neighborhood_P41_P43_P46.csv") %>%
  dplyr::select(
    contrast, target, ncells_1, ncells_2, fold_change, p.value
  ) %>%
  mutate(
    Patient = "ALL",                         # Label all rows with a unified patient tag
    log2FC = log2(fold_change),              # Compute log2 fold change
    padj = p.adjust(p.value, method = "BH")  # Apply Benjamini–Hochberg correction
  )


#----------------------------------------------------------------------------------------------------
# Define gene categories of interest
#----------------------------------------------------------------------------------------------------

gene_categories <- c(
  # Anti-microbial and immune response genes
  "CHI3L1" = "antimicrobial_immune_response",
  "IL22RA2" = "antimicrobial_immune_response",

  # Metabolism and lipid homeostasis
  "ACSL1" = "metabolism_lipid_homeostasis",
  "ACSL5" = "metabolism_lipid_homeostasis",
  "FADS2" = "metabolism_lipid_homeostasis",
  "FASN" = "metabolism_lipid_homeostasis",
  "FABP3" = "metabolism_lipid_homeostasis",
  "FABP1" = "metabolism_lipid_homeostasis",
  "SCD" = "metabolism_lipid_homeostasis",
  "SREBF1" = "metabolism_lipid_homeostasis",
  "LDLR" = "metabolism_lipid_homeostasis",
  "LDLRAP1" = "metabolism_lipid_homeostasis",
  "G6PD" = "metabolism_lipid_homeostasis",
  "PLIN5" = "metabolism_lipid_homeostasis",
  "THBS1" = "metabolism_lipid_homeostasis",
  "IRS1" = "metabolism_lipid_homeostasis",
  "IRS2" = "metabolism_lipid_homeostasis",
  "PCSK9" = "metabolism_lipid_homeostasis",
  "ABCA2" = "metabolism_lipid_homeostasis",
  "PREP" = "metabolism_lipid_homeostasis",
  "ETFDH" = "metabolism_lipid_homeostasis",
  "CPT1A" = "metabolism_lipid_homeostasis",
  "XBP1" = "metabolism_lipid_homeostasis",
  "ABCD1" = "metabolism_lipid_homeostasis",
  "APOE" = "metabolism_lipid_homeostasis",
  "ADIPOR1" = "metabolism_lipid_homeostasis",
  "CPNE3" = "metabolism_lipid_homeostasis",

  # Stress response, apoptosis, and autophagy
  "ATF4" = "stress_apoptosis_autophagy",
  "ATF6" = "stress_apoptosis_autophagy",
  "XBP1" = "stress_apoptosis_autophagy",
  "ERN1" = "stress_apoptosis_autophagy",
  "DDIT3" = "stress_apoptosis_autophagy",
  "HSPA1A" = "stress_apoptosis_autophagy",
  "HSPB1" = "stress_apoptosis_autophagy",
  "P4HB" = "stress_apoptosis_autophagy",
  "NFE2L2" = "stress_apoptosis_autophagy",
  "GPX1" = "stress_apoptosis_autophagy",
  "SOD1" = "stress_apoptosis_autophagy",
  "SOD2" = "stress_apoptosis_autophagy",
  "PRDX1" = "stress_apoptosis_autophagy",
  "SESN2" = "stress_apoptosis_autophagy",
  "TXNDC2" = "stress_apoptosis_autophagy",
  "OXR1" = "stress_apoptosis_autophagy",
  "DNAJB11" = "stress_apoptosis_autophagy",
  "HSPA9" = "stress_apoptosis_autophagy",
  "TP53" = "stress_apoptosis_autophagy",
  "ATM" = "stress_apoptosis_autophagy",
  "CHEK1" = "stress_apoptosis_autophagy",
  "CHEK2" = "stress_apoptosis_autophagy",
  "GADD45A" = "stress_apoptosis_autophagy",
  "GADD45B" = "stress_apoptosis_autophagy",
  "CASP9" = "stress_apoptosis_autophagy",
  "CASP8" = "stress_apoptosis_autophagy",
  "PMAIP1" = "stress_apoptosis_autophagy",
  "BCL2" = "stress_apoptosis_autophagy",
  "BAK1" = "stress_apoptosis_autophagy",
  "MCL1" = "stress_apoptosis_autophagy",
  "BAX" = "stress_apoptosis_autophagy",
  "BID" = "stress_apoptosis_autophagy",
  "FADD" = "stress_apoptosis_autophagy",
  "DIABLO" = "stress_apoptosis_autophagy",
  "EIF4E" = "stress_apoptosis_autophagy",
  "EIF2AK3" = "stress_apoptosis_autophagy",
  "DUSP6" = "stress_apoptosis_autophagy",
  "ATG4D" = "stress_apoptosis_autophagy",
  "MAP1LC3B" = "stress_apoptosis_autophagy",
  "CUL2" = "stress_apoptosis_autophagy",
  "FBXO9" = "stress_apoptosis_autophagy",
  "UBA5" = "stress_apoptosis_autophagy",
  "UBE2G1" = "stress_apoptosis_autophagy",

  # Angiogenesis
  "BMP7" = "angiogenesis",
  "CD164" = "angiogenesis",
  "EPHA7" = "angiogenesis",
  "FGFR2" = "angiogenesis",
  "PTPRF" = "angiogenesis",
  "TMSB10" = "angiogenesis",
  "THBS1" = "angiogenesis",
  "VEGFA" = "angiogenesis",
  "MDK" = "angiogenesis",
  "LGMN" = "angiogenesis",
  "PRKD1" = "angiogenesis",
  "PRKD2" = "angiogenesis"
)

# Set clean labels
clean_labels <- c(
  "antimicrobial_immune_response" = "Anti–microbial and immune response",
  "metabolism_lipid_homeostasis" = "Metabolism and lipid homeostasis",
  "stress_apoptosis_autophagy" = "Stress response, apoptosis/autophagy",
  "angiogenesis" = "Angiogenesis"
)

# Factor level ordering for categories
category_levels <- c(
  "Anti–microbial and immune response",
 "Metabolism and lipid homeostasis",
 "Stress response, apoptosis/autophagy",
  "Angiogenesis"
)

#----------------------------------------------------------------------------------------------------
# Filter significant genes and annotate
#----------------------------------------------------------------------------------------------------

# Keep genes with nominal p < 0.05 and log2FC > 0.58
brm_dge_neighborhood_all_sig <- brm_dge_neighborhood_all %>%
  filter(p.value < 0.05 & log2FC > 0.58) %>% 
  mutate(Order = rep(1:49))

# Map gene → category_key → clean label
brm_dge_neighborhood_all_sig$category_key <- gene_categories[brm_dge_neighborhood_all_sig$target]
brm_dge_neighborhood_all_sig$category <- clean_labels[brm_dge_neighborhood_all_sig$category_key]
brm_dge_neighborhood_all_sig$category <-  factor(brm_dge_neighborhood_all_sig$category, levels = category_levels)

brm_dge_neighborhood_all_sig <- brm_dge_neighborhood_all_sig %>%
  filter( !is.na(category))

#----------------------------------------------------------------------------------------------------
# Dot plot
#----------------------------------------------------------------------------------------------------

ggplot(brm_dge_neighborhood_all_sig,
       aes(
         x = Patient,
         y = reorder(target, Order),
         fill = padj,
         size = log2FC
       )) +
  geom_point(shape = 21, color = "grey30", stroke = 0.3) +
  scale_fill_distiller(palette = "PuBu", direction = 1, trans = "reverse")+
  scale_size_continuous(limits = c(0.58, 4)) +
  facet_grid(~ category, scales = "free", space = "free_x",labeller = labeller(category = label_wrap_gen(20))) +
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
  labs(x = NULL, y = NULL, size = "log2FC", fill = "FDR") +
  guides(
    fill = guide_colorbar(order = 1),
    size = guide_legend(order = 2)
  )

# Save plot
  ggsave(
    filename = "Output_files/SMI/Figures/Figure5E2.pdf"),
    width = 10,
    height = 1.8
  )
