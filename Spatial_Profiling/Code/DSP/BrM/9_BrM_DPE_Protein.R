#----------------------------------------------------------------------------------------------------
# Code based on GeoMXTools Protein vignette (nanoString)
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# Source files
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")
source("src/Functions.R")
load("Structure/DSP/BrM/BrM_ProteinDatanorm.RData")

#----------------------------------------------------------------------------------------------------
# Transform data
#----------------------------------------------------------------------------------------------------

# Log samples
assayDataElement(object = BrM_ProteinData, elt = "log_q_neg_norm") <-
    assayDataApply(BrM_ProteinData, 2, FUN = log, base = 2, elt = "neg_norm")

#---------------------------------------------------------------------------------------------------
# Save object
#---------------------------------------------------------------------------------------------------
save(BrM_ProteinData, file = "Structure/DSP/BrM/BrM_ProteinData_PCA.RData")

#---------------------------------------------------------------------------------------------------
# Differential protein expression (Linear Mix Model)
#---------------------------------------------------------------------------------------------------

# 25-75% (Tumor only)
results <- c()

for (compartment in "BrM") {
    ind <- pData(BrM_ProteinData)$Type == compartment
    mixedOutmc <-
        mixedModelDE(BrM_ProteinData[, ind],
            elt = "log_q_neg_norm",
            modelFormula = ~ test_16SrRNA_status_25_75 + Primary_Tumor +
                (1 + test_16SrRNA_status_25_75 | Patient), # Controlling for Primary tumor and 16S random slope
            groupVar = "test_16SrRNA_status_25_75",
            nCores = parallel::detectCores(),
            multiCore = FALSE
        )
    r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
    tests <- rownames(r_test)
    r_test <- as.data.frame(r_test)
    r_test$Contrast <- tests
    r_test$Protein <-
        unlist(lapply(
            colnames(mixedOutmc),
            rep, nrow(mixedOutmc["lsmeans", ][[1]])
        ))
    r_test$Subset <- compartment
    r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "BH")
    r_test <- r_test[, c(
        "Protein", "Subset", "Contrast", "Estimate",
        "Pr(>|t|)", "FDR"
    )]
    results <- rbind(results, r_test)
}
BrM_DPE_by_16S_25_75 <- as.data.frame(results) %>%
    dplyr::rename(log2FC = Estimate) # rename for better readability

# Save object
save(BrM_DPE_by_16S_25_75, file = "Structure/DSP/BrM/BrM_DPE_by_16S_25_75.RData")

# Write csv
write.csv(BrM_DPE_by_16S_25_75, "Output_files/DSP/BrM/BrM_DPE_by_16S_25_75.csv")

#---------------------------------------------------------------------------------------------------
# Protein volcano plot
#---------------------------------------------------------------------------------------------------

# Plot volcano using add_volcano_protein function
add_volcano_protein(BrM_DPE_by_16S_25_75, theme_size = 12) +
    labs(x = expression("16S rRNA-Low <-- log"[2] ~ "(Fold Change) --> 16S rRNA-High")) +
    geom_text_repel(
        data = subset(BrM_DPE_by_16S_25_75, abs(log2FC) > 0.58 & FDR < 0.05),
        size = 1.8, max.overlaps = 20, nudge_y = 0.15, segment.size = 0.15,
        color = "black", fontface = "bold", alpha = 0.6
    )

# Save pdf
ggsave("Figures/FigureS4D.pdf", width = 6, height = 4)

#----------------------------------------------------------------------------------------------------
# Protein dotplot                                                                                                   
#----------------------------------------------------------------------------------------------------

# Select biologically relevant proteins
proteins <- c(
    "IRF3", "NFkB p105 / p50", "IKK gamma/NEMO", "TLR9", "NF-kB p65", "IKB alpha",
    "TRAF2", "IRF5", "MyD88", "CD276", "Neutrophil Elastase", "LAG-3", "HLA G", "IFNGR1",
    "IL-12A", "CD94", "Lactate Dehydrogenase", "FABP4", "IDH1", "ENO1 + ENO2 + ENO3",
    "Hexokinase II", "Aldolase", "Caspase-3 p12", "DR5", "ATG7"
)

# Transform data for plotting
BrM_DPE_by_16S_25_75_dotplot <- BrM_DPE_by_16S_25_75 %>%
    filter(FDR < 0.05 & Protein %in% proteins) %>%
    dplyr::rename("Target" = "Protein") %>%
    arrange(factor(Target, levels = proteins)) %>%
    mutate(
        Analyte = rep("Protein", 25),
        Category = (rep(c(
            "Anti-microbial response", "Immune response", "Metabolism and lipid homeostasis",
            "Stress, apoptosis, autophagy"
        ), times = c(9, 7, 6, 3)))
    ) %>%
    group_by(Category) %>%
    arrange(log2FC, .by_group = TRUE) %>%
    ungroup() %>%
    mutate(Order = rep(1:25))

# Plot dotplot
ggplot(BrM_DPE_by_16S_25_75_dotplot, aes(
    x = reorder(Target, Order, decreasing = TRUE),
    y = Analyte, color = FDR, size = log2FC
)) +
    geom_point() +
    scale_color_distiller(palette = "OrRd", guide = guide_colorbar(
        barwidth = 5
    )) +
    theme_bw(base_size = 9) +
    theme(
        legend.key.size = unit(0.3, "cm"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank(), strip.placement = "outside",
        strip.background = element_rect(colour = "white", fill = "white"),
        axis.text.x = element_text(angle = 45, hjust = 1), strip.text = element_text(face = "bold"),
        legend.margin = margin(), legend.position = "bottom"
    ) +
    facet_grid(~Category,
        scales = "free", space = "free_x",
        labeller = labeller(Category = label_wrap_gen(25))
    ) +
    scale_size_continuous(limits = c(0.5, 2.5)) +
    labs(x = " ", y = " ", size = "log2FC", color = "FDR")

# Save pdf
ggsave("Figures/Figure3D.pdf", width = 8, height = 2.25)

