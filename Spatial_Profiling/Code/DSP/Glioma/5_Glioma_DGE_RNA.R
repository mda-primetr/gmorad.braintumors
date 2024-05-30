#----------------------------------------------------------------------------------------------------
# Code based on GeoMXTools vignette by Jason Reeves, Prajan Divakar, Nicole Ortogero,
# Maddy Griswold, Zhi Yang, Stephanie Zimmerman, Rona Vitancol and David Henderson (nanoString)
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# Source files                                                                                                   
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")
load("Structure/DSP/Glioma/Glioma_target_Data_WTA_BIS_scores.RData")

#----------------------------------------------------------------------------------------------------
# Transform data                                                                                                   
#----------------------------------------------------------------------------------------------------

# Log samples
assayDataElement(object = Glioma_target_Data_WTA, elt = "log_q_bgsub") <-
    assayDataApply(Glioma_target_Data_WTA, 2, FUN = log, base = 2, elt = "q_norm_bgsub_1")

# Remove adjacent normal
Glioma_target_Data_WTA <- Glioma_target_Data_WTA[, which(Glioma_target_Data_WTA$Type
!= "Adjacent normal")]

#----------------------------------------------------------------------------------------------------
# Create factors                                                                                                   
#----------------------------------------------------------------------------------------------------

# Create 16S categories
pData(Glioma_target_Data_WTA) <- pData(Glioma_target_Data_WTA) %>%
    mutate(`16SrRNA_status_25_75` = case_when(
        `16SrRNA_score` < quantile(`16SrRNA_score`, 0.25) ~ "Low",
        `16SrRNA_score` > quantile(`16SrRNA_score`, 0.75) ~ "High"
    ))

# Create DGE factors
pData(Glioma_target_Data_WTA) <- pData(Glioma_target_Data_WTA) %>%
    mutate(
        Patient = factor(pData(Glioma_target_Data_WTA)$Patient),
        test_16SrRNA_status_25_75 = factor(pData(Glioma_target_Data_WTA)$`16SrRNA_status_25_75`,
            levels = c("High", "Low")
        ),
        IDH = factor(pData(Glioma_target_Data_WTA)$IDH),
        Recurrence = factor(pData(Glioma_target_Data_WTA)$Recurrence)
    )

#----------------------------------------------------------------------------------------------------
# Subset                                                                                                   
#----------------------------------------------------------------------------------------------------

# Subset ROIs that have no NAs in test_16SrRNA_status_25_75
Glioma_target_Data_WTA_16S_25_75 <- Glioma_target_Data_WTA[
    , which(pData(Glioma_target_Data_WTA)$test_16SrRNA_status_25_75 != "NA")
]

#----------------------------------------------------------------------------------------------------
# Save objects                                                                                                   
#----------------------------------------------------------------------------------------------------

# All tumor ROIs
save(Glioma_target_Data_WTA, file = "Structure/DSP/Glioma/Glioma_target_Data_WTA_PCA.RData")

# Only those w/no NAs in test_16SrRNA_status_25_75
save(Glioma_target_Data_WTA_16S_25_75,
    file = "Structure/DSP/Glioma/Glioma_target_Data_WTA_16S_25_75.RData"
)

#----------------------------------------------------------------------------------------------------
# Bacterial plot by High vs Low                                                                                                   
#----------------------------------------------------------------------------------------------------

# Create df with scores and 16S categories
df_BIS_scores <- as.data.frame(pData(Glioma_target_Data_WTA)) %>%
    mutate(test_16SrRNA_status_25_75 = factor(test_16SrRNA_status_25_75, levels = c("Low", "High"))) %>%
    filter(test_16SrRNA_status_25_75 != "NA")

# Plot 16S scores by category
ggplot(df_BIS_scores, aes(
    x = test_16SrRNA_status_25_75, y = `16SrRNA_score`,
    fill = test_16SrRNA_status_25_75
)) +
    geom_boxplot(alpha = 0.5, show.legend = FALSE, width = 0.3) +
    geom_dotplot(
        dotsize = 1, binaxis = "y", stackdir = "center",
        position = position_dodge(0.75), show.legend = FALSE
    ) +
    scale_fill_manual(values = c("darkblue", "darkred")) +
    theme_bw(base_size = 12) +
    theme(
        legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    labs(x = "", y = "16S rRNA score")

# Save pdf
ggsave("Figures/FigureS4B1.pdf", width = 1.85, height = 2.2)

#----------------------------------------------------------------------------------------------------
# Differential gene expression (Linear Mix Model)
#----------------------------------------------------------------------------------------------------

# 25-75% (Tumor only)
results <- c()

for (compartment in "Glioma") {
    ind <- pData(Glioma_target_Data_WTA_16S_25_75)$Type == compartment
    mixedOutmc <-
        mixedModelDE(Glioma_target_Data_WTA_16S_25_75[, ind],
            elt = "log_q_bgsub",
            modelFormula = ~ test_16SrRNA_status_25_75 + IDH + Recurrence +
                (1 | Patient), # Controlling for IDH and Recurrence
            groupVar = "test_16SrRNA_status_25_75",
            nCores = parallel::detectCores(),
            multiCore = FALSE
        )
    r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
    tests <- rownames(r_test)
    r_test <- as.data.frame(r_test)
    r_test$Contrast <- tests
    r_test$Gene <-
        unlist(lapply(
            colnames(mixedOutmc),
            rep, nrow(mixedOutmc["lsmeans", ][[1]])
        ))
    r_test$Subset <- compartment
    r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "BH")
    r_test <- r_test[, c(
        "Gene", "Subset", "Contrast", "Estimate",
        "Pr(>|t|)", "FDR"
    )]
    results <- rbind(results, r_test)
}
Glioma_DGE_by_16S_25_75 <- as.data.frame(results) %>%
    dplyr::rename(log2FC = Estimate) # rename for better readability

# Save object
save(Glioma_DGE_by_16S_25_75, file = "Structure/DSP/Glioma/Glioma_DGE_by_16S_25_75.RData")
# Write csv
write.csv(Glioma_DGE_by_16S_25_75, "Output_files/DSP/Glioma/Glioma_DGE_by_16S_25_75.csv")

#----------------------------------------------------------------------------------------------------
# Dot plots                                                                                                   
#----------------------------------------------------------------------------------------------------

# Select biologically relevant genes
genes <- c(
    "ECSIT", "IL32", "CADM1", "GSDME", "FKBIL1", "TRIM5", "TRIM24", "TRAF3IP1", "TRIM26",
    "TNFAIP1", "SP110", "IL11RA", "CBLB", "IFI35", "PILRB", "TRADD", "TAP1", "STAT5B", "TNFRSF10B", 
    "IL17RB", "NKIRAS2", "LRBA"
)

# Transform data for plotting
Glioma_DGE_by_16S_25_75_dotplot <- Glioma_DGE_by_16S_25_75 %>%
    filter(FDR < 0.05 & Gene %in% genes) %>%
    dplyr::rename("Target" = "Gene") %>%
    arrange(factor(Target, levels = genes)) %>%
    mutate(
        Analyte = rep("RNA", 21),
        Order = (rep(1:21))
    )

# Plot genes
ggplot(Glioma_DGE_by_16S_25_75_dotplot, aes(
    x = Analyte, y = reorder(Target, Order, decreasing = TRUE),
    color = FDR, size = log2FC
)) +
    geom_point() +
    scale_color_distiller(palette = "PuBu") +
    theme_bw(base_size = 12) +
    theme(
        legend.key.size = unit(0.3, "cm"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x = element_blank()
    ) +
    labs(x = " ", y = " ", size = "log2FC", color = "FDR") +
    guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2))

# Save pdf
ggsave("Figures/FigureS4E.pdf", width = 2.8, height = 6)
