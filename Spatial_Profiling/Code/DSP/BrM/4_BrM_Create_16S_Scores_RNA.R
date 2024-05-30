#----------------------------------------------------------------------------------------------------
# Source files
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")
load("Structure/DSP/BrM/BrM_target_Data_BIS_norm.RData")
load("Structure/DSP/BrM/BrM_target_Data_WTA_norm.RData")

#----------------------------------------------------------------------------------------------------
# Create 16S score
#----------------------------------------------------------------------------------------------------

# Remove custom negative set 1
BIS_probes_to_be_removed <- c("Custom Negative Set 1")
BrM_target_Data_BIS <- subset(
    BrM_target_Data_BIS,
    !fData(BrM_target_Data_BIS)$TargetName %in% BIS_probes_to_be_removed
)

# Calculate geoMean of 16S probes (score)
BIS_probe_geoMean <- data.frame(assayDataApply(BrM_target_Data_BIS,
    MARGIN = 2, FUN = ngeoMean, elt = "neg_norm_bgsub"
))
colnames(BIS_probe_geoMean) <- "BIS_probe_geoMean"

# Save score in pData of WTA and BIS
pData(BrM_target_Data_BIS) <- pData(BrM_target_Data_BIS) %>%
    mutate(`16SrRNA_score` = BIS_probe_geoMean$BIS_probe_geoMean)
pData(BrM_target_Data_WTA) <- pData(BrM_target_Data_WTA) %>%
    mutate(`16SrRNA_score` = BIS_probe_geoMean$BIS_probe_geoMean)

#----------------------------------------------------------------------------------------------------
# Bacterial scores by tissue type
#----------------------------------------------------------------------------------------------------

# Create df with 16S scores and phenoData
df_assay_BIS <- pData(BrM_target_Data_BIS) %>%
    mutate(Type = factor(Type, levels = c("BrM", "Tumor border")))

# Plot scores
ggplot(df_assay_BIS, aes(x = Type, y = `16SrRNA_score`, fill = Type)) +
    geom_boxplot(alpha = 0.5, show.legend = FALSE, width = 0.4) +
    geom_dotplot(
        dotsize = 1, binaxis = "y", stackdir = "center",
        position = position_dodge(0.75), show.legend = FALSE
    ) +
    scale_fill_manual(values = c("darkred", "darkblue")) +
    theme_bw(base_size = 12) +
    theme(
        legend.key.size = unit(0.4, "cm"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    labs(x = "", y = "16S rRNA score") +
    scale_x_discrete(, labels = function(x) str_wrap(x, width = 5))

# Save pdf
ggsave("Figures/FigureS4A2.pdf", width = 1.85, height = 2.2)

# Add ROI # to BIS
pData(BrM_target_Data_BIS) <- pData(BrM_target_Data_BIS) %>%
    mutate(ROI = protocolData(BrM_target_Data_BIS)$roi)

# Add ROI # to WTA
pData(BrM_target_Data_WTA) <- pData(BrM_target_Data_WTA) %>%
    mutate(ROI = protocolData(BrM_target_Data_WTA)$roi)

#----------------------------------------------------------------------------------------------------
# Save objects
#----------------------------------------------------------------------------------------------------

save(BrM_target_Data_BIS, file = "Structure/DSP/BrM/BrM_target_Data_BIS_BIS_scores.RData")
save(BrM_target_Data_WTA, file = "Structure/DSP/BrM/BrM_target_Data_WTA_BIS_scores.RData")
