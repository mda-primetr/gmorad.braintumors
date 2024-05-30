#----------------------------------------------------------------------------------------------------
# Source files                                                                                                   
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")
source("src/Functions.R")
load("Structure/DSP/Glioma/Glioma_target_Data_WTA_PCA.RData")

#----------------------------------------------------------------------------------------------------
# Load data and transform                                                                                                   
#----------------------------------------------------------------------------------------------------

df_data <- as.data.frame(t(assayDataElement(Glioma_target_Data_WTA, elt = "log_q_bgsub"))) %>%
    dplyr::select(-"NegProbe-WTX")

#----------------------------------------------------------------------------------------------------
# Run prcomp                                                                                                   
#----------------------------------------------------------------------------------------------------

data_pca <- prcomp(df_data, center = TRUE, scale. = TRUE)

# PCA variance
summary(data_pca)
# 0.1997  0.09557

# Add factors
df_data <- df_data %>%
    mutate(
        PC1 = data_pca$x[, 1], PC2 = data_pca$x[, 2],
        Patient = factor(pData(Glioma_target_Data_WTA)$Patient),
        test_16SrRNA_status_25_75 = factor(pData(Glioma_target_Data_WTA)$test_16SrRNA_status_25_75,
            levels = c("Low", "High")
        ),
        IDH = factor(pData(Glioma_target_Data_WTA)$IDH),
        Recurrence = factor(pData(Glioma_target_Data_WTA)$Recurrence)
    )

#----------------------------------------------------------------------------------------------------
# Plot by 16S category                                                                                                   
#----------------------------------------------------------------------------------------------------

ggplot(subset(df_data, test_16SrRNA_status_25_75 != "NA"), aes(
    x = PC1, y = PC2,
    color = test_16SrRNA_status_25_75
)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_manual(values = c("darkblue", "darkred")) +
    theme_bw(base_size = 12) +
    theme(
        legend.key.size = unit(0.4, "cm"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    stat_ellipse() +
    labs(x = "PC1", y = "PC2", color = "16S category")

# Save pdf
ggsave("Figures/Figure3B1.pdf", width = 4, height = 2.5)
