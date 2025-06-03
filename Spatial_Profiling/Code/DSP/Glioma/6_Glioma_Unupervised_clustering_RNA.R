#----------------------------------------------------------------------------------------------------
# Source files                                                                                                   
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")
source("src/Functions.R")
load("Processed_data/Glioma/Glioma_target_Data_WTA_PCA.RData")

# Set reproducibility seed
set.seed(42)

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
# Filter out NA values for test_16SrRNA_status_25_75
#----------------------------------------------------------------------------------------------------

df_data_test <- df_data %>%
    filter(!is.na(test_16SrRNA_status_25_75))

#----------------------------------------------------------------------------------------------------
# Run PERMANOVA (adonis2) 
#----------------------------------------------------------------------------------------------------

adonis_res <- adonis2(df_data_test[, c("PC1", "PC2")] ~ test_16SrRNA_status_25_75,
    data = df_data_test,
    method = "euclidean"
)

# Format stats label
p_val <- adonis_res$`Pr(>F)`[1]
r2_val <- adonis_res$R2[1]
adonis_label <- paste0(
    "R² = ", round(r2_val, 2),
    ", p = ", format.pval(p_val, digits = 3, eps = 0.001)
)

#----------------------------------------------------------------------------------------------------
# Plot PCA 
#----------------------------------------------------------------------------------------------------

ggplot(df_data_test, aes(x = PC1, y = PC2, color = test_16SrRNA_status_25_75)) +
    geom_point(alpha = 0.6, size = 2) +
    stat_ellipse(type = "norm", linetype = "dashed") +
    scale_color_manual(values = c("darkblue", "darkred")) +
    theme_bw(base_size = 12) +
    theme(
        legend.key.size = unit(0.4, "cm"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    labs(x = "PC1", y = "PC2", color = "16S category") +
    annotate("text",
        x = Inf, y = -Inf, hjust = 1.1, vjust = -0.5,
        label = adonis_label, size = 2
    )

# Save pdf
ggsave("Output_files/DSP/Figures/Figure4B1.pdf", width = 4, height = 2.5)
