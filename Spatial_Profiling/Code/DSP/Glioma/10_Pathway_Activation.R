#----------------------------------------------------------------------------------------------------
# Source files
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")
load("Processed_data/DSP/Glioma/Glioma_ProteinData_PCA.RData")

# Set reproducibility seed
set.seed(42)

#----------------------------------------------------------------------------------------------------
# Define signature markers
#----------------------------------------------------------------------------------------------------

signature_proteins <- list(
    bacterial_response = c(
        "HMGB1", "HMGB2"
    ),
    metabolic_response = c(
        "GLB1/Beta-galactosidase", "Fatty Acid Synthase"
    ),
    general_inflammation = c(
        "CD68", "HLA-DR", "CD44", "IL-6R", "MMP9", "C Reactive Protein", "STAT3" 
    )
)

#----------------------------------------------------------------------------------------------------
# Compute signature scores
#----------------------------------------------------------------------------------------------------

# Create empty list to store scores
signature_score_list <- list()

# Loop over each signature
for (sig_name in names(signature_proteins)) {
    sig_proteins <- signature_proteins[[sig_name]]
    sig_present <- sig_proteins[sig_proteins %in% rownames(Glioma_ProteinData)]

    # Subset to only available proteins in the data
    sig_data <- Glioma_ProteinData[rownames(Glioma_ProteinData) %in% sig_present, ]

    # Calculate geometric mean (on original scale)
    score <- assayDataApply(sig_data, MARGIN = 2, FUN = ngeoMean, elt = "neg_norm")

    # Store log2-transformed scores
    signature_score_list[[sig_name]] <- log2(score)
}

# Combine all into one data frame
score_matrix <- do.call(cbind, signature_score_list)

score_matrix <- as.matrix(t(score_matrix))

#----------------------------------------------------------------------------------------------------
# Create object
#----------------------------------------------------------------------------------------------------

# Create feature_data slot
feature_data <- data.frame(Signature = names(signature_proteins))
rownames(feature_data) <- names(signature_proteins)

signature_set <- NanoStringGeoMxSet(
    assayData = score_matrix,
    phenoData = AnnotatedDataFrame(pData(Glioma_ProteinData)),
    protocolData = protocolData(Glioma_ProteinData),
    featureData = AnnotatedDataFrame(feature_data),
    featureType = "Signature",
    check = FALSE
)

#----------------------------------------------------------------------------------------------------
# Run Linear Mix Model
#----------------------------------------------------------------------------------------------------

model_results <- mixedModelDE(
    signature_set,
    elt = "exprs",
    modelFormula = ~ test_16SrRNA_status_25_75 + IDH + Recurrence + (1 | Patient),
    groupVar = "test_16SrRNA_status_25_75",
    nCores = parallel::detectCores(),
    multiCore = FALSE
)

#----------------------------------------------------------------------------------------------------
# Extract and format results
#----------------------------------------------------------------------------------------------------

results_df <- map_dfr(names(signature_proteins), function(sig) {
    lsm <- model_results["lsmeans", ][[sig]]
    data.frame(
        signature = sig,
        estimate = lsm[1, "Estimate"],
        p_value = lsm[1, "Pr(>|t|)"]
    )
}) %>%
    mutate(
        fdr = p.adjust(p_value, method = "BH"),
        inflammation_ratio = 2^(estimate - estimate[signature == "general_inflammation"]),
        signature_type = ifelse(signature == "general_inflammation", "Reference signature", "16S-associated signature"),
        sig_label = case_when(
            fdr < 0.001 ~ "***",
            fdr < 0.01 ~ "**",
            fdr < 0.05 ~ "*",
            TRUE ~ "ns"
        )
    )

#----------------------------------------------------------------------------------------------------
# Plot
#----------------------------------------------------------------------------------------------------

x_max <- max(results_df$estimate)
x_range <- x_max
padding <- x_range * 0.1

ggplot(results_df, aes(x = estimate, y = reorder(signature, estimate))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(size = -log10(fdr), color = signature_type)) +
    geom_text(aes(x = estimate + 0.05, label = sig_label), hjust = -0.5, size = 4) +
    scale_color_manual(
        values = c("16S-associated signature" = "#203E71", "Reference signature" = "#CC3D3C"),
        name = "Signature type"
    ) +
    scale_size_continuous(name = "-log10(FDR)", range = c(3, 8)) +
    scale_x_continuous(
        limits = c(0 - padding, x_max + padding + 0.1)
    ) +
    scale_y_discrete(labels = c(
        "Tumor\ninflammation",
        "Anti-microbial",
        "Metabolic"
    )) +
    labs(
        x = expression("16S-High vs 16S-Low (log"[2] * " fold change)"),
        y = NULL,
    ) +
    theme_bw(base_size = 11) +
    theme(
        legend.key.size = unit(0.4, "cm"),
        panel.grid = element_blank(),
        plot.caption = element_text(size = 8, hjust = 0),
        legend.position = "right"
    ) +
    guides(
        color = guide_legend(order = 1),
        size = guide_legend(order = 2)
    )

# Save pdf
ggsave("Output_files/DSP/Figures/FigureS4G1.pdf", width = 5, height = 3)
