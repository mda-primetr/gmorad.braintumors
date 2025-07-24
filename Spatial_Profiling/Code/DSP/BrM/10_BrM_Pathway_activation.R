#----------------------------------------------------------------------------------------------------
# Source files
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")
load("Processed_data/DSP/BrM/BrM_ProteinData_PCA.RData")

# Set reproducibility seed
set.seed(42)

#----------------------------------------------------------------------------------------------------
# Define signature markers
#----------------------------------------------------------------------------------------------------

signature_proteins <- list(
    bacterial_response = c(
        "TLR9", "MyD88", "IRF3", "IRF5",
        "NFkB p105 / p50", "NF-kB p65",
        "IKK gamma/NEMO", "IKB alpha"
    ),
    metabolic_response = c(
        "Hexokinase II", "Aldolase", "IDH1",
        "ENO1 + ENO2 + ENO3", "Lactate Dehydrogenase"
    ),
    general_inflammation = c(
        "CD68", "CD14", "CD44", "MMP9", "C Reactive Protein", "STAT3"
    )
)

#----------------------------------------------------------------------------------------------------
# Compute log2 geometric mean scores for each signature
#----------------------------------------------------------------------------------------------------

signature_score_list <- list()

for (sig_name in names(signature_proteins)) {
    # Filter proteins available in the dataset
    present_proteins <- signature_proteins[[sig_name]][signature_proteins[[sig_name]] %in% rownames(BrM_ProteinData)]

    # Subset to available markers
    sig_data <- BrM_ProteinData[present_proteins, ]

    # Compute geometric mean (original scale) across all markers per sample
    score <- assayDataApply(sig_data, MARGIN = 2, FUN = ngeoMean, elt = "neg_norm")

    # Log2-transform and store scores
    signature_score_list[[sig_name]] <- log2(score)
}

# Combine into matrix with samples as rows and signatures as columns
score_matrix <- t(do.call(cbind, signature_score_list))

#----------------------------------------------------------------------------------------------------
# Wrap as NanoStringGeoMxSet for LMM
#----------------------------------------------------------------------------------------------------

# Create feature data for signature labels
feature_data <- data.frame(Signature = names(signature_proteins))
rownames(feature_data) <- names(signature_proteins)

# Assemble NanoStringGeoMxSet
signature_set <- NanoStringGeoMxSet(
    assayData = score_matrix,
    phenoData = AnnotatedDataFrame(pData(BrM_ProteinData)),
    protocolData = protocolData(BrM_ProteinData),
    featureData = AnnotatedDataFrame(feature_data),
    featureType = "Signature",
    check = FALSE
)

#----------------------------------------------------------------------------------------------------
# Run linear mixed model
#----------------------------------------------------------------------------------------------------

model_results <- mixedModelDE(
    signature_set,
    elt = "exprs",
    modelFormula = ~ test_16SrRNA_status_25_75 + (1 | Patient),
    groupVar = "test_16SrRNA_status_25_75",
    nCores = parallel::detectCores(),
    multiCore = FALSE
)

#----------------------------------------------------------------------------------------------------
# Format results
#----------------------------------------------------------------------------------------------------

results_df <- map_dfr(names(signature_proteins), function(sig) {
    lsm <- model_results["lsmeans", ][[sig]]
    data.frame(
        signature = sig,
        estimate = lsm[1, "Estimate"], # log2 fold-change
        p_value = lsm[1, "Pr(>|t|)"] # p-value
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
# Plot signature results
#----------------------------------------------------------------------------------------------------

x_max <- max(results_df$estimate)
padding <- x_max * 0.1

ggplot(results_df, aes(x = estimate, y = reorder(signature, estimate))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(aes(size = -log10(fdr), color = signature_type)) +
    geom_text(aes(x = estimate + 0.05, label = sig_label), hjust = -0.5, size = 4) +
    scale_color_manual(
        values = c("16S-associated signature" = "#203E71", "Reference signature" = "#CC3D3C"),
        name = "Signature type"
    ) +
    scale_size_continuous(name = "-log10(FDR)", range = c(3, 8)) +
    scale_x_continuous(limits = c(0 - padding, x_max + padding + 0.1)) +
    scale_y_discrete(labels = c(
        "general_inflammation" = "Tumor\ninflammation",
        "bacterial_response" = "Anti-microbial",
        "metabolic_response" = "Metabolic"
    )) +
    labs(
        x = expression("16S-High vs. 16S-Low (log"[2] * " fold change)"),
        y = NULL
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
ggsave("Output_files/DSP/Figures/FigureS4G2.pdf", width = 5, height = 3)
