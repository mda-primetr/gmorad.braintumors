#----------------------------------------------------------------------------------------------------
# Code based on GeoMXTools vignette by Jason Reeves, Prajan Divakar, Nicole Ortogero,
# Maddy Griswold, Zhi Yang, Stephanie Zimmerman, Rona Vitancol and David Henderson (nanoString)
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# Soirce files                                                                                                   
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")
load("Structure/DSP/BrM/BrM_target_Data_WTA.RData")
load("Structure/DSP/BrM/BrM_target_Data_BIS.RData")

#----------------------------------------------------------------------------------------------------
# Signal:Background                                                                                                   
#----------------------------------------------------------------------------------------------------

# Calculate Q3 and background values
ann_of_interest <- "Type"
Stat_data <-
    data.frame(
        row.names = colnames(exprs(BrM_target_Data_WTA)),
        ROI = colnames(exprs(BrM_target_Data_WTA)),
        Annotation = pData(BrM_target_Data_WTA)[, ann_of_interest],
        Q3 = unlist(apply(exprs(BrM_target_Data_WTA), 2, # Q3 values
            quantile, 0.75,
            na.rm = TRUE
        )),
        NegProbe = exprs(BrM_target_Data_WTA)["NegProbe-WTX", ]
    ) # Background values

# Transform data for plotting
Stat_data_m <- Stat_data %>%
    pivot_longer(cols = c("Q3", "NegProbe"), names_to = "Statistic", values_to = "Value")

# Plot Q3 and background histograms
ggplot(
    Stat_data_m,
    aes(x = Value, fill = Statistic)
) +
    geom_histogram(bins = 40, alpha = 0.7) +
    theme_bw(base_size = 12) +
    theme(
        legend.key.size = unit(0.3, "cm"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    scale_x_continuous(trans = "log2") +
    scale_fill_brewer(palette = 3, type = "qual") +
    labs(x = "Counts", y = "ROI, #") +
    facet_wrap(~Annotation, ncol = 2)

# Plot Q3 against Background
ggplot(
    Stat_data,
    aes(x = NegProbe, y = Q3, color = Annotation)
) +
    geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
    geom_point(alpha = 0.7) +
    theme_bw(base_size = 12) +
    theme(
        legend.key.size = unit(0.3, "cm"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    scale_color_manual(values = c("darkblue", "darkred")) +
    scale_x_continuous(trans = "log2") +
    scale_y_continuous(trans = "log2") +
    labs(x = "Negative Probe GeoMean", y = "Q3 Value")

# Plot Signal:Background
ggplot(Stat_data, aes(x = Q3 / NegProbe, fill = Annotation)) +
    geom_histogram(bins = 50, alpha = 0.7) +
    theme_bw(base_size = 12) +
    theme(
        legend.key.size = unit(0.3, "cm"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    scale_fill_npg() +
    facet_wrap(~Annotation, ncol = 2) +
    labs(x = "Signal:Background", y = "ROI, #", fill = "Type")

# Plot Signal:Background against Background
ggplot(
    Stat_data,
    aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)
) +
    geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
    geom_point(alpha = 0.7) +
    theme_bw(base_size = 12) +
    theme(
        legend.key.size = unit(0.3, "cm"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    scale_color_npg() +
    scale_x_continuous(trans = "log2") +
    scale_y_continuous(trans = "log2") +
    labs(
        x = "Negative Probe GeoMean", y = "Q3/NegProbe",
        color = "Type"
    )

#----------------------------------------------------------------------------------------------------
# WTA probe normalization
#----------------------------------------------------------------------------------------------------

# Subtract background to account for signal:noise discrepancies
BrM_target_Data_WTA <- normalize(BrM_target_Data_WTA,
    norm_method = "subtractBackground",
    fromElt = "exprs",
    toElt = "bgsub"
)

# Q3 normalization (background subtracted counts)
BrM_target_Data_WTA <- normalize(BrM_target_Data_WTA,
    norm_method = "quant",
    desiredQuantile = .75,
    fromElt = "bgsub",
    toElt = "q_norm_bgsub"
)

# Threshold data to max value <1 for log2 transformation
# Create new assay slot for Q3 norm with background subtraction + threshold
assayDataElement(object = BrM_target_Data_WTA, elt = "q_norm_bgsub_1") <- assayDataElement(
    object =
        BrM_target_Data_WTA, elt = "q_norm_bgsub"
)

# Find max value <1 for threshold
max_q_norm_bgsub <- max(assayDataElement(BrM_target_Data_WTA,
    elt =
        "q_norm_bgsub"
)[which(assayDataElement(BrM_target_Data_WTA, "q_norm_bgsub") < 1)])

# Make counts <1 the max value <1
assayDataElement(BrM_target_Data_WTA, "q_norm_bgsub_1")[assayDataElement(
    BrM_target_Data_WTA,
    "q_norm_bgsub_1"
) < 1] <- max_q_norm_bgsub

#----------------------------------------------------------------------------------------------------
# BIS probe normalization
#----------------------------------------------------------------------------------------------------

# Subtract background to account for signal:noise discrepancies
BrM_target_Data_BIS <- normalize(BrM_target_Data_BIS,
    norm_method = "subtractBackground",
    fromElt = "exprs",
    toElt = "bgsub"
)
# Calculate background normalization factor
pData(BrM_target_Data_BIS)$neg_normFactor <-
    pData(BrM_target_Data_BIS)$NegGeoMean_BiS_pilot_1_RnD_v1.2 /
        ngeoMean(pData(BrM_target_Data_BIS)$NegGeoMean_BiS_pilot_1_RnD_v1.2)

# Normalize to background
assayDataElement(BrM_target_Data_BIS, "neg_norm_bgsub") <-
    sweep(assayDataElement(BrM_target_Data_BIS, "bgsub"), 2L,
        pData(BrM_target_Data_BIS)$neg_normFactor,
        FUN = "/"
    )

# Threshold data to max value <1 for log2 transformation
# Create new assay slot for Q3 norm with background subtraction + threshold
assayDataElement(object = BrM_target_Data_BIS, elt = "neg_norm_bgsub_1") <-
    assayDataElement(object = BrM_target_Data_BIS, elt = "neg_norm_bgsub")

# Find max <1 value for threshold
max_neg_norm_bgsub <- max(assayDataElement(BrM_target_Data_BIS,
    elt = "neg_norm_bgsub"
)[which(assayDataElement(BrM_target_Data_BIS, "neg_norm_bgsub") < 1)])

# Make counts <1 the max value <1
assayDataElement(BrM_target_Data_BIS, "neg_norm_bgsub_1")[assayDataElement(
    BrM_target_Data_BIS,
    "neg_norm_bgsub_1"
) < 1] <- max_neg_norm_bgsub

#----------------------------------------------------------------------------------------------------
# Save objects
#----------------------------------------------------------------------------------------------------

save(BrM_target_Data_WTA, file = "Structure/DSP/BrM/BrM_target_Data_WTA_norm.RData")
save(BrM_target_Data_BIS, file = "Structure/DSP/BrM/BrM_target_Data_BIS_norm.RData")
