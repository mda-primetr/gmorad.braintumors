#----------------------------------------------------------------------------------------------------
# Code based on GeoMXTools Protein vignette (nanoString)
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# Source files
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")
load("Structure/DSP/Glioma/Glioma_ProteinData.RData")

#----------------------------------------------------------------------------------------------------
# Choose best normalization method
#----------------------------------------------------------------------------------------------------

# Get housekeeping proteins
hk.names <- c("GAPDH", "Histone H3", "RPS6")

# Get IgGs
igg.names <- as.data.frame(fData(Glioma_ProteinData)) %>%
    filter(CodeClass == "Negative") %>%
    dplyr::select(TargetName) %>%
    deframe()

# Plot IgG concordance
plotConcordance(object = Glioma_ProteinData, targetList = igg.names, plotFactor = "test_16SrRNA_status_25_75")

# Plot HK concordance
plotConcordance(object = Glioma_ProteinData, targetList = hk.names, plotFactor = "test_16SrRNA_status_25_75")

# Plot norm factor concordance
normfactors <- computeNormalizationFactors(
    object = Glioma_ProteinData,
    igg.names = igg.names,
    hk.names = hk.names,
    area = "area",
    nuclei = "nuclei"
)

# Plot norm factor concordance
plotNormFactorConcordance(
    object = Glioma_ProteinData, plotFactor = "test_16SrRNA_status_25_75",
    normfactors = normfactors
)

# Given the strong correlation between IgGs and the weak correlation between HKs, we will
# perform background normalization using the IgGs.

#----------------------------------------------------------------------------------------------------
# Background normalization
#----------------------------------------------------------------------------------------------------

# Subset object to contain IgGs
iggsubset <- subset(Glioma_ProteinData, subset = TargetName %in% igg.names)

# Calculate geoMean of neg
iggs <- apply(exprs(iggsubset), 2, function(x) ngeoMean(x))

# Calculate neg normalization factor and save in pData
pData(Glioma_ProteinData)$iggFactors_exprs <- iggs / ngeoMean(iggs)

# Normalize to neg
assayDataElement(Glioma_ProteinData, "neg_norm") <-
    sweep(assayDataElement(Glioma_ProteinData, "exprs"), 2L, iggs / ngeoMean(iggs),
        FUN = "/"
    )

#----------------------------------------------------------------------------------------------------
# Save object
#----------------------------------------------------------------------------------------------------

save(Glioma_ProteinData, file = "Structure/DSP/Glioma/Glioma_ProteinDatanorm.RData")
