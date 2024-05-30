#----------------------------------------------------------------------------------------------------
# Code based on GeoMXTools Protein vignette (nanoString)
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# Source files                                                                                                   
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")
load("Structure/DSP/BrM/BrM_ProteinData.RData")

#----------------------------------------------------------------------------------------------------
# Choose best normalization method                                                                                                   
#----------------------------------------------------------------------------------------------------

# Get housekeeping proteins
hk.names <- c("GAPDH", "Histone H3", "RPS6")

# Get IgGs
igg.names <- as.data.frame(fData(BrM_ProteinData)) %>%
    filter(CodeClass == "Negative") %>%
    dplyr::select(TargetName) %>%
    deframe()

# Plot IgG concordance
plotConcordance(object = BrM_ProteinData, targetList = igg.names, plotFactor = "test_16SrRNA_status_25_75")

# Plot HK concordance
plotConcordance(object = BrM_ProteinData, targetList = hk.names, plotFactor = "test_16SrRNA_status_25_75")

# Plot norm factor concordance
normfactors <- computeNormalizationFactors(
    object = BrM_ProteinData,
    igg.names = igg.names,
    hk.names = hk.names,
    area = "area",
    nuclei = "nuclei"
)

# Plot norm factor concordance
plotNormFactorConcordance(
    object = BrM_ProteinData, plotFactor = "test_16SrRNA_status_25_75",
    normfactors = normfactors
)

# Given the strong correlation between IgGs and the weak correlation between HKs, we will
# perform background normalization using the IgGs.

#----------------------------------------------------------------------------------------------------
# Backgroun normalization                                                                                                   
#----------------------------------------------------------------------------------------------------

# Subset object to contain IgGs
iggsubset <- subset(BrM_ProteinData, subset = TargetName %in% igg.names)

# Calculate geoMean of neg
iggs <- apply(exprs(iggsubset), 2, function(x) ngeoMean(x))

# Calculate neg normalization factor and save in pData
pData(BrM_ProteinData)$iggFactors_exprs <- iggs / ngeoMean(iggs)

# Normalize to background
assayDataElement(BrM_ProteinData, "neg_norm") <-
    sweep(assayDataElement(BrM_ProteinData, "exprs"), 2L, iggs / ngeoMean(iggs),
        FUN = "/"
    )

#----------------------------------------------------------------------------------------------------
# Save object                                                                                                   
#----------------------------------------------------------------------------------------------------

save(BrM_ProteinData, file = "Structure/DSP/BrM/BrM_ProteinDatanorm.RData")
