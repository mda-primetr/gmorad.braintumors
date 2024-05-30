#----------------------------------------------------------------------------------------------------
# Code based on GeoMXTools Protein vignette (nanoString)
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# Source files
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")
source("src/Functions.R")
load("Structure/DSP/BrM/BrM_ProteinData_raw.RData")
load("Structure/DSP/BrM/BrM_target_Data_WTA_16S_25_75.RData")

#----------------------------------------------------------------------------------------------------
# Match RNA ROIs
#----------------------------------------------------------------------------------------------------

BrM_ProteinData <- BrM_ProteinData[, which(pData(BrM_ProteinData)$RNA_ID %in%
     pData(BrM_target_Data_WTA_16S_25_75)$RNA_ID)]

# Create DPE factors
pData(BrM_ProteinData) <- pData(BrM_ProteinData) %>%
     mutate(
          Patient = factor(pData(BrM_target_Data_WTA_16S_25_75)$Patient),
          test_16SrRNA_status_25_75 = factor(pData(BrM_target_Data_WTA_16S_25_75)$`16SrRNA_status_25_75`,
               levels = c("High", "Low")
          ),
          Primary_Tumor = factor(pData(BrM_target_Data_WTA_16S_25_75)$Primary_Tumor)
     )

#----------------------------------------------------------------------------------------------------
# Target QC
#----------------------------------------------------------------------------------------------------

# Get housekeeping proteins
hk.names <- c("GAPDH", "Histone H3", "RPS6")

# Get IgGs
igg.names <- as.data.frame(fData(BrM_ProteinData)) %>%
     filter(CodeClass == "Negative") %>%
     dplyr::select(TargetName) %>%
     deframe()

# Select proteins that have at least 10% ROIs with Signal-Noise ratio > 2
SNR <- as.data.frame(snrOrder(BrM_ProteinData, igg.names)) %>%
     rownames_to_column() %>%
     dplyr::rename("Protein" = "rowname") %>%
     pivot_longer(cols = starts_with("DSP")) %>%
     mutate(SNR = case_when(
          value > 2 ~ "Above 2",
          value < 2 ~ "Below 2"
     )) %>%
     filter(SNR == "Above 2" | Protein %in% igg.names | Protein %in% hk.names) %>%
     group_by(Protein) %>%
     filter(n() > 4 | Protein %in% igg.names | Protein %in% hk.names) %>%
     select(Protein) %>%
     distinct()

# Remove low performing proteins
BrM_ProteinData <- BrM_ProteinData[SNR$Protein, ]

# Plot remaining proteins
fig <- qcProteinSignal(object = BrM_ProteinData, neg.names = igg.names)
fig()
rect(
     xleft = 0, xright = 5,
     ybottom = -2, ytop = 2, density = 0, col = "darkblue", lwd = 1
)
rect(
     xleft = 5, xright = 338,
     ybottom = 2, ytop = 8.9, density = 0, col = "darkred", lwd = 1
)

#----------------------------------------------------------------------------------------------------
# Save object
#----------------------------------------------------------------------------------------------------

save(BrM_ProteinData, file = "Structure/DSP/BrM/BrM_ProteinData.RData")
