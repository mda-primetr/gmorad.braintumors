#----------------------------------------------------------------------------------------------------
# Code based on GeoMXTools vignette by Jason Reeves, Prajan Divakar, Nicole Ortogero,
# Maddy Griswold, Zhi Yang, Stephanie Zimmerman, Rona Vitancol and David Henderson (nanoString)
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# Source Files
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")
load("Structure/DSP/BrM/BrM_RNAData_raw.RData")
load("Structure/DSP/BrM/BrM_Modules.RData")

#----------------------------------------------------------------------------------------------------
# Pre-processing
#----------------------------------------------------------------------------------------------------

# Shift counts to one
BrM_RNAData <- shiftCountsOne(BrM_RNAData, useDALogic = TRUE)

#----------------------------------------------------------------------------------------------------
# ROI QC
#----------------------------------------------------------------------------------------------------

# QC parameters
QC_params <-
  list(
    minSegmentReads = 1000,
    percentTrimmed = 80,
    percentStitched = 80,
    percentAligned = 80,
    percentSaturation = 50,
    minNegativeCount = 1,
    maxNTCCount = 1000,
    minNuclei = 50,
    minArea = 1000
  )

# Set ROI QC flags
BrM_RNAData <- setSegmentQCFlags(BrM_RNAData, qcCutoffs = QC_params)

# QC Results
QCResults <- protocolData(BrM_RNAData)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(
  Pass = colSums(!QCResults[, flag_columns]),
  Warning = colSums(QCResults[, flag_columns])
)
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})
QC_Summary["TOTAL FLAGS", ] <- c(
  sum(QCResults[, "QCStatus"] == "PASS"),
  sum(QCResults[, "QCStatus"] == "WARNING")
)

# Print QC summary
QC_Summary %>%
  kbl() %>%
  kable_minimal(full_width = F)

# Show all NTC values
table(NTC_Count = sData(BrM_RNAData)$NTC)

# Remove ROIs that failed QC
BrM_RNAData <- BrM_RNAData[, QCResults$QCStatus == "PASS"]

#----------------------------------------------------------------------------------------------------
# Probe QC
#----------------------------------------------------------------------------------------------------

# Set probe QC flags
BrM_RNAData <- setBioProbeQCFlags(BrM_RNAData,
  qcCutoffs = list(
    minProbeRatio = 0.1,
    percentFailGrubbs = 20
  ),
  removeLocalOutliers = TRUE
)

ProbeQCResults <- fData(BrM_RNAData)[["QCFlags"]]

# Exclude outlier probes
ProbeQCPassed <- subset(
  BrM_RNAData,
  fData(BrM_RNAData)[["QCFlags"]][, c("LowProbeRatio")] == FALSE &
    fData(BrM_RNAData)[["QCFlags"]][, c("GlobalGrubbsOutlier")] == FALSE
)

BrM_RNAData <- ProbeQCPassed

# Remove Bugs-in-space (BIS) probes not present in both GeoMx and CosMx studies (Total 16S probes)
BIS_probes_to_be_removed <- c("RTS1003672", "RTS1003675", "RTS1003676", "RTS1003680")

BrM_RNAData <- subset(BrM_RNAData, !fData(BrM_RNAData)$RTS_ID %in% BIS_probes_to_be_removed)

#----------------------------------------------------------------------------------------------------
# Collapse probes
#----------------------------------------------------------------------------------------------------

# Create gene-level count data
BrM_target_Data <- aggregateCounts(BrM_RNAData)

#----------------------------------------------------------------------------------------------------
# LOQ
#----------------------------------------------------------------------------------------------------

# Define LOQ SD threshold and minimum value
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module (WTA and BIS)
LOQ <- data.frame(row.names = colnames(BrM_target_Data))

for (module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"), module)
  if (all(vars[1:2] %in% colnames(pData(BrM_target_Data)))) {
    LOQ[, module] <- pmax(
      minLOQ,
      pData(BrM_target_Data)[, vars[1]] *
        pData(BrM_target_Data)[, vars[2]]^cutoff
    )
  }
}

# Add LOQ to pData
pData(BrM_target_Data)$LOQ <- LOQ

# Determine genes > LOQ per module (WTA and BIS)
LOQ_Mat <- c()
for (module in modules) {
  ind <- fData(BrM_target_Data)$Module == module
  Mat_i <- t(esApply(BrM_target_Data[ind, ],
    MARGIN = 1,
    FUN = function(x) {
      x > LOQ[, module]
    }
  ))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}

LOQ_Mat <- LOQ_Mat[fData(BrM_target_Data)$TargetName, ]

#----------------------------------------------------------------------------------------------------
# Probe QC
#----------------------------------------------------------------------------------------------------

# Calculate gene detection rate per ROI
pData(BrM_target_Data) <- pData(BrM_target_Data) %>%
  mutate(
    GenesDetected = colSums(LOQ_Mat, na.rm = TRUE),
    GeneDetectionRate = GenesDetected / nrow(BrM_target_Data)
  )

# Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(BrM_target_Data)$DetectionThreshold <-
  cut(pData(BrM_target_Data)$GeneDetectionRate,
    breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
    labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%")
  )

# Plot gene detection per ROI
ggplot(
  pData(BrM_target_Data),
  aes(x = DetectionThreshold)
) +
  geom_bar(aes(fill = Type), alpha = 0.7, width = 0.5) +
  geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5, size = 2) +
  scale_fill_npg() +
  theme_bw(base_size = 12) +
  theme(
    legend.key.size = unit(0.3, "cm"), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "Gene detection rate",
    y = "ROI, #",
    fill = "Type"
  )

# Remove ROIs with < 10% gene detection rate
BrM_target_Data <- BrM_target_Data[, pData(BrM_target_Data)$GeneDetectionRate >= 0.1]

# Calculate gene detection across ROIs
LOQ_Mat <- LOQ_Mat[, colnames(BrM_target_Data)]
fData(BrM_target_Data) <- fData(BrM_target_Data) %>%
  mutate(
    DetectedSegments = rowSums(LOQ_Mat, na.rm = TRUE),
    DetectionRate = DetectedSegments / nrow(pData(BrM_target_Data))
  )

# Determine thresholds: 1%, 5%, 10%, 20%, 30% and 50%
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))

# Number of probes detected
plot_detect$Number <-
  unlist(lapply(
    c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
    function(x) {
      sum(fData(BrM_target_Data)$DetectionRate >= x)
    }
  ))

# Percent of probes detected
plot_detect$Rate <- plot_detect$Number / nrow(fData(BrM_target_Data))
rownames(plot_detect) <- plot_detect$Freq

# Plot percent of probes detected across ROIs
ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
    vjust = 1.6, color = "black", size = 2.5
  ) +
  scale_fill_gradient2(
    low = "orange2", mid = "lightblue",
    high = "dodgerblue3", midpoint = 0.65,
    limits = c(0, 1),
    labels = scales::percent
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.key.size = unit(0.3, "cm"), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(
    labels = scales::percent, limits = c(0, 1),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(x = "% of ROIs", y = str_wrap("Genes Detected, % of Panel > LOQ",
    width = 16
  ))

# Get negative probes
negativeProbefData <- subset(fData(BrM_target_Data), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)

# Subset to genes detected in 10% of ROIs. Manually add negative probes
BrM_target_Data <- BrM_target_Data[fData(BrM_target_Data)$DetectionRate >= 0.1 |
  fData(BrM_target_Data)$TargetName %in% neg_probes, ]

#----------------------------------------------------------------------------------------------------
# Subset object by probe set
#----------------------------------------------------------------------------------------------------

# Subset WTA counts
BrM_target_Data_WTA <- subset(BrM_target_Data, Module == "Hs_R_NGS_WTA_v1.0")
save(BrM_target_Data_WTA, file = "Structure/DSP/BrM/BrM_target_Data_WTA.RData")

# Subset BIS counts
BrM_target_Data_BIS <- subset(BrM_target_Data, Module == "BiS_pilot_1_RnD_v1.2")
save(BrM_target_Data_BIS, file = "Structure/DSP/BrM/BrM_target_Data_BIS.RData")
