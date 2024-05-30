#----------------------------------------------------------------------------------------------------
# Code based on GeoMXTools vignette by Jason Reeves, Prajan Divakar, Nicole Ortogero,
# Maddy Griswold, Zhi Yang, Stephanie Zimmerman, Rona Vitancol and David Henderson (nanoString)
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# Source files
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")

# Set working directory
datadir <- file.path("Input_files/DSP/BrM")

#----------------------------------------------------------------------------------------------------
# RNA NanoStringGeoMxSet object
#----------------------------------------------------------------------------------------------------

# Files location
DCCFiles <- dir(file.path(datadir, "DCC"), pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
PKCFiles <- dir(file.path(datadir, "PKC"), pattern = ".pkc$", full.names = TRUE, recursive = TRUE)
AnnotationFile <- file.path(datadir, "BrM_Annotations.xlsx")

# Create RNA object
BrM_RNAData <- suppressWarnings(readNanoStringGeoMxSet(
    dccFiles = DCCFiles,
    pkcFiles = PKCFiles,
    phenoDataFile = AnnotationFile,
    phenoDataSheet = "Annotations_RNA",
    phenoDataDccColName = "Sample_ID",
    protocolDataColNames = c("aoi", "roi"),
    configFile = NULL,
    analyte = "RNA",
    phenoDataColPrefix = "",
    experimentDataColNames = NULL
))

# Modules used
pkcs <- annotation(BrM_RNAData)
modules <- gsub(".pkc", "", pkcs)

#----------------------------------------------------------------------------------------------------
# Protein NanoStringGeoMxSet object
#----------------------------------------------------------------------------------------------------

# Load feature data
BrM_IPA_featureData <- read_csv("Input_files/DSP/BrM/BrM_IPA_featureData.csv") %>%
    arrange(ProbeName) %>% # sort by ProbeName
    column_to_rownames("Target_Name") %>% # rename rows with protein name
    as.matrix() # convert to matrix

# Get Probe and Target Names from feature data
BrM_Targets <- BrM_IPA_featureData %>%
    as.data.frame() %>%
    select(ProbeName, TargetName) # only select ProbeName and TargetName

# AssayData
BrM_IPA_assayData <- read_csv("Input_files/DSP/BrM/BrM_IPA_assayData.csv") %>%
    left_join(BrM_Targets, by = "ProbeName") %>% # combine assayData and Targets
    column_to_rownames("TargetName") %>% # rename rows with protein name
    select(-ProbeName) %>% # remove extra ProbeName column
    as.matrix() # convert to matrix

# ProtocolData
BrM_IPA_protocolData <- read_csv("Input_files/DSP/BrM/BrM_IPA_protocolData.csv") %>%
    column_to_rownames("SampleID") %>% # rename rows with SampleID
    as.matrix() # convert to matrix

# PhenoData
BrM_IPA_phenoData <- read_csv("Input_files/DSP/BrM/BrM_IPA_phenoData.csv") %>%
    column_to_rownames("SampleID") %>% # rename rows with SampleID
    as.matrix() # convert to matrix

# Create Protein object
BrM_ProteinData <- suppressWarnings(NanoStringGeoMxSet(
    assayData = BrM_IPA_assayData,
    phenoData = AnnotatedDataFrame(as.data.frame(BrM_IPA_phenoData)),
    protocolData = AnnotatedDataFrame(
        as.data.frame(BrM_IPA_protocolData)
    ),
    featureData = AnnotatedDataFrame(
        as.data.frame(BrM_IPA_featureData)
    ),
    featureType = "Target",
    analyte = "Protein",
    check = FALSE
))

#----------------------------------------------------------------------------------------------------
# Save objects                                                                                                   
#----------------------------------------------------------------------------------------------------

save(modules, file = "Structure/DSP/BrM/BrM_Modules.RData")
save(BrM_RNAData, file = "Structure/DSP/BrM/BrM_RNAData_raw.RData")
save(BrM_ProteinData, file = "Structure/DSP/BrM/BrM_ProteinData_raw.RData")
