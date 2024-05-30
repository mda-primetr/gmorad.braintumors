#----------------------------------------------------------------------------------------------------
# Code based on GeoMXTools vignette by Jason Reeves, Prajan Divakar, Nicole Ortogero,
# Maddy Griswold, Zhi Yang, Stephanie Zimmerman, Rona Vitancol and David Henderson (nanoString)
#----------------------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------------------
# Source files
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")

# Set working directory
datadir <- file.path("Input_files/DSP/Glioma")

#----------------------------------------------------------------------------------------------------
# RNA NanoStringGeoMxSet object
#----------------------------------------------------------------------------------------------------

# Files location
DCCFiles <- dir(file.path(datadir, "DCC"), pattern = ".dcc$", full.names = TRUE, recursive = TRUE)
PKCFiles <- dir(file.path(datadir, "PKC"), pattern = ".pkc$", full.names = TRUE, recursive = TRUE)
AnnotationFile <- file.path(datadir, "Glioma_Annotations.xlsx")

# Create RNA object
Glioma_RNAData <- suppressWarnings(readNanoStringGeoMxSet(
    dccFiles = DCCFiles,
    pkcFiles = PKCFiles,
    phenoDataFile = AnnotationFile,
    phenoDataSheet = "Annotation-RNA",
    phenoDataDccColName = "Sample_ID",
    protocolDataColNames = c("aoi", "roi"),
    configFile = NULL,
    analyte = "RNA",
    phenoDataColPrefix = "",
    experimentDataColNames = NULL
))

# Modules used
pkcs <- annotation(Glioma_RNAData)
modules <- gsub(".pkc", "", pkcs)

#----------------------------------------------------------------------------------------------------
# Protein NanoStringGeoMxSet object
#----------------------------------------------------------------------------------------------------

# Load featureData
Glioma_IPA_featureData <- read_csv("Input_files/DSP/Glioma/Glioma_IPA_featureData.csv") %>%
    arrange(ProbeName) %>% # sort by ProbeName
    column_to_rownames("Target_Name") %>% # rename rows with protein name
    as.matrix() # convert to matrix

# Get Probe and Target Names from feature data
Glioma_Targets <- Glioma_IPA_featureData %>%
    as.data.frame() %>%
    select(ProbeName, TargetName) # only select ProbeName and TargetName

# AssaydData
Glioma_IPA_assayData <- read_csv("Input_files/DSP/Glioma/Glioma_IPA_assayData.csv") %>%
    left_join(Glioma_Targets, by = "ProbeName") %>% # combine assayData and Targets
    column_to_rownames("TargetName") %>% # rename rows with protein name
    select(-ProbeName) %>% # remove extra ProbeName column
    as.matrix() # convert to matrix

# ProtocolData
Glioma_IPA_protocolData <- read_csv("Input_files/DSP/Glioma/Glioma_IPA_protocolData.csv") %>%
    column_to_rownames("SampleID") %>% # rename rows with SampleID
    as.matrix() # convert to matrix

# PhenoData
Glioma_IPA_phenoData <- read_csv("Input_files/DSP/Glioma/Glioma_IPA_phenoData.csv") %>%
    column_to_rownames("SampleID") %>% # rename rows with SampleID
    as.matrix() # convert to matrix

# Create Protein object
Glioma_ProteinData <- suppressWarnings(NanoStringGeoMxSet(
    assayData = Glioma_IPA_assayData,
    phenoData = AnnotatedDataFrame(as.data.frame(Glioma_IPA_phenoData)),
    protocolData = AnnotatedDataFrame(
        as.data.frame(Glioma_IPA_protocolData)
    ),
    featureData = AnnotatedDataFrame(
        as.data.frame(Glioma_IPA_featureData)
    ),
    featureType = "Target",
    analyte = "Protein",
    check = FALSE
))

#---------------------------------------------------------------------------------------------------
# Save objects
#---------------------------------------------------------------------------------------------------

save(modules, file = "Structure/DSP/Glioma/Glioma_Modules.RData")
save(Glioma_RNAData, file = "Structure/DSP/Glioma/Glioma_RNAData_raw.RData")
save(Glioma_ProteinData, file = "Structure/DSP/Glioma/Glioma_ProteinData_raw.RData")
