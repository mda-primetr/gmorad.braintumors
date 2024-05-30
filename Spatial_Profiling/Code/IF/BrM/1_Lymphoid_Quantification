#----------------------------------------------------------------------------------------------------
# Source files
#----------------------------------------------------------------------------------------------------

source("src/Libraries.R")
source("src/Functions.R")
load("Structure/DSP/BrM/BrM_target_Data_WTA_PCA.RData")

#----------------------------------------------------------------------------------------------------
# Load and transform IF data
#----------------------------------------------------------------------------------------------------

# Load Lymphoid IF data
BrM_Lymphoid <- read_csv("Input_files/IF/BrM_Lymphoid_IF.csv") %>%
    na.omit()

# Load BrM phenoData from DSP
df_pData <- pData(BrM_target_Data_WTA) %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    dplyr::rename("Sample_ID" = "rowname")

# Join IF data with phenoData
df_data_Lymphoid <- df_pData %>%
    inner_join(BrM_Lymphoid, by = "ROI") %>%
    filter(test_16SrRNA_status_25_75 != "NA") %>%
    as.data.frame() %>%
    column_to_rownames("Sample_ID")

# Select columns to be evaluated (log)
m_data_Lymphoid_log <- df_data_Lymphoid %>%
    select(27:43) %>%
    mutate(across(.cols = c(1:17), .fns = function(x) x + 1)) %>%
    mutate(across(.cols = c(1:17), .fns = function(x) log(x, 2))) %>%
    as.matrix() %>%
    t()

# Select columns to be evaluated (non-log)
m_data_Lymphoid <- df_data_Lymphoid %>%
    select(27:43) %>%
    as.matrix() %>%
    t()

#----------------------------------------------------------------------------------------------------
# Create NanoStringGeoMxSet object to use Linear Mix Model
#----------------------------------------------------------------------------------------------------

# Subset BrM DSP object to use as basis for new object containing IF data
BrM_target_Data_WTA_16S_25_75 <- BrM_target_Data_WTA[
    , which(pData(BrM_target_Data_WTA)$test_16SrRNA_status_25_75 != "NA" &
        pData(BrM_target_Data_WTA)$RNA_ID %in% df_data_Lymphoid$RNA_ID)
]

# Create NanoStringGeoMxSet with log transformed counts
Lymphoid_IF_log <- NanoStringGeoMxSet(
    assayData = m_data_Lymphoid_log,
    phenoData = AnnotatedDataFrame(pData(BrM_target_Data_WTA_16S_25_75)),
    protocolData = protocolData(BrM_target_Data_WTA_16S_25_75),
    featureData = AnnotatedDataFrame(as.data.frame(m_data_Lymphoid_log)),
    featureType = "IF",
    check = FALSE
)

# Create NanoStringGeoMxSet with untransformed counts
Lymphoid_IF <- NanoStringGeoMxSet(
    assayData = m_data_Lymphoid,
    phenoData = AnnotatedDataFrame(pData(BrM_target_Data_WTA_16S_25_75)),
    protocolData = protocolData(BrM_target_Data_WTA_16S_25_75),
    featureData = AnnotatedDataFrame(as.data.frame(m_data_Lymphoid)),
    featureType = "IF",
    check = FALSE
)

#---------------------------------------------------------------------------------------------------
# Differential expression (Linear Mix Model)
#---------------------------------------------------------------------------------------------------

# 25-75% log
results <- c()

for (compartment in "BrM") {
    ind <- pData(Lymphoid_IF_log)$Type == compartment
    mixedOutmc <-
        mixedModelDE(Lymphoid_IF_log[, ind],
            elt = "exprs",
            modelFormula = ~ test_16SrRNA_status_25_75 + Primary_Tumor +
                (1 + test_16SrRNA_status_25_75 | Patient), # Same as DSP
            groupVar = "test_16SrRNA_status_25_75",
            nCores = parallel::detectCores(),
            multiCore = FALSE
        )
    r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
    tests <- rownames(r_test)
    r_test <- as.data.frame(r_test)
    r_test$Contrast <- tests
    r_test$Cell <-
        unlist(lapply(
            colnames(mixedOutmc),
            rep, nrow(mixedOutmc["lsmeans", ][[1]])
        ))
    r_test$Subset <- compartment
    r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "BH")
    r_test <- r_test[, c(
        "Cell", "Subset", "Contrast", "Estimate",
        "Pr(>|t|)", "FDR"
    )]
    results <- rbind(results, r_test)
}
BrM_Lymphoid_IF_by_16S_25_75_log <- as.data.frame(results) %>%
    dplyr::rename(log2FC = Estimate) # rename for better readability

# Save object
save(BrM_Lymphoid_IF_by_16S_25_75_log,
    file =
        "Structure/IF/BrM_Lymphoid_IF_by_16S_25_75_log.RData"
)

# Write csv
write.csv(
    BrM_Lymphoid_IF_by_16S_25_75_log,
    "Output_files/IF/BrM_Lymphoid_IF_by_16S_25_75_log.csv"
)

# 25-75% non-log
results <- c()

for (compartment in "BrM") {
    ind <- pData(Lymphoid_IF)$Type == compartment
    mixedOutmc <-
        mixedModelDE(Lymphoid_IF[, ind],
            elt = "exprs",
            modelFormula = ~ test_16SrRNA_status_25_75 + Primary_Tumor +
                (1 + test_16SrRNA_status_25_75 | Patient),
            groupVar = "test_16SrRNA_status_25_75",
            nCores = parallel::detectCores(),
            multiCore = FALSE
        )
    r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
    tests <- rownames(r_test)
    r_test <- as.data.frame(r_test)
    r_test$Contrast <- tests
    r_test$Cell <-
        unlist(lapply(
            colnames(mixedOutmc),
            rep, nrow(mixedOutmc["lsmeans", ][[1]])
        ))
    r_test$Subset <- compartment
    r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "BH")
    r_test <- r_test[, c(
        "Cell", "Subset", "Contrast", "Estimate",
        "Pr(>|t|)", "FDR"
    )]
    results <- rbind(results, r_test)
}
BrM_Lymphoid_IF_by_16S_25_75_non_log <- as.data.frame(results)

# Save object
save(BrM_Lymphoid_IF_by_16S_25_75_non_log,
    file =
        "Structure/IF/BrM_Lymphoid_IF_by_16S_25_75_non_log.RData"
)

# Write csv
write.csv(
    BrM_Lymphoid_IF_by_16S_25_75_non_log,
    "Output_files/IF/BrM_Lymphoid_IF_by_16S_25_75_non_log.csv"
)
#----------------------------------------------------------------------------------------------------
# Plot immune cells
#----------------------------------------------------------------------------------------------------

# Non-log
# Change order of 16S category
df_data_Lymphoid_plot <- as.data.frame(t(assayData(Lymphoid_IF)$exprs)) %>%
    mutate(test_16SrRNA_status_25_75 = factor(pData(Lymphoid_IF)$test_16SrRNA_status_25_75,
        levels = c("Low", "High")
    ))

# Plot CD16+ CD56- GZMB- Cells
ggplot(df_data_Lymphoid_plot, aes(
    x = test_16SrRNA_status_25_75, y = `CD16+ CD56- GZMB- Cells`,
    fill = test_16SrRNA_status_25_75
)) +
    geom_boxplot(alpha = 0.5, show.legend = FALSE, width = 0.3) +
    geom_dotplot(
        dotsize = 0.8, binaxis = "y", stackdir = "center",
        position = position_dodge(0.75), show.legend = FALSE
    ) +
    scale_fill_manual(values = c("darkblue", "darkred")) +
    theme_bw(base_size = 12) +
    theme(
        legend.position = "none", panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
    ) +
    labs(x = "", y = "CD16+ CD56- GZMB-")

# Save pdf
ggsave("Figures/Figure3F.pdf", width = 1.85, height = 2.2)
