#----------------------------------------------------------------------------------------------------
# Install and Load Required Packages
#----------------------------------------------------------------------------------------------------

# Define required packages
cran_packages <- c(
 "tidyverse", "ggforce", "ggpubr", "ggrepel", "scales", "ggrastr", "ggsci", "patchwork", "gridExtra", 
  "data.table", "RColorBrewer", "viridis", "kableExtra", "cowplot", "DT", "spdep", "ggsignif", "MASS", 
  "stringr", "forcats", "sf", "purrr", "tidyr", "rlang")

bioc_packages <- c(
  "ComplexHeatmap", "circlize", "NanoStringNCTools", "GeomxTools", "GeoMxWorkflows", "Seurat", "SeuratObject",
  "lmerTest", "lme4", "clusterProfiler", "org.Hs.eg.db", "enrichplot")

# Install CRAN packages if missing
cran_missing <- setdiff(cran_packages, rownames(installed.packages()))
if (length(cran_missing) > 0) {
  install.packages(cran_missing, dependencies = TRUE)
}

# Install Bioconductor packages if missing
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
bioc_missing <- setdiff(bioc_packages, rownames(installed.packages()))
if (length(bioc_missing) > 0) {
  BiocManager::install(bioc_missing, update = FALSE, ask = FALSE)
}

# Load all packages
all_packages <- c(cran_packages, bioc_packages)
invisible(lapply(all_packages, function(pkg) {
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}))

