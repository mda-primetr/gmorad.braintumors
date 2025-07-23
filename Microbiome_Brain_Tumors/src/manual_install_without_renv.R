# This script installs specific versions of R packages without using renv.
# Please ensure you have R version 4.4.3 with the necessary dependencies installed.
# Also note that we have to use Bioconductor version 3.20 for compatibility with the packages.


# Install remotes if needed
if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes")
}

# Deactivate any existing renv environment

# Check if renv is available and active before deactivating
if (requireNamespace("renv", quietly = TRUE)) {
    if (renv::status()$active) {
        message("Deactivating renv environment...")
        renv::deactivate()
    } else {
        message("No active renv environment found.")
    }
} else {
    message("renv package is not installed.")
}

# Install packages with specific versions using remotes
remotes::install_version("tidyverse", version = "2.0.0")
remotes::install_version("broom", version = "1.0.8")
remotes::install_version("cardx", version = "0.2.5")
remotes::install_version("ggalluvial", version = "0.12.5")
remotes::install_version("ggchangepoint", version = "0.1.0")
remotes::install_version("ggh4x", version = "0.3.0")
remotes::install_version("ggrepel", version = "0.9.6")
remotes::install_version("gtsummary", version = "2.1.0")
remotes::install_version("janitor", version = "2.2.1")
remotes::install_version("patchwork", version = "1.3.0")
remotes::install_version("tidytext", version = "0.4.2")
remotes::install_version("vegan", version = "2.7-1")

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Bioconductor packages with specific versions
BiocManager::install(version = "3.20")

BiocManager::install("ANCOMBC", version = "2.8.1")
BiocManager::install("Maaslin2", version = "1.20.0")
BiocManager::install("phyloseq", version = "1.50.0")


# Install packages from GitHub with specific commits
remotes::install_github("david-barnett/microViz@1dadd117bd3c5bd724c0be66e9572b4aa6443191")
remotes::install_github("mikemc/speedyseq@0057652ff7a4244ccef2b786dca58d901ec2fc62")
remotes::install_github("Shenhav-and-Korem-labs/SCRuB@fcbb8524190f0b27b7ad52cde232c8c4f59810e0")