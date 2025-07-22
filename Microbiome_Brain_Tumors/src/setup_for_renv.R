# setup.R
cat("Setting up R environment...\n")


# Define packages to install with versions
packages <- c(
    "ANCOMBC@2.8.1",
    "tidyverse@2.0.0",
    "broom@1.0.8",
    "cardx@0.2.5",
    "ggalluvial@0.12.5",
    "ggchangepoint@0.1.0",
    "ggh4x@0.3.0",
    "ggrepel@0.9.6",
    "gtsummary@2.1.0",
    "janitor@2.2.1",
    "Maaslin2@1.20.0",
    "patchwork@1.3.0", 
    "phyloseq@1.50.0",
    "tidytext@0.4.2",
    "vegan@2.7-1",
    "david-barnett/microViz@1dadd117bd3c5bd724c0be66e9572b4aa6443191",  #MicroViz@0.12.6
    "mikemc/speedyseq@0057652ff7a4244ccef2b786dca58d901ec2fc62", # speedyseq@0.5.3.9021
     "Shenhav-and-Korem-labs/SCRuB@fcbb8524190f0b27b7ad52cde232c8c4f59810e0"  #"SCRuB@0.0.1"
)

# Install each package with error handling
for (pkg in packages) {
    cat(paste0("Installing ", pkg, "...\n"))
    tryCatch(
        {
            renv::install(pkg, prompt = TRUE)
            cat(paste0("Successfully installed ", pkg, "\n"))
        },
        error = function(e) {
            cat(paste0("Error installing ", pkg, ": ", conditionMessage(e), "\n"))
        },
        warning = function(w) {
            cat(paste0("Warning while installing ", pkg, ": ", conditionMessage(w), "\n"))
        }
    )
}

renv::snapshot()

