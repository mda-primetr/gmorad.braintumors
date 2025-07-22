# Check and install required packages
required_packages <- c("jsonlite", "renv", "dplyr","purrr","magrittr", "tidyr")

for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
}

# Load necessary libraries
# Load each required package
for (pkg in required_packages) {
    library(pkg, character.only = TRUE)
}





# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

# Restore other packages from renv.lock
if (!file.exists("renv.lock")) {
    stop("renv.lock file not found. Please run renv::snapshot() first.")
}
cat("Restoring packages from renv.lock...\n")
renv::restore()

cat("Setup complete!\n")


# Read the renv.lock.orig file since this file contains the original package information
# renv.lock may have been modified by renv::snapshot() function
lock_data <- fromJSON("renv.lock.orig")


# Extract package information from the lock data

df_renv_packages <- names(lock_data$Packages) %>%
    map(function(pkg_name) {
        pkg <- lock_data$Packages[[pkg_name]]
        data.frame(
            package = pkg_name,
            version = pkg$Version,
            repo = ifelse(
                is.null(pkg$Repository),
                NA_character_,
                pkg$Repository
            ),
            stringsAsFactors = FALSE
        )
    }) %>%
    list_rbind() %>%
    separate(
        repo,
        into = c("repo_type", "repo_version"),
        sep = " ",
        fill = "right"
    )   %T>% 
    # Check package name, version and repo for manual inspection
    write.csv("renv_packages.csv", row.names = FALSE)

# Create a data frame of installed packages
# This will be used to check which packages are already installed
installed <- installed.packages()

installed_df <- data.frame(
    package = installed[, "Package"],
    version = installed[, "Version"],
    stringsAsFactors = FALSE
)



# List of required packages with specific versions  -----
list_of_required_packages <- df_renv_packages %>% 
    filter(
        package %in% c(
            "ANCOMBC", "tidyverse", "broom", "cardx", "ggalluvial",
            "ggchangepoint", "ggh4x", "ggrepel", "gtsummary",
            "janitor", "Maaslin2", "patchwork", "phyloseq",
            "tidytext", "vegan","microViz", "speedyseq","SCRuB"
        )
    )  %>% 
    pull(package)

# List of packages that may not be easily installed from CRAN or Bioconductor
list_custom_sources <- list(

    # MicroViz
    # If there is an issue with the renv snapshot restore, one can install the package using this R 
    # install.packages("src/local_packages/microViz_0.12.6.tar.gz", repos = NULL, type = "source")
    "microViz",
    
    

    # SCRuB 
    # If there is an issue with the renv snapshot restore, one can install the package using this R
    # install.packages("src/local_packages/SCRuB_0.0.1.tar.gz", repos = NULL, type = "source")

    "SCRuB",
  
    # SpeedySeq release version is behind the latest commit
     # If there is an issue with the renv snapshot restore, one can install the package directly
    # install.packages("src/local_packages/speedyseq_0.5.3.9021.tar.gz", repos = NULL, type = "source")
    "speedyseq"

)


# List of packages that need to be installed
(packages_to_install <- df_renv_packages %>% 
    dplyr::anti_join(installed_df, by = c("package", "version")) %>% 
    mutate(
        source_type = ifelse(
            package %in% list_of_required_packages,
            "required",
            "standard"
        )
    )) 
    

# If there is no mismatch between the installed packages and the required packages
if (nrow(packages_to_install) == 0) {
    cat("All required and standard packages are already installed.\n")
}