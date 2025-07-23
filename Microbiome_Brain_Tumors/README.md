# Code for Microbiome analysis 


## R packages:

Please use `renv` R package and packages listed below as R list. Make sure to use R version 4.4.3 from (https://cran.r-project.org/bin) and Bioconductor release version 3.20. This analysis was conducted using macOS 15.5 with  ARMv9 CPU architecture. 


```
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
```

## Initial setup:
This initial setup is **optional** if the required packages are installed correctly using `renv::restore()` function. However, if any issue arises, one can recreate the environment using this script `setup_for_renv.R` located within `src` folder. Within terminal type these commands. 

```
pwd # to make sure one is in Microbiome_Brain_Tumors folder

Rscript src/setup_for_renv.R
```

This script will make sure the above required packages are installed and and snapshot is created with these package. It will update the renv.lock file so we have included the renv.lock.orig file to compare for resolving any issues. 


## Helper script:
Please use `check_packages.R` script for restoring/installing or   checking packages installed in the renv.lock file. 

```
pwd # to make sure one is in Microbiome_Brain_Tumors folder
Rscript src/check_packages.R
```

If all the installed packages are same as the one in renv.lock.orig file then it should print


> Restoring packages from renv.lock...
> 
> The library is already synchronized with the lockfile.
> 
> Setup complete!
> 
> [1] package      version      repo_type    repo_version source_type 
> 
> <0 rows> (or 0-length row.names)
> 
> All required and standard packages are already installed.


## Manual Install:
Use this option if you are comfortable installing the packages globally and have the required R version. One can also use Rig (https://github.com/r-lib/rig) package to install and manage multiple versions of R. Install script is available within `src` folder. 

```
Rscript src/manual_install_without_renv.R
```


## Script structure

If the installation part looks ok. Proceed to the microbiome folder if not already in that folder by using these commands. 

```
cd Microbiome_Brain_Tumors
```

Run R scripts in order using this one liner. 

```
for dir in $(ls -d ./scripts/[0-9]*_* | sort -n); do echo "Directory: $dir"; find "$dir" -name "*.R" | sort | xargs -I{} sh -c 'echo "Running: {}"; Rscript "{}"'; done
```


The `src folder` contains all the scripts for common function and loading libraries.

The `script folder` contains all the scripts under numbered folder to guide the users for order of executing scripts.

`data`  folder contains

 * `raw_data` :  raw data files which are used to generate figures and files in `processed_data` folder. 
 * `metadata` : metadata related to the sequencing information.
 * `processed_data` : Data files generated while running this analysis script.

`output` folder contains

* `figures` : figures generated during the analysis.
* `tables` : Summary and stats table generated during the analysis.

``session_info.txt`` : This file contains all the packages that are loaded in the R environment during the analysis. 

`WGS` and `16S` folders are created where necessary to indicate which sequencing data was used to generate the output. 


## Script workflow 
![Schema of script order](script_workflow.png)

