# src

This folder contains shared R scripts for package loading and reusable plotting functions used across the DSP and SMI analysis pipelines.

## Contents

1. Libraries.R

Loads all required CRAN and Bioconductor packages for downstream analysis, including:
	•	Core tidyverse tools (dplyr, ggplot2, purrr)
	•	Spatial and differential analysis libraries (Seurat, lme4, GeoMxTools, etc)
	•	Visualization utilities (ComplexHeatmap, viridis, etc)

Missing packages are automatically installed.

2. Functions.R

Reusable helper functions for:
- Volcano Plots
- Protein QC
- Spatial Visualization
	
These functions are sourced across both GeoMx DSP and CosMx SMI figure scripts to standardize output.
