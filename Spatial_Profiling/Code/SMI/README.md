# CosMx SMI Analysis – Code Overview

This folder contains R scripts used to generate key figures for the CosMx Spatial Molecular Imaging (SMI) analysis of glioma and brain metastasis (BrM) samples.

All input data are pre-processed and stored as intermediate .rds` files.

## Figure Types

- **UMAPs**  
  Visualization of cell types and normalized 16S signal.

- **Dotplots for Cell Typing**  
  Percentage of marker gene expression per cell type.

- **Differential Expression Dotplots**  
  16S-associated genes shown per patient and neighborhood.

- **AUCell Violin Plots**  
  Module scores custom programs across 16S+ tumor cells and their neighborhoods.

- **Spatial Plots**  
  Example tissue-level visualization of key gene programs.

## Input Summary

- Slimmed Seurat .rds` objects 
- Spatial coordinates and cell metadata (including 16S classification)
- Stored in `Input_files/SMI`

## Output Summary

- **Figures:**  
  `.pdf` plots saved to `Output_files/SMI/Figures/`
