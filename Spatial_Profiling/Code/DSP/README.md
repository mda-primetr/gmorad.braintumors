# GeoMx DSP Analysis – Code Overview

This folder contains R scripts used to analyze GeoMx Digital Spatial Profiler (DSP) data across various tumor types. The pipeline supports RNA and protein data analysis from raw input through normalization, 16S classification, PCA, and differential expression.

## Script Overview

Scripts are numbered to indicate the recommended order of execution.

### RNA Pipeline

**1_Create_DSP_Object.R**  
Loads raw ROI-level RNA data and creates a `NanoStringGeoMxSet` object.  
*Output:* `.RData` file for downstream use.

**2_QC_RNA.R**  
Performs RNA probe- and ROI-level quality control. Filters out low-quality regions.  
*Output:* Filtered `.RData` object.

**3_Normalization_RNA.R**  
Applies normalization methods (e.g., Q3, background adjustment).  
*Output:* Normalized `.RData` object.

**4_Create_16S_Scores_RNA.R**  
Classifies ROIs as 16S-high or 16S-low using 25th–75th percentile thresholds.  
*Output:* Annotated `.RData` object with 16S classification.

**5_DGE_RNA.R**  
Runs differential gene expression analysis between 16S-high and 16S-low groups using linear mixed models.  
*Outputs:*  
- `.csv` results in `Output_files/DSP/Tables/`  
- Updated `.RData` object  

**6_Unsupervised_Clustering_RNA.R**  
Performs PCA to explore variance in spatial regions driven by 16S.  
*Output:* `.RData` object

### Protein Pipeline

**7_QC_Protein.R**  
Quality control of protein panel data, filtering out low-expression or non-informative markers.  
*Output:* QC-passed `.RData` object.

**8_Normalization_Protein.R**  
Normalizes protein expression data for downstream analysis.  
*Output:* Normalized `.RData` object.

**9_DPE_Protein.R**  
Performs differential protein expression analysis by 16S status using linear mixed models.  
*Outputs:*  
- `.csv` results in `Output_files/DSP/Tables/`  
- Updated `.RData` object  

**10_Pathway_activation.R**  
Computes module scores from normalized protein data.  
Generates  dotplots.  
*Output:* `.pdf` plots in `Output_files/DSP/Figures/`

## Output Summary

**Intermediate Files**  
`.RData` files saved in `Processed_data/<TumorType>/`

**Tables**  
Statistical outputs as `.csv` files in `Output_files/DSP/Tables/`

**Figures**  
Plots saved as `.pdf` in `Output_files/DSP/Figures/`

## Notes

- Tumor-type-specific folders (e.g., `BrM/`, `Glioma/`) contain copies of this workflow with matched metadata.
