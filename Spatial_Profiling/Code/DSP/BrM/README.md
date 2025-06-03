# DSP Analysis – Brain Metastasis (BrM)

This folder contains R scripts used to analyze GeoMx Digital Spatial Profiler (DSP) data from brain metastasis (BrM) samples. The pipeline covers RNA and protein analysis from raw input through normalization, 16S scoring, PCA, and differential expression.

## Script Overview

Scripts are numbered to indicate the recommended run order.

### RNA Pipeline

1. **`1_BrM_Create_DSP_Object.R`**  
   Loads raw DSP ROI-level data and creates a NanoString GeoMx object. Saves a `.RData` file for downstream use.

2. **`2_BrM_QC_RNA.R`**  
   Applies RNA probe- and ROI-level quality control. Filters out low-quality data. Saves `.RData` for downstream use.

3. **`3_BrM_Normalization_RNA.R`**  
   Performs RNA normalization (e.g., Q3, background-adjusted). Saves normalized `.RData` object.

4. **`4_BrM_Create_16S_Scores_RNA.R`**  
   Classifies ROIs as 16S-high or 16S-low based on the 25th–75th quantiles. Saves `.RData` object.

5. **`5_BrM_DGE_RNA.R`**  
   Performs differential gene expression analysis between 16S-high and 16S-low regions using linear mixed models.  
   - Output: `.csv` tables in `Output_files/DSP/Tables/`  
   - Saves `.RData` for downstream use.

6. **`6_BrM_Unsupervised_Clustering_RNA.R`**  
   Performs PCA to assess 16S-mediated variance across spatial regions. Saves `.RData` object.

### Protein Pipeline

7. **`7_BrM_QC_Protein.R`**  
   Conducts QC on protein panel data. Filters out non-informative or low-expressing markers. Saves `.RData`.

8. **`8_BrM_Normalization_Protein.R`**  
   Normalizes protein expression data. Saves `.RData` object for downstream use.

9. **`9_BrM_DPE_Protein.R`**  
   Performs differential protein expression analysis by 16S status.  
   - Output: `.csv` tables in `Output_files/DSP/Tables/`  
   - Saves `.RData` object.

10. **`10_Pathway_activation.R`**  
    Computes module scores from normalized protein expression. Generates dotplots for visualization. 

## Output Summary

- **Intermediate files:**  
  `.RData` objects saved to `Processed_data/BrM/`

- **Tables:**  
  `.csv` outputs saved to `Output_files/DSP/Tables/`

- **Figures:**  
  `.pdf` plots saved to `Output_files/DSP/Figures/`
