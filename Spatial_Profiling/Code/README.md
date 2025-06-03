# Spatial Analysis Code

This folder contains all analysis scripts used in the spatial profiling of brain metastasis (BrM) and glioma samples. Data were generated using the following platforms:

- **DSP** – GeoMx Digital Spatial Profiler (Bruker Spatial Biology)
- **SMI** – CosMx Spatial Molecular Imager (Bruker Spatial Biology)
- **IF** – Multiplex immunofluorescence imaging

Scripts are organized by platform and tumor type: **brain metastasis (BrM)** and **glioma**.

## Directory Structure

```
Code/
├── DSP/            # Full DSP analysis starting from raw input
│   ├── BrM/        # Scripts for brain metastasis samples
│   └── Glioma/     # Scripts for glioma samples
│
├── SMI/            # SMI analysis using slimmed .rds or .csv files
│   ├── BrM/        # SMI analysis for brain metastasis
│   └── Glioma/     # SMI analysis for glioma
│
└── IF/             # Quantification and analysis of multiplex IF data
```

## Platform Details

### DSP (GeoMx)
- Starts from **raw DSP output files** (DCC and PKC files).
- Includes QC, normalization, 16S region scoring, differential expression, and PCA.
- Intermediate `.RData` objects are saved to `Processed_data/` to reduce rerun time.

### SMI (CosMx)
- Scripts start from **slimmed Seurat objects** or **pre-processed `.csv` files** in `Input_files/SMI/`.
- Includes AUC module scoring, tumor and neighborhood annotations, heatmaps and spatial polygon visualization.
- Downstream figuresare written to `Output_files/SMI/`.

### IF (Multiplex Immunofluorescence)
- Scripts generate plots and evaluate statistically significant differences across regions.

## Output Organization

Each pipeline outputs results to:
- `Processed_data/` — intermediate `.RData` files
- `Output_files/[Platform]/Figures/` — publication-ready plots (PDF)
- `Output_files/[Platform]/Tables/` — CSV tables with DE results, scores, and annotations

## Tumor Types
Each platform’s scripts are organized into:
- **BrM/** – Brain metastasis
- **Glioma/** – Primary gliomas

## How to Use
1. Choose a platform folder: `DSP/`, `SMI/`, or `IF/`
2. Choose a tumor type: `BrM/` or `Glioma/`
3. Run scripts in numerical order

Scripts are modular but designed to work in sequence, with each step building on the previous one.

## R Package Dependencies

All required R packages are automatically installed when missing. Each script includes code to check for and install:
- CRAN packages
- Bioconductor packages

No manual installation is required. Packages will be installed from CRAN or Bioconductor the first time a script is run.

