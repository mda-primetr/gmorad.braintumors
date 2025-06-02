# Spatial Analysis Repository
This repository contains code, data, and outputs for analyzing brain metastasis (BrM) and glioma samples across multiple spatial profiling platforms.


## Project Structure

```
Spatial_Profiling/
│
├── Code/                 # Final scripts used to generate figures and results
│   ├── DSP/              # Digital Spatial Profiling
│   │   ├── BrM/          # Code for brain metastasis samples
│   │   └── Glioma/       # Code for glioma samples
│   ├── SMI/              # Spatial Molecular Imaging 
│   │   ├── BrM/          # Code for brain metastasis samples
│   │   └── Glioma/       # Code for glioma samples
│   └── IF/               # Code for multiplex immunofluorescence imaging
│
│
├── Input_files/          # Input files for DSP, SMI, and IF 
│   ├── DSP/              
│   │   ├── BrM/          # Raw GeoMx files and metadata for BrM
│   │   └── Glioma/       # Raw GeoMx files and metadata for glioma
│   ├── SMI/              
│   │   ├── BrM/          # Slimmed Seurat objects or CSVs for BrM figure generation
│   │   └── Glioma/       # Slimmed Seurat objects or CSVs for Glioma figure generation
│   └── IF/               # Count files for immunofluorescence imaging 
│
├── Output_files/         # Final output results by platform
│   ├── DSP/
│   │   ├── Tables/       # CSVs: differential expression (DGE/DPE)
│   │   └── Figures/      # Final plots from DSP analysis (PDF)
│   ├── SMI/
│   │   └── Figures/      # Final plots from SMI analysis (PDF)
│   └── IF/
│       ├── Tables/       # CSVs: quantification stats
│       └── Figures/      # Final plots from IF analysis (PDF)
│
├── Processed_data/       # Intermediate .RData files for DSP analysis
│
└── src/                  # Shared libraries, custom functions, and plotting tools
```
## Platforms & Sample Types

- **DSP (GeoMx)** – spatial transcriptomics & proteomics profiling
- **SMI (CosMx)** – spatial transcriptomics using a 6K RNA panel and pan-16S probes
- **IF** – Multiplex immunofluorescence

Samples are categorized by:
- **BrM** – Brain metastasis
- **Glioma** – Primary gliomas
