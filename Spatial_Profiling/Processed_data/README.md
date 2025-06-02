# Processed_data

This folder contains intermediate `.RData`files generated during the analysis of brain metastasis (BrM) and glioma samples using the DSP platform. These files are used to speed up downstream analysis and figure generation by avoiding the need to reprocess raw inputs.

## Folder Structure
```
├── Processed_data/       # Intermediate .RData files for DSP analysis
│   ├── BrM/              # Intermediate .RData files for BrM analysis
│   └── Glioma/           # Intermediate .RData files for glioma analysis
```

## How to Use

These `.RData` files are loaded in the following:
- **Figure generation scripts** in `Code/DSP/BrM/`

To regenerate, re-run the preprocessing scripts located in `Code/DSP/BrM/`.
  
## Note

Do not edit these files manually. They are intermediate products automatically created during the analysis pipeline and are reused across multiple steps.
