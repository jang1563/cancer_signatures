# Cancer Signatures Analysis - R Scripts

This directory contains annotated R scripts for analyzing cancer-related gene expression signatures in spaceflight studies.

## Repository

GitHub: https://github.com/[username]/cancer_signatures

## Scripts Overview

| Script | Description |
|--------|-------------|
| `01_volcano_plot_cancer_genes.R` | Volcano plots for differential expression of oncogenes/TSGs |
| `02_gsea_nasa_twins.R` | GSEA visualization for NASA Twins Study data |
| `03_cancer_genes_heatmap.R` | Heatmaps of cancer gene expression across cell types |
| `04_gsea_hallmark_c4_c6.R` | GSEA analysis for Hallmark, C4, and C6 gene sets |

## Data Sources

- **Differential Expression**: scRNA-seq analysis comparing post-flight to pre-flight conditions
- **Gene Sets**: MSigDB (Hallmark, C4 computational cancer, C6 oncogenic signatures)
- **Gene Lists**:
  - Curated oncogene/TSG list
  - OncoKB cancer gene database

## Timepoints

- **R+1**: Immediately post-flight (return + 1 day)
- **R+45&R+82**: Long-term post-flight recovery

## Cell Types Analyzed

- PBMC (bulk)
- CD4 T cells
- CD8 T cells
- Other T cells
- B cells
- NK cells
- CD14 Monocytes
- CD16 Monocytes
- Dendritic cells (DC)
- Other

## Output Formats

All visualizations are saved in both:
- PDF (vector format for publication)
- TIFF (300 DPI, LZW compression)

## Dependencies

```r
# Core packages
library(ggplot2)
library(dplyr)
library(tidyverse)

# Bioconductor
library(DESeq2)
library(fgsea)
library(ComplexHeatmap)
library(msigdbr)

# Visualization
library(ggrepel)
library(pheatmap)
library(circlize)
library(forcats)
library(patchwork)
```

## Usage

Run scripts in order (01-04) to reproduce all analyses. Ensure input data paths are correctly configured for your local environment.

## Contact

[Your Name]
[Your Email]
