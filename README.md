# Integrated Analysis of Differential Gene Expression, Enrichment, and Deconvolution in Luminal A and Basal-Like Breast Cancer

Cassi Chen and Fang (Amy) Liu

This project reproduces and extends a comparative analysis of Luminal A and Basal-like breast cancer subtypes using TCGA-BRCA RNA-seq data. We perform differential expression and enrichment analyses to identify subtype-specific biological patterns recovering key themes such as hormone-related signaling in Luminal A and cell cycle and extracellular matrix processes in Basal-like tumors. In addition, we extend the analysis with deconvolution to estimate tumor microenvironment composition, providing a complementary, cell-type–level perspective. Together, this workflow integrates gene-level, pathway-level, and cellular insights to better characterize differences between breast cancer subtypes.

# Data Loading and Preprocessing Files

## Scripts

- **`tcga_prepare.R`** — Main script to obtain TCGA-BRCA bulk samples from the Genomic Data Commons (GDC).

- **`tcga_v2.R`** — Main script to preprocess TCGA-BRCA bulk samples to have only Normal, Luminal A, and Basal samples.


# Differential Expression Analysis Files

## Scripts

- **`DE_analysis.R`** — Main analysis script that performs differential expression analysis on the TCGA-BRCA bulk samples.

## Outputs

- **`Volcano_LuminalA_vs_Normal.png`** — Volcano plots for Luminal A vs Normal, showing upregulated and downregulated genes.

- **`Volcano_Basal_vs_Normal.png`** — Volcano plots for Basal vs Normal, showing upregulated and downregulated genes.


# Enrichment Analysis Files

## Scripts

- **`enrichment_analysis.R`** — Main analysis script that runs enrichment analysis on the differentially expressed genes.

## Outputs

- **`Enrichment_LuminalA_vs_Normal.png`** — Enrichment Plots (Biological Process, Cell Component, Molecular Function, KEGG Pathway) for Luminal A vs Normal.

- **`Enrichment_Basal_vs_Normal.png`** — Enrichment Plots (Biological Process, Cell Component, Molecular Function, KEGG Pathway) for Basal vs Normal.


# MuSiC Deconvolution Files

## Scripts

- **`music_deconvolution.R`** — Main analysis script that runs MuSiC deconvolution on the TCGA-BRCA bulk samples using the Wu 2021 single-cell reference and validates the output against TCGA tumor-purity scores.

- **`music_figures.R`** — Plotting script that loads the saved deconvolution output and produces the purity-validation scatter plots and the per-subtype stacked composition figure.

## Outputs

- **`music_proportions.csv`** — Tidy table of MuSiC weighted cell-type proportions per tumor sample, joined with PAM50 subtype metadata for downstream analysis.

- **`music_stacked_composition.png`** — Stacked bar chart of mean MuSiC-estimated cell-type composition for Luminal A and Basal-like tumors.

- **`music_purity_validation_LumA.png`** — Scatter plots of MuSiC compartment proportions versus TCGA ABSOLUTE and ESTIMATE purity scores for Luminal A tumors, with per-panel Spearman correlations.

- **`music_purity_validation_Basal.png`** — Same purity-validation scatter plots as above but restricted to Basal-like tumors.

## Reference data

- **NOTE:** Single-cell RNA-seq reference atlas can be downloaded from Wu et al. 2021 (GSE176078); will contain 10X-formatted barcodes, genes, count matrix, and cell-level metadata
