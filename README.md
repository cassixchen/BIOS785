RAURRRRRRRR


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
