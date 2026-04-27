# ============================================================
# MuSiC: Bulk RNA-seq deconvolution using a single-cell reference
# Extension analysis for TCGA-BRCA - tumor-only (LumA vs Basal).
# Normal samples excluded because the Wu 2021 reference was derived
# from primary tumors and lacks adipocyte populations that dominate
# normal breast tissue - Normal-sample deconvolution is unreliable.
#
# Reference: Wang et al. (2019) Nat Commun
# Tutorial:  https://xuranw.github.io/MuSiC/articles/MuSiC.html
# ============================================================

# ---- One-time install ----
# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::install(c("TOAST", "SingleCellExperiment", "SummarizedExperiment",
#                        "Biobase", "limma", "org.Hs.eg.db", "GEOquery", "DropletUtils"), force = T)
# if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
# remotes::install_github("xuranw/MuSiC")

library(MuSiC)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(org.Hs.eg.db)
library(DropletUtils)
library(GEOquery)
library(TCGAbiolinks)
library(dplyr)
library(tidyr)
library(ggplot2)

setwd("/Users/amyfliu/Documents/bio785_final")

# ============================================================
# STEP 1: Prepare the bulk input
# MuSiC wants a RAW counts matrix: genes (rows) x samples (cols).
# Gene identifiers MUST match the single-cell reference.
# ============================================================
load("data_filtered.rda")   # SummarizedExperiment from tcga_v2.R

# ---- Tumor-only filter ----
# Keep LumA + Basal primary tumors; drop the 113 Solid Tissue Normal samples.
# Reference lacks adipocytes, so Normal deconvolution misattributes ~50-70% of
# adipose signal to CAF/endothelial/PVL. Keep Normals out of the main analysis.
data_filtered <- data_filtered[, colData(data_filtered)$condition %in% c("LuminalA", "Basal")]
colData(data_filtered)$condition <- droplevels(colData(data_filtered)$condition)

# Raw unstranded counts
bulk_counts <- assay(data_filtered, "unstranded")

# The STAR - Counts workflow returns versioned ENSEMBL IDs (e.g. ENSG...14).
# Strip the version, then map to HGNC symbols 
rownames(bulk_counts) <- sub("\\..*$", "", rownames(bulk_counts))

symbols <- mapIds(org.Hs.eg.db,
                  keys     = rownames(bulk_counts),
                  column   = "SYMBOL",
                  keytype  = "ENSEMBL",
                  multiVals = "first") #keep only first match if multiple matches found

keep <- !is.na(symbols) & !duplicated(symbols)
bulk_counts <- bulk_counts[keep, ]
rownames(bulk_counts) <- symbols[keep]

# Keep sample-level metadata for plotting / stats later
bulk_meta <- as.data.frame(colData(data_filtered))[, c("barcode", "patient", "condition")]

cat(sprintf("Bulk matrix: %d genes x %d samples\n",
            nrow(bulk_counts), ncol(bulk_counts)))
# After tumor-only filter: ~36650 genes x 768 samples (571 LumA + 197 Basal)

# ============================================================
# STEP 2: Load the single-cell reference
#
# Reference: Wu et al. 2021, Nat Genet
#   "A single-cell and spatially resolved atlas of human breast cancers"
#   ~130k cells from 26 patients, covers tumor + TME
#   GSE176078  /  Broad SCP: https://singlecell.broadinstitute.org/single_cell/study/SCP1039
#
# Expected object: a SingleCellExperiment with
#   assay(sce, "counts")       raw counts
#   colData(sce)$celltype      cell-type label  (pass this to `clusters`)
#   colData(sce)$sampleID      donor / subject  (pass this to `samples`)
# ============================================================

sce <- read10xCounts("Wu_etal_2021_BRCA_scRNASeq", col.names = TRUE)
meta <- read.csv("Wu_etal_2021_BRCA_scRNASeq/metadata.csv", row.names = 1)
meta <- meta[colnames(sce), ]  # align metadata to sce
colData(sce)$celltype <- meta$celltype_major
colData(sce)$celltype_minor <- meta$celltype_minor
colData(sce)$sampleID <- meta$orig.ident
colData(sce)$subtype <- meta$subtype

# Quick sanity checks before running deconvolution:
# stopifnot("counts" %in% assayNames(sce))
# stopifnot(all(c("celltype", "sampleID") %in% colnames(colData(sce))))
# table(sce$celltype)
# length(unique(sce$sampleID))   # MuSiC needs multiple subjects to estimate cross-subject variance
# length(intersect(rownames(bulk_counts), rownames(sce))) #~18800 overlapping genes 

# ============================================================
# STEP 3: Run MuSiC deconvolution
# `select.ct` lets you restrict to meaningful, well-represented cell types.
# Low-abundance types (< ~20 cells per subject) tend to be noisy - drop them.
# ============================================================
celltypes_to_use <- c(
  "Cancer Epithelial", "Normal Epithelial", "CAFs",
  "Endothelial", "T-cells", "B-cells",
  "Myeloid", "Plasmablasts", "PVL"
)

# Mean library size per cell type from the sc reference
# Cancer epithelial cells and CAFs have 5–10× more mRNA per cell than resting T/B cells
# So a bulk sample that's 15% T-cells by count only contributes ~3% of total mRNA from T-cells
cell_size_df <- data.frame(
  celltype = names(table(sce$celltype)),
  size     = sapply(split(colSums(assay(sce, "counts")), sce$celltype), mean)
)

music_props <- music_prop(
  bulk.mtx  = bulk_counts,
  sc.sce    = sce,
  clusters  = "celltype",
  samples   = "sampleID",
  select.ct = celltypes_to_use,
  cell_size = cell_size_df,      # <- size correction
  verbose   = TRUE
)

# Used 18728 common genes...
# Used 9 cell types in deconvolution...

# Returned list:
#   $Est.prop.weighted   - MuSiC weighted estimates (use these)
#   $Est.prop.allgene    - NNLS baseline (analogous to CIBERSORT-lite)
#   $Weight.gene         - gene weights
#   $r.squared.full      - per-sample R^2 between observed bulk and reconstructed bulk

saveRDS(music_props, file = "music_deconvolution.rds")

# ============================================================
# PURITY CROSS-CHECK
# ============================================================
# Pull ABSOLUTE purity + ESTIMATE scores for the samples
purity <- TCGAbiolinks::Tumor.purity

props <- as.data.frame(music_props$Est.prop.weighted)
props$barcode <- rownames(props)
props <- merge(props, bulk_meta, by = "barcode")
write.csv(props, "music_proportions.csv", row.names = FALSE)
aggregate(. ~ condition, data = props[, c(celltypes_to_use, "condition")], FUN = mean)

#LumA should have high Cancer Epithelial, moderate CAFs, modest immune. 
#Basal should have elevated T-cells / Myeloid and lower Normal Epithelial

# Tumor.purity is pan-cancer and can have multiple rows per patient
# (primary + metastatic + duplicate vials). Filter to BRCA, match on a
# 15-char sample-level ID (TCGA-XX-XXXX-NN) instead of 12-char patient ID,
# and drop any remaining duplicates so the merge doesn't inflate row counts.
if ("Cancer.type" %in% colnames(purity)) {
  purity <- purity[purity$Cancer.type == "BRCA", ]
}
purity$sample_id <- substr(purity$Sample.ID, 1, 15)
purity <- purity[!duplicated(purity$sample_id), ]

props$sample_id <- substr(props$barcode, 1, 15)

#merge the tables so that we have purity estimates for each sample
df <- merge(props, purity, by = "sample_id")

# Correlate estimated cancer-epithelial fraction with ABSOLUTE tumor purity
# ESTIMATE comparison - bulk-derived, so less platform mismatch 
#     usually gives a stronger signal than ABSOLUTE
cor.test(df[["Cancer Epithelial"]], as.numeric(df$ABSOLUTE),   method = "spearman", exact = F) #0.305
cor.test(df[["Cancer Epithelial"]], as.numeric(df$ESTIMATE),   method = "spearman", exact = F) #0.712

# Stromal compartment should correlate NEGATIVELY with purity
df$stromal <- df[["CAFs"]] + df[["Endothelial"]] + df[["PVL"]]
cor.test(df$stromal, as.numeric(df$ABSOLUTE), method = "spearman", exact = F) #-0.269
cor.test(df$stromal, as.numeric(df$ESTIMATE), method = "spearman", exact = F) #-0.514

# ============================================================
# STEP 4: Tidy output + per-condition comparison
# ============================================================

props_long <- props %>%
  pivot_longer(all_of(celltypes_to_use),
               names_to = "celltype", values_to = "prop")

# Boxplots: cell-type proportions stratified by PAM50 subtype
p <- ggplot(props_long, aes(x = condition, y = prop, fill = condition)) +
  geom_boxplot(outlier.size = 0.4) +
  facet_wrap(~ celltype, scales = "free_y") +
  theme_bw(base_size = 11) +
  labs(x = NULL, y = "Estimated cell-type proportion",
       title = "MuSiC deconvolution: TCGA-BRCA tumors (LumA vs Basal)")

p

# ggsave("music_proportions_by_condition.png", p,
#        width = 10, height = 7, dpi = 300)
# 
# # Quick hypothesis tests: LumA vs Basal per cell type (Wilcoxon, BH-adjusted)
# tests <- props_long %>%
#   filter(condition %in% c("LuminalA", "Basal")) %>%
#   group_by(celltype) %>%
#   summarise(p_value = wilcox.test(prop ~ condition)$p.value, .groups = "drop") %>%
#   mutate(padj = p.adjust(p_value, method = "BH")) %>%
#   arrange(padj)

# write.csv(tests, "music_celltype_tests.csv", row.names = FALSE)
# cat("\nDone. Outputs:\n",
#     "  music_deconvolution.rds       - full MuSiC result\n",
#     "  music_proportions.csv         - per-sample estimated proportions + metadata\n",
#     "  music_proportions_by_condition.png\n",
#     "  music_celltype_tests.csv      - LumA vs Basal Wilcoxon + BH-adjusted p\n",
#     sep = "")
