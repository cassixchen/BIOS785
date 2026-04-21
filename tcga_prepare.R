library(TCGAbiolinks)
library(SummarizedExperiment)

# windows only: increase memory limit to 28GB
memory.limit(size = 28000)

# ============================================================
# STEP 1: Query GDC
# ============================================================
query_brca <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal"),
  access = "open"
)

# Verify sample counts (expect ~1111 tumor, 113 normal)
output_query_brca <- getResults(query_brca)
cat("Sample counts from query:\n")
print(table(output_query_brca$sample_type))

# ============================================================
# STEP 2: Download files
# ============================================================
GDCdownload(query_brca, method = "api", files.per.chunk = 20)

# ============================================================
# STEP 3: Prepare and save full dataset
# ============================================================
brca_data <- GDCprepare(
  query_brca,
  save = TRUE,
  save.filename = "brca_data_raw.rda"
)
cat("brca_data_raw.rda saved successfully.\n")

# ============================================================
# STEP 4: Get PAM50 subtypes and define sample groups
# ============================================================
subtypes <- TCGAquery_subtype(tumor = "BRCA")
cat("PAM50 subtype counts:\n")
print(table(subtypes$BRCA_Subtype_PAM50))

# Normal samples
normal_samples <- colData(brca_data)$barcode[
  colData(brca_data)$sample_type == "Solid Tissue Normal"
]

# LumA and Basal patient IDs (with NA guard)
lumA_patients <- subtypes$patient[
  !is.na(subtypes$BRCA_Subtype_PAM50) & subtypes$BRCA_Subtype_PAM50 == "LumA"
]
basal_patients <- subtypes$patient[
  !is.na(subtypes$BRCA_Subtype_PAM50) & subtypes$BRCA_Subtype_PAM50 == "Basal"
]

# Translate patient IDs to sample barcodes
lumA_samples <- colData(brca_data)$barcode[
  colData(brca_data)$patient %in% lumA_patients &
  colData(brca_data)$sample_type == "Primary Tumor"
]
basal_samples <- colData(brca_data)$barcode[
  colData(brca_data)$patient %in% basal_patients &
  colData(brca_data)$sample_type == "Primary Tumor"
]

# Combine into final barcode list
final_barcodes <- c(normal_samples, lumA_samples, basal_samples)

# Sanity checks
stopifnot(length(normal_samples) > 0)
stopifnot(length(lumA_samples) > 0)
stopifnot(length(basal_samples) > 0)
cat(sprintf("Sample counts -> Normal: %d | LumA: %d | Basal: %d\n",
            length(normal_samples), length(lumA_samples), length(basal_samples)))

# ============================================================
# STEP 5: Save final_barcodes and free memory
# ============================================================
save(final_barcodes, normal_samples, lumA_samples, basal_samples, file = "final_barcodes.rda")
cat("final_barcodes.rda saved successfully.\n")

# Free the large object from RAM — no longer needed in this script
rm(brca_data)
gc()
cat("Done. You can now run tcga_v2.R using load() for both .rda files.\n")
