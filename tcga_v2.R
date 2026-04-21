setwd("/Users/amyfliu/Documents/bio785_final")

library(TCGAbiolinks)
library(SummarizedExperiment)

# ============================================================
# STEP 1: Load pre-prepared data (1224 samples)
# Run tcga_prepare.R first if these files don't exist yet
# ============================================================
load("brca_data_raw.rda")
load("final_barcodes.rda")  # loads: final_barcodes, normal_samples, lumA_samples, basal_samples
brca_data <- data   # assign to new name
rm(data) 

# ============================================================
# STEP 2: Filter to the three groups of interest (881 samples)
# ============================================================
data_filtered <- brca_data[, colData(brca_data)$barcode %in% final_barcodes]

# Free full dataset from RAM immediately after filtering
rm(brca_data)
gc()

# ============================================================
# STEP 3: Create the 'condition' metadata column
# ============================================================
colData(data_filtered)$condition <- "Normal"
colData(data_filtered)$condition[colData(data_filtered)$barcode %in% lumA_samples] <- "LuminalA"
colData(data_filtered)$condition[colData(data_filtered)$barcode %in% basal_samples] <- "Basal"

# Convert to factor with Normal as reference level for DE analysis
colData(data_filtered)$condition <- factor(colData(data_filtered)$condition,
                                           levels = c("Normal", "LuminalA", "Basal"))

# ============================================================
# STEP 4: Verify final sample counts
# ============================================================
cat("Final sample counts:\n")
print(table(colData(data_filtered)$condition))

# Normal LuminalA    Basal 
# 113      571      197 

length(unique(colData(data_filtered)$patient))
#we have 881 samples from 786 unique patients

save(data_filtered, file = "data_filtered.rda")


# ============================================================
# STEP 5: Deduplicate to one sample per patient
# ============================================================

# Find patients with multiple samples
# patient_counts <- table(colData(data_filtered)$patient)
# dup_patients <- names(patient_counts[patient_counts > 1])
# 
# cat(sprintf("Number of patients with multiple samples: %d\n", length(dup_patients)))
# 
# # Pick one patient to inspect
# example_patient <- dup_patients[1]
# example_samples <- data_filtered[, colData(data_filtered)$patient == example_patient]
# 
# # Show their sample metadata
# colData(example_samples)[, c("barcode", "patient", "sample_type", "condition")]

# keeps the first occurrence of each patient and drops the rest
# keep <- !duplicated(colData(data_filtered)$patient)
# data_filtered <- data_filtered[, keep]

# Normal LuminalA    Basal 
# 77      530      179 

# Verify deduplication
# cat("Sample counts after deduplication:\n")
# print(table(colData(data_filtered)$condition))
# cat(sprintf("Unique patients: %d\n", length(unique(colData(data_filtered)$patient))))
