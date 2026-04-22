# setwd("/Users/amyfliu/Documents/bio785_final")
setwd("/Users/cassixchen/Desktop/unc/BIOS785/Final")

library(edgeR)
library(limma)
library(SummarizedExperiment)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)

# ============================================================
# STEP 1: Load filtered data from tcga_v2.R
# Make sure you have run tcga_v2.R first
# ============================================================
load("data_filtered.rda")

# ============================================================
# STEP 2: Extract raw counts and create DGEList object
# DGEList is edgeR's container for count data
# ============================================================
counts <- assay(data_filtered, "unstranded") # row=gene, column=sample
condition <- colData(data_filtered)$condition

dge <- DGEList(counts = counts, group = condition)

# ============================================================
# STEP 3: Filter low-count genes
# Keep genes sufficiently expressed in at least the smallest group size
# ============================================================
samples <- table(condition) # 133 Normal, 571 LuminalA, 197 Basal
cat(sprintf("Total samples: %d\n", length(condition))) # 881
cat(sprintf("Normal: %d\nLuminalA: %d\nBasal: %d\n",
            samples["Normal"],
            samples["LuminalA"],
            samples["Basal"]))
min_samples <- min(samples) 
# a gene needs at least 10 counts in roughly 113 samples 
keep_genes <- filterByExpr(dge, group = condition, min.count = 10)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]
cat(sprintf("Genes remaining after filtering: %d\n", nrow(dge))) # 26199

# ============================================================
# STEP 4: TMM normalization
# Corrects for composition bias between samples
# ============================================================
# ensures the counts are comparable across samples before any modeling
dge <- calcNormFactors(dge, method = "TMM")

# ============================================================
# STEP 5: Build design matrix
# Normal is the reference level (set in tcga_v2.R)
# ============================================================
design <- model.matrix(~ condition)
colnames(design)  # verify: intercept=Normal, conditionLuminalA, conditionBasal

# ============================================================
# STEP 6: voom transformation
# Converts counts to logCPM with precision weights
# ============================================================
v <- voom(dge, design, plot = TRUE)

# ============================================================
# STEP 7: Fit linear model and apply empirical Bayes
# Emperical Bayes stabilizes gene-wise variance estimates
# ============================================================
fit <- lmFit(v, design)
fit <- eBayes(fit) 


# ============================================================
# STEP 8: Extract DE
# ============================================================
get_de <- function(fit, coef_name, label) {
  
  res <- topTable(fit, coef = coef_name,
                  number = Inf,
                  sort.by = "none")
  
  res$ENSEMBL <- sub("\\..*", "", rownames(res))
  
  # dynamic threshold
  thr <- mean(abs(res$logFC)) + 2*sd(abs(res$logFC))
  
  sig <- res[abs(res$logFC) > thr & res$P.Value < 0.05, ]
  
  cat("\n", label,
      "\nTotal DE:", nrow(sig),
      "\nUp:", sum(sig$logFC > 0),
      "\nDown:", sum(sig$logFC < 0), "\n")
  
  list(res = res, sig = sig, threshold = thr)
}

lumA <- get_de(fit, "conditionLuminalA", "LuminalA vs Normal") # 1409 Total, 402 Up, 1007 Down
basal <- get_de(fit, "conditionBasal", "Basal vs Normal") # 1349 Total, 421 Up, 928 Down


# ============================================================
# STEP 9: Annotate Data
# ============================================================

annotate <- function(df) {
  
  map <- bitr(df$ENSEMBL,
              fromType = "ENSEMBL",
              toType = c("SYMBOL", "GENENAME"),
              OrgDb = org.Hs.eg.db)
  
  df <- left_join(df, map, by = "ENSEMBL")
  
  df
}

# annonate data
res_lumA <- annotate(lumA$res) # 22.02% fail to map
res_basal <- annotate(basal$res) # 22.02% fail to map
sig_lumA <- annotate(lumA$sig) # 13.63% fail to map
sig_basal <- annotate(basal$sig) # 13.79% fail to map

# delete na
res_lumA <- res_lumA[!is.na(res_lumA$SYMBOL), ] 
res_basal <- res_basal[!is.na(res_basal$SYMBOL), ]
sig_lumA <- sig_lumA[!is.na(sig_lumA$SYMBOL), ]
sig_basal <- sig_basal[!is.na(sig_basal$SYMBOL), ]

# delete duplicates
# res_lumA <- res_lumA[!duplicated(res_lumA$SYMBOL), ]
# res_basal <- res_basal[!duplicated(res_basal$SYMBOL), ]
# sig_lumA <- sig_lumA[!duplicated(sig_lumA$SYMBOL), ]
# sig_basal <- sig_basal[!duplicated(sig_basal$SYMBOL), ]


# ============================================================
# STEP 10: Directional Analysis
# ============================================================
up_lumA <- sig_lumA$SYMBOL[sig_lumA$logFC > 0]
down_lumA <- sig_lumA$SYMBOL[sig_lumA$logFC < 0]

up_basal <- sig_basal$SYMBOL[sig_basal$logFC > 0]
down_basal <- sig_basal$SYMBOL[sig_basal$logFC < 0]

all_lumA <- union(up_lumA, down_lumA)
all_basal <- union(up_basal, down_basal)

shared <- intersect(all_lumA, all_basal)
unique_lumA <- setdiff(all_lumA, all_basal)
unique_basal <- setdiff(all_basal, all_lumA)

opposite <- union(
  intersect(up_lumA, down_basal),
  intersect(down_lumA, up_basal)
)

cat("Shared:", length(shared), "\n") # 569
cat("Unique LumA:", length(unique_lumA), "\n") # 664
cat("Unique Basal:", length(unique_basal), "\n") # 610
cat("Opposite:", length(opposite), "\n") # 11


# ============================================================
# STEP 11: Save Data
# ============================================================

write.csv(sig_lumA, "DE_LuminalA_vs_Normal.csv", row.names = FALSE)
write.csv(sig_basal, "DE_Basal_vs_Normal.csv", row.names = FALSE)
write.csv(res_lumA, "All_LuminalA_vs_Normal.csv", row.names = FALSE)
write.csv(res_basal, "All_Basal_vs_Normal.csv", row.names = FALSE)
cat("\nResults saved to DE_LuminalA_vs_Normal.csv, DE_Basal_vs_Normal.csv, 
    All_LuminalA_vs_Normal.csv, and All_Basal_vs_Normal.csv\n")

# ============================================================
# STEP 12: Plot Data
# Volcano Plots
# ============================================================

plot_volcano <- function(res, threshold, title, file_name) {
  
  res$group <- "Equal"
  res$group[res$logFC > threshold & res$P.Value < 0.05] <- "Up"
  res$group[res$logFC < -threshold & res$P.Value < 0.05] <- "Down"
  
  colors <- c("Up" = "red", "Equal" = "black", "Down" = "green")
  
  png(file_name, width = 1200, height = 1000, res = 150)
  
  plot(res$logFC,
       -log10(res$P.Value),
       pch = 20,
       col = colors[res$group],
       xlab = "log2 Fold Change",
       ylab = "-log10(P-value)",
       main = title)
  
  abline(h = -log10(0.05), col = "grey70")
  abline(v = c(-threshold, threshold), col = "grey70")
  
  legend("topright",
         legend = names(colors),
         col = colors,
         pch = 20,
         bty = "o")
  
  dev.off()
}

# plots
plot_volcano(res_lumA,
             lumA$threshold,
             "Luminal A vs Normal Volcano Plot",
             "Volcano_LuminalA_vs_Normal.png")

plot_volcano(res_basal,
             basal$threshold,
             "Basal vs Normal Volcano Plot",
             "Volcano_Basal_vs_Normal.png")