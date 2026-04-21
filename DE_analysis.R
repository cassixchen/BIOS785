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
cat(sprintf("Genes remaining after filtering: %d\n", nrow(dge))) #26199

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
# STEP 8: Extract results using dynamic threshold
# Matching paper's criteria:
# |logFC| > mean(|logFC|) + 2*sd(|logFC|) AND P value < 0.05
# ============================================================

# --- a) LuminalA vs Normal ---
res_lumA <- topTable(fit,
                     coef = "conditionLuminalA",
                     number = Inf,
                     sort.by = "none")

# Dynamic logFC threshold
logfc_lumA <- abs(res_lumA$logFC)
threshold_lumA <- mean(logfc_lumA) + 2 * sd(logfc_lumA)
cat(sprintf("\nLumA logFC threshold: %.3f\n", threshold_lumA))

sig_lumA <- res_lumA[abs(res_lumA$logFC) > threshold_lumA & res_lumA$P.Value < 0.05, ]
sig_lumA <- sig_lumA[order(sig_lumA$P.Value), ]
cat(sprintf("Significant DE genes (LumA vs Normal): %d\n", nrow(sig_lumA)))
cat(sprintf("Upregulated: %d\n",   sum(sig_lumA$logFC > 0)))
cat(sprintf("Downregulated: %d\n", sum(sig_lumA$logFC < 0)))

# --- b) Basal vs Normal ---
res_basal <- topTable(fit,
                      coef = "conditionBasal",
                      number = Inf,
                      sort.by = "none")

# Dynamic logFC threshold
logfc_basal <- abs(res_basal$logFC)
threshold_basal <- mean(logfc_basal) + 2 * sd(logfc_basal)
cat(sprintf("\nBasal logFC threshold: %.3f\n", threshold_basal))

sig_basal <- res_basal[abs(res_basal$logFC) > threshold_basal & res_basal$P.Value < 0.05, ]
sig_basal <- sig_basal[order(sig_basal$P.Value), ]
cat(sprintf("Significant DE genes (Basal vs Normal): %d\n", nrow(sig_basal)))
cat(sprintf("Upregulated: %d\n",   sum(sig_basal$logFC > 0)))
cat(sprintf("Downregulated: %d\n", sum(sig_basal$logFC < 0)))

# ============================================================
# STEP 9: clusterProfiler and org.H.eg.db
# Convert IDs into gene names
# ============================================================

# --- a) LuminalA vs Normal ---

sig_lumA$gene_id <- sub("\\..*", "", rownames(sig_lumA)) # remove Ensembl version numbers

# 13.63% fail to map
gene_map_lumA <- bitr(sig_lumA$gene_id,
                      fromType = "ENSEMBL",
                      toType = "SYMBOL",
                      OrgDb = org.Hs.eg.db)

sig_lumA <- merge(sig_lumA, gene_map_lumA,
                        by.x = "gene_id",
                        by.y = "ENSEMBL",
                        all.x = TRUE)


# --- b) Basal vs Normal ---

sig_basal$gene_id <- sub("\\..*", "", rownames(sig_basal)) # remove Ensembl version numbers

# 13.79% fail to map
gene_map_basal<- bitr(sig_basal$gene_id,
                      fromType = "ENSEMBL",
                      toType = "SYMBOL",
                      OrgDb = org.Hs.eg.db)

sig_basal <- merge(sig_basal, gene_map_basal,
                  by.x = "gene_id",
                  by.y = "ENSEMBL",
                  all.x = TRUE)

# ============================================================
# STEP 10: Volcano plots
# Visualize DE and non DE genes
# ============================================================

# --- a) LuminalA vs Normal ---

# get log fold change group
res_lumA$group <- "Equal"
res_lumA$group[res_lumA$logFC > threshold_lumA & res_lumA$P.Value < 0.05] <- "Up"
res_lumA$group[res_lumA$logFC < -threshold_lumA & res_lumA$P.Value < 0.05] <- "Down"

colors <- c("Up" = "red", "Equal" = "black", "Down" = "green")

# png
png("Volcano_LuminalA_vs_Normal.png", width = 1200, height = 1000, res = 150)

# scatterplot of log fold change vs statistical significance
plot(res_lumA$logFC,
     -log10(res_lumA$P.Value),
     pch = 20,
     col = colors[res_lumA$group],
     xlab = "log2(FC)",
     ylab = "-log10(P-value)",
     main = "Luminal A vs Normal Volcano Plot")

# significance thresholds
abline(h = -log10(0.05), col = "grey70")
abline(v = c(-threshold_lumA, threshold_lumA), col = "grey70")

# legend
legend("topright",
       legend = names(colors),
       col = colors,
       pch = 20,
       bty = "o")

dev.off()
cat("\nResults saved to Volcano_LuminalA_vs_Normal_limma.png\n")

# --- b) Basal vs Normal ---

res_basal$group <- "Equal"
res_basal$group[res_basal$logFC > threshold_basal & res_basal$P.Value < 0.05] <- "Up"
res_basal$group[res_basal$logFC < -threshold_basal & res_basal$P.Value < 0.05] <- "Down"

colors <- c("Up" = "red", "Equal" = "black", "Down" = "green")

# png
png("Volcano_Basal_vs_Normal.png", width = 1200, height = 1000, res = 150)

# scatterplot of log fold change vs statistical significance
plot(res_basal$logFC,
     -log10(res_basal$P.Value),
     pch = 20,
     col = colors[res_basal$group],
     xlab = "log2(FC)",
     ylab = "-log10(P-value)",
     main = "Basal vs Normal Volcano Plot")

# significance thresholds
abline(h = -log10(0.05), col = "grey70")
abline(v = c(-threshold_basal, threshold_basal), col = "grey70")

# legend
legend("topright",
       legend = names(colors),
       col = colors,
       pch = 20,
       bty = "o")

dev.off()
cat("\nResults saved to Volcano_Basal_vs_Normal.png\n")

# ============================================================
# STEP 11: dplyr or alternative
# Identify common, unique, opposite DE genes
# ============================================================

# delete genes with NA symbol
sig_lumA <- sig_lumA[!is.na(sig_lumA$SYMBOL), ]
sig_basal <- sig_basal[!is.na(sig_basal$SYMBOL), ]

# get direction of DE genes
up_lumA   <- subset(sig_lumA, logFC > 0)$SYMBOL
down_lumA <- subset(sig_lumA, logFC < 0)$SYMBOL
up_basal   <- subset(sig_basal, logFC > 0)$SYMBOL
down_basal <- subset(sig_basal, logFC < 0)$SYMBOL

# DE genes with direction
all_lumA <- union(up_lumA, down_lumA)
all_basal <- union(up_basal, down_basal)

# common DE genes
shared_genes <- intersect(all_lumA, all_basal)
cat(sprintf("Shared DE genes between LumA and Basal: %d\n", length(shared_genes))) # 569

# unique DE genes
unique_lumA <- setdiff(all_lumA, all_basal)
cat(sprintf("Unique DE genes to LumA: %d\n", length(unique_lumA))) # 664 
unique_basal <- setdiff(all_basal, all_lumA)
cat(sprintf("Unique DE genes to Basal: %d\n", length(unique_basal))) # 610

# opposite DE genes
opposite_genes <- union(
  intersect(up_lumA, down_basal),
  intersect(down_lumA, up_basal)
)
cat(sprintf("Opposite DE genes between LumA and Basal: %d\n", length(opposite_genes))) # 11 


# ============================================================
# STEP 12: Save DE results
# ============================================================
write.csv(sig_lumA,  file = "DE_LuminalA_vs_Normal.csv")
write.csv(sig_basal, file = "DE_Basal_vs_Normal.csv")
cat("\nResults saved to DE_LuminalA_vs_Normal.csv and DE_Basal_vs_Normal.csv\n")