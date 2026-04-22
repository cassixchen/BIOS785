setwd("/Users/cassixchen/Desktop/unc/BIOS785/Final")

library(clusterProfiler)  
library(org.Hs.eg.db)
library(tidyverse)
library(enrichplot)
library(patchwork)

# ============================================================
# STEP 1: Load data, get gene and universe
# Make sure you have run DE_analysis.R first
# ============================================================

# get entrez column
get_entrez <- function(ensembl_vec) {
  ids <- bitr(ensembl_vec,
              fromType = "ENSEMBL",
              toType = "ENTREZID",
              OrgDb = org.Hs.eg.db)
  
  unique(ids$ENTREZID)
}

# get gene and universe
get_gene_universe <- function(all_file, sig_file) {
  
  all_df <- read.csv(all_file)
  sig_df <- read.csv(sig_file)
  
  universe <- get_entrez(all_df$ENSEMBL)
  gene <- get_entrez(sig_df$ENSEMBL)
  
  list(gene = gene, universe = universe)
}


# ---- basal ----
sets_basal <- get_gene_universe(
  "All_Basal_vs_Normal.csv",
  "DE_Basal_vs_Normal.csv"
)

gene_basal <- sets_basal$gene
universe_basal <- sets_basal$universe

# ---- lumA ----
sets_lumA <- get_gene_universe(
  "All_LuminalA_vs_Normal.csv",
  "DE_LuminalA_vs_Normal.csv"
)

gene_lumA <- sets_lumA$gene
universe_lumA <- sets_lumA$universe

# ============================================================
# STEP 2: Run GO and KEGG
# ============================================================

# run GO enrichment analysis
run_GO <- function(gene, universe, ont) {
  
  ego <- enrichGO(
    gene = gene,
    universe = universe,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ont,
    pvalueCutoff = 0.05
  )
  
  return(ego)
}

# run KEGG path analysis
run_KEGG <- function(gene, universe = NULL, organism = "hsa") {
  
  ekegg <- enrichKEGG(
    gene = gene,
    universe = universe,
    organism = organism,
    pvalueCutoff = 0.05
  )
  
  return(ekegg)
}

# ---- basal ----

bp_basal  <- run_GO(gene_basal, universe_basal, ont = "BP")
cc_basal  <- run_GO(gene_basal, universe_basal, ont = "CC")
mf_basal  <- run_GO(gene_basal, universe_basal, ont = "MF")
kegg_basal <- run_KEGG(gene_basal, universe_basal)

# ---- lumA ----

bp_lumA  <- run_GO(gene_lumA, universe_lumA, ont = "BP")
cc_lumA  <- run_GO(gene_lumA, universe_lumA, ont = "CC")
mf_lumA  <- run_GO(gene_lumA, universe_lumA, ont = "MF")
kegg_lumA <- run_KEGG(gene_lumA, universe_lumA)


# ============================================================
# STEP 3: Plot GO and KEGG
# ============================================================

make_dotplot <- function(res, title, showCategory = 10) {
  dotplot(res, showCategory = showCategory) +
    ggtitle(title)
}

plot_enrichment_panel <- function(bp, cc, mf, kegg,
                                  file_name, main_title) {
  p1 <- make_dotplot(bp, "DEGs: Biological Process")
  p2 <- make_dotplot(cc, "DEGs: Cell Component")
  p3 <- make_dotplot(mf, "DEGs: Molecular Function")
  p4 <- make_dotplot(kegg, "KEGG pathway of DEGs")
  
  png(file_name, width = 3000, height = 2000, res = 200)
  
  print(
    (p1 | p2) / (p3 | p4) +
      plot_annotation(
        title = main_title,
        theme = theme(
          plot.title = element_text(
            size = 20,
            face = "bold",
            hjust = 0.5   # 👈 center it
          )
        )
      )
  )
  
  dev.off()
}


# ---- basal ----
plot_enrichment_panel(
  bp_basal,
  cc_basal,
  mf_basal,
  kegg_basal,
  file_name = "Enrichment_Basal_vs_Normal.png",
  main_title = "Basal vs Normal Enrichment Plots"
)

# ---- basal ----

plot_enrichment_panel(
  bp_lumA,
  cc_lumA,
  mf_lumA,
  kegg_lumA,
  file_name = "Enrichment_LuminalA_vs_Normal.png",
  main_title = "Luminal A vs Normal Enrichment Plots"
)
