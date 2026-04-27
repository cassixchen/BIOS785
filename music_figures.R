# ============================================================
# MuSiC deconvolution figures
#   Fig 1  (music_purity_validation_{LumA,Basal}.png) - purity validation scatter
#   Fig 2  (music_stacked_composition.png)            - mean composition per subtype
#
# Need music_proportions.csv (from music_deconvolution.R)
# ============================================================
setwd("/Users/amyfliu/Documents/bio785_final")

library(ggplot2)
library(dplyr)
library(tidyr)
library(TCGAbiolinks)

# ---- Load saved MuSiC output ----
props <- read.csv("music_proportions.csv", check.names = FALSE,
                  stringsAsFactors = FALSE)

celltypes_to_use <- c(
  "Cancer Epithelial", "Normal Epithelial", "CAFs",
  "Endothelial", "T-cells", "B-cells",
  "Myeloid", "Plasmablasts", "PVL"
)

# Subtype factor order (LumA, Basal) - match throughout
props$condition <- factor(props$condition, levels = c("LuminalA", "Basal"))

# ============================================================
# FIG 1 -- Purity validation scatter 
# ============================================================
purity <- TCGAbiolinks::Tumor.purity

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

df <- merge(props, purity, by = "sample_id")
df$stromal  <- df[["CAFs"]] + df[["Endothelial"]] + df[["PVL"]]

# Purity columns in TCGAbiolinks::Tumor.purity are character with comma decimals
# in some locales - coerce safely.
df$ABSOLUTE <- as.numeric(sub(",", ".", df$ABSOLUTE))
df$ESTIMATE <- as.numeric(sub(",", ".", df$ESTIMATE))

scatter_df <- df %>%
  select(barcode, condition,
         `Cancer Epithelial`, stromal,
         ABSOLUTE, ESTIMATE) %>%
  pivot_longer(cols = c(`Cancer Epithelial`, stromal),
               names_to = "compartment", values_to = "proportion") %>%
  pivot_longer(cols = c(ABSOLUTE, ESTIMATE),
               names_to = "purity_metric", values_to = "purity_score") %>%
  filter(!is.na(purity_score), !is.na(proportion))

scatter_df$compartment <- recode(scatter_df$compartment,
                                 "stromal" = "Stromal (CAFs + Endo + PVL)")
scatter_df$compartment <- factor(scatter_df$compartment,
                                 levels = c("Cancer Epithelial",
                                            "Stromal (CAFs + Endo + PVL)"))

# Per-facet Spearman rho + n for annotation
# Group by subtype since panels are now split LumA / Basal
cor_labels <- scatter_df %>%
  group_by(compartment, purity_metric, condition) %>%
  summarise(
    rho = cor(proportion, purity_score, method = "spearman",
              use = "complete.obs"),
    n   = sum(!is.na(proportion) & !is.na(purity_score)),
    .groups = "drop"
  ) %>%
  mutate(label = sprintf("rho = %.2f\nn = %d", rho, n))

# Build one plot per subtype (easier to read than a single 2x4 grid!)
build_purity_plot <- function(subtype_name, point_color) {
  sub_data   <- scatter_df %>% filter(condition == subtype_name)
  sub_labels <- cor_labels %>% filter(condition == subtype_name)
  n_tumors   <- length(unique(sub_data$barcode))

  ggplot(sub_data, aes(x = purity_score, y = proportion)) +
    geom_point(alpha = 0.45, size = 1.0, color = point_color) +
    geom_smooth(method = "lm", se = FALSE,
                color = "black", linewidth = 0.6,
                fullrange = FALSE) +
    geom_text(data = sub_labels,
              aes(x = -Inf, y = Inf, label = label),
              hjust = -0.08, vjust = 1.25,
              inherit.aes = FALSE, size = 3.6, lineheight = 0.95) +
    facet_grid(compartment ~ purity_metric, scales = "free_x") +
    scale_y_continuous(limits = c(0, NA),
                       expand = expansion(mult = c(0.01, 0.08))) +
    labs(x = "TCGA purity score",
         y = "MuSiC-estimated compartment proportion",
         title = sprintf(
           "MuSiC proportions vs TCGA tumor-purity estimates - %s",
           subtype_name),
         subtitle = sprintf(
           "n = %d tumors | Spearman correlations per panel", n_tumors)) +
    theme_bw(base_size = 11) +
    theme(strip.background = element_rect(fill = "grey95"),
          plot.title = element_text(face = "bold"))
}

p1_lumA  <- build_purity_plot("LuminalA", "#2c7bb6")
p1_basal <- build_purity_plot("Basal",    "#d7191c")

# Export figures
ggsave("music_purity_validation_LumA.png",  p1_lumA,
       width = 6, height = 5, dpi = 300)
ggsave("music_purity_validation_Basal.png", p1_basal,
       width = 6, height = 5, dpi = 300)

# ============================================================
# FIG 2 -- Mean cell-type composition per PAM50 subtype (stacked bar)
#   Biological story: LumA is epithelial-dominant; Basal has a
#   larger stromal / immune share; 
#   Note: Proportions within each bar sum to 1 because MuSiC 
#   normalizes across the 9 selected cell types.
# ============================================================

# Stacking order in the bars (bottom -> top): epithelial -> stromal -> immune.
# ggplot default position = "stack" places factor level 1 at the bottom.
celltype_stack_order <- c(
  "Cancer Epithelial", "Normal Epithelial",
  "CAFs", "Endothelial", "PVL",
  "T-cells", "B-cells", "Myeloid", "Plasmablasts"
)

# Legend order includes 3 "header" rows that act as section titles for each
# lineage group. They are added as fake factor levels with white tile fills,
# so they appear as labeled section breaks in the legend without showing as
# colored entries.
celltype_legend_order <- c(
  "Epithelial",
  "Cancer Epithelial", "Normal Epithelial",
  "Stromal",
  "CAFs", "Endothelial", "PVL",
  "Immune",
  "T-cells", "B-cells", "Myeloid", "Plasmablasts"
)

celltype_colors <- c(
  "Epithelial"        = "white",     # header (invisible against white legend bg)
  "Cancer Epithelial" = "#B2182B",   # dark red    - tumor epithelial
  "Normal Epithelial" = "#F4A582",   # salmon      - normal epithelial
  "Stromal"           = "white",     # header
  "CAFs"              = "#762A83",   # purple      - stroma
  "Endothelial"       = "#FDB863",   # orange      - stroma
  "PVL"               = "#8C510A",   # brown       - stroma
  "Immune"            = "white",     # header
  "T-cells"           = "#2166AC",   # blue        - immune
  "B-cells"           = "#92C5DE",   # light blue  - immune
  "Myeloid"           = "#1B7837",   # dark green  - immune
  "Plasmablasts"      = "#A6DBA0"    # light green - immune
)

# Bold + italic the header labels via plotmath; indent the cell-type labels
# with leading whitespace so the lineage hierarchy is visible at a glance.
celltype_legend_labels <- setNames(celltype_legend_order, celltype_legend_order)
header_idx <- celltype_legend_order %in% c("Epithelial", "Stromal", "Immune")
celltype_legend_labels[header_idx]  <- toupper(celltype_legend_labels[header_idx])
celltype_legend_labels[!header_idx] <- paste0("   ",
                                              celltype_legend_labels[!header_idx])

mean_comp <- props %>%
  select(condition, all_of(celltypes_to_use)) %>%
  pivot_longer(-condition, names_to = "celltype", values_to = "prop") %>%
  group_by(condition, celltype) %>%
  summarise(mean_prop = mean(prop), .groups = "drop") %>%
  mutate(celltype = factor(celltype, levels = celltype_stack_order)) %>%
  arrange(condition, celltype) %>%
  group_by(condition) %>%
  # Label position: midpoint of each stacked segment (measured from 0 up).
  mutate(label_pos = cumsum(mean_prop) - mean_prop / 2,
         label     = ifelse(mean_prop >= 0.03,
                            sprintf("%.0f%%", mean_prop * 100), "")) %>%
  ungroup()

# Promote the celltype factor to include the 3 header levels so the legend
# can show them. The data still has only 9 actual rows per subtype (no rows
# for the header levels), so no extra tiles are drawn in the bars.
mean_comp$celltype <- factor(as.character(mean_comp$celltype),
                             levels = celltype_legend_order)

# Add sample size under each subtype label
n_per_subtype <- props %>% count(condition)
subtype_labels <- setNames(
  sprintf("%s\n(n = %d)", n_per_subtype$condition, n_per_subtype$n),
  as.character(n_per_subtype$condition)
)

p2 <- ggplot(mean_comp,
             aes(x = condition, y = mean_prop, fill = celltype)) +
  geom_col(width = 0.55, color = "white", linewidth = 0.3) +
  geom_col(position = position_stack(reverse = T)) +
  geom_text(aes(y = label_pos, label = label),
            color = "white", size = 3.5, fontface = "bold") +
  scale_x_discrete(labels = subtype_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02)),
                     labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values  = celltype_colors,
                    breaks  = celltype_legend_order,
                    limits  = celltype_legend_order,
                    labels  = celltype_legend_labels,
                    drop    = FALSE) +
  labs(x = NULL, y = "Mean estimated proportion",
       fill = "Cell type",
       title = "MuSiC cell-type composition - TCGA-BRCA by subtype",
       subtitle = "Mean across tumors in each subtype") +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold"),
    axis.text.x        = element_text(size = 11),
    legend.key         = element_rect(color = NA, fill = NA),
    legend.text        = element_text(size = 9),
    legend.spacing.y   = unit(0.05, "cm"),
    legend.title=element_blank()
  )

# Export figure
ggsave("music_stacked_composition.png", p2,
       width = 6, height = 5, dpi = 300)
