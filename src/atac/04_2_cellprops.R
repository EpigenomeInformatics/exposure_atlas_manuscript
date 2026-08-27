#!/usr/bin/env Rscript

#####################################################################
# 04_2_cellprops.R
# created on 2024-08-24 by Irem Gunduz
# Plot cell type proportions across samples
#####################################################################

# load libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(ggpubr)
  library(ArchR)
})
set.seed(12)

outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/"
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)
fig_dir <- "/scratch/icbb/igunduz/irem_github/exposure_atlas_manuscript/figures/"

if (!dir.exists(fig_dir)) {
  dir.create(fig_dir, recursive = TRUE)
}

# Convert `project@cellColData` to a data frame and filter out unnecessary groups
cell_data <- as.data.frame(project@cellColData) %>%
  dplyr::filter(!sample_exposure_group %in% c("BA_na", "BA_vac")) %>%
  dplyr::mutate(NewSample = paste0(sample_exposure_group, "_", Sample))

# Define variables for grouping and coloring
each.pt <- "Sample"
group_by <- "sample_exposure_group"
color_by <- "ClusterCellTypes"

# Calculate proportions for each cell type by sample
fq <- prop.table(table(cell_data$ClusterCellTypes, cell_data$Sample), 2)
df <- reshape2::melt(fq, value.name = "freq", varnames = c("cell.type", each.pt))
colnames(df) <- c("ClusterCellTypes", "Sample", "freq")

# Merge with metadata
meta.include <- unique(c(each.pt, group_by))
ei <- unique(cell_data[, meta.include])
df <- merge(df, ei, by = each.pt)
df <- cbind(df, null.group = paste("1"))
df[, each.pt] <- as.factor(df[, each.pt])

# Define an order for `sample_exposure_group`
order <- c(
  "C19_ctrl", "C19_mild", "C19_mod", "C19_sev",
  "HIV_ctrl", "HIV_acu", "HIV_chr",
  "Influenza_ctrl", "Influenza_d3", "Influenza_d6", "Influenza_d30",
  "OP_low", "OP_med", "OP_high"
)

# Set factor levels for the correct ordering
df$sample_exposure_group <- factor(df$sample_exposure_group, levels = order)

# Define specific comparisons for each cohort (control vs other conditions)
comparisons <- list(
  c("C19_ctrl", "C19_mild"),
  c("C19_ctrl", "C19_mod"),
  c("C19_ctrl", "C19_sev"),
  c("HIV_ctrl", "HIV_acu"),
  c("HIV_ctrl", "HIV_chr"),
  c("Influenza_ctrl", "Influenza_d3"),
  c("Influenza_ctrl", "Influenza_d6"),
  c("Influenza_ctrl", "Influenza_d30"),
  c("OP_low", "OP_med"),
  c("OP_low", "OP_high")
)

# Define the custom color palette based on your legend
custom_colors <- c(
  "C19_Ctrl" = "#6CA6CD",
  "C19_mild" = "#8FBC8F",
  "C19_mod" = "#2C948F",
  "C19_sev" = "#006400",
  "HIV_acu" = "#893368",
  "HIV_chr" = "#825CA3",
  "HIV_ctrl" = "#4F619D",
  "Influenza_ctrl" = "#4F619D",
  "Influenza_d3" = "#FFD34E",
  "Influenza_d6" = "#EBB332",
  # day 28 is the true final timepoint; "d30" is a legacy mislabel
  "Influenza_d30" = "#E1861A",
  "Influenza_d28" = "#E1861A",
  "OP_low" = "#BADBF4",
  "OP_med" = "#94A9D3",
  "OP_high" = "#16528A"
)
p <- ggplot(df, aes(x = sample_exposure_group, y = freq, fill = sample_exposure_group)) +
  labs(y = "Proportion of PBMCs (%)", x = "Exposure") +
  theme_classic() + # Clean background
  theme(
    strip.background = element_blank(), # Clean strip background
    strip.text = element_text(face = "bold"),
    axis.ticks.x = element_blank(),
    axis.text = element_text(color = "black"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    aspect.ratio = 1,
    text = element_text(size = 20),
    legend.position = "bottom",
    legend.direction = "vertical",
    panel.border = element_rect(color = "black", fill = NA) # Panel border
  ) +
  facet_wrap("ClusterCellTypes", scales = "free_y", ncol = 4) +
  # outlier.shape = NA so outliers are not drawn twice once the points are added
  geom_boxplot(outlier.shape = NA, alpha = 1, color = "black") +
  # Individual samples shown as points (reviewer 1, minor comment 1): with only a
  # few donors per group the distribution is not meaningful on its own, so every
  # sample is displayed. Seeded jitter keeps the figure reproducible.
  geom_point(
    position = position_jitter(width = 0.18, height = 0, seed = 12),
    size = 1.1, alpha = 0.75, colour = "black", show.legend = FALSE
  ) +
  scale_fill_manual(values = custom_colors) + # Vibrant colors
  stat_compare_means(
    comparisons = comparisons,
    p.adjust.method = "bonferroni",
    method = "wilcox.test",
    exact = FALSE,
    label = "p.signif",
    tip.length = 0.01,
    step.increase = 0.1
  )

# Save the plot
ggsave(p, filename = paste0(fig_dir, "cell_prop_plots.pdf"), width = 20, height = 20)

# Version without the significance brackets. The brackets add roughly a third of
# the panel height and, as noted below, show RAW p-values; the adjusted values
# are reported in cell_prop_wilcox_bonferroni.csv. Use this version if space in
# the supplementary figure is tight.
p_nobrackets <- p
p_nobrackets$layers <- p_nobrackets$layers[
  !vapply(p_nobrackets$layers, function(l) inherits(l$stat, "StatCompareMeans"), logical(1))
]
ggsave(p_nobrackets,
  filename = paste0(fig_dir, "cell_prop_plots_nobrackets.pdf"),
  width = 20, height = 16
)

#####################################################################
# Explicit per-sample Wilcoxon rank-sum tests with Bonferroni correction
# (R2 / R3.2). NOTE: stat_compare_means(comparisons = ...) draws each bracket's
# RAW pairwise p-value; its p.adjust.method is ignored for bracket comparisons.
# We therefore compute the adjusted p-values we report here, at the sample
# (donor) level, so the reported padj values are reproducible.
#####################################################################
prop_test_df <- do.call(rbind, lapply(comparisons, function(cmp) {
  do.call(rbind, lapply(unique(as.character(df$ClusterCellTypes)), function(ct) {
    sub <- df[as.character(df$ClusterCellTypes) == ct &
      df$sample_exposure_group %in% cmp, ]
    g1 <- sub$freq[sub$sample_exposure_group == cmp[1]]
    g2 <- sub$freq[sub$sample_exposure_group == cmp[2]]
    if (length(g1) < 1 || length(g2) < 1) {
      return(NULL)
    }
    wt <- suppressWarnings(wilcox.test(g1, g2, exact = FALSE))
    data.frame(
      cell_type = ct, group1 = cmp[1], group2 = cmp[2],
      n1 = length(g1), n2 = length(g2),
      median1 = median(g1), median2 = median(g2),
      p_value = wt$p.value, stringsAsFactors = FALSE
    )
  }))
}))

# Bonferroni across the full family of tests. To match a narrower family (e.g.
# per cohort or per cell type), group_by before p.adjust instead.
prop_test_df$p_adj <- p.adjust(prop_test_df$p_value, method = "bonferroni")
prop_test_df <- prop_test_df[order(prop_test_df$p_adj), ]
write.csv(prop_test_df,
  file = paste0(fig_dir, "cell_prop_wilcox_bonferroni.csv"), row.names = FALSE
)
message("Significant cell-type composition changes (Bonferroni p_adj < 0.05):")
print(subset(prop_test_df, p_adj < 0.05))

#####################################################################
# NOTE: this script now only writes the raw statistics CSV above
# (cell_prop_wilcox_bonferroni.csv). Insertion of the composition table into the
# supplementary workbook and ALL supplementary-table numbering are handled in one
# place — 12_assemble_supplementary_tables.R — so there is a single source of
# truth for the numbering. Run that script last, after this one.
#####################################################################

#####################################################################
