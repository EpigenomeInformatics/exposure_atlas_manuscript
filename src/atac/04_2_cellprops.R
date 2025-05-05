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
fig_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/figures/"

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
  "Influenza_d30" = "#E1861A",
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
  geom_boxplot(outlier.shape = NA, alpha = 1, color = "black") + # Add black outlines
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

#####################################################################
