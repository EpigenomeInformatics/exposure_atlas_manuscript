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
  # show individual samples as jittered points (seeded for reproducibility)
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

# version without the significance brackets
p_nobrackets <- p
p_nobrackets$layers <- p_nobrackets$layers[
  !vapply(p_nobrackets$layers, function(l) inherits(l$stat, "StatCompareMeans"), logical(1))
]
ggsave(p_nobrackets,
  filename = paste0(fig_dir, "cell_prop_plots_nobrackets.pdf"),
  width = 20, height = 16
)

#####################################################################
# Per-sample Wilcoxon rank-sum tests with Bonferroni correction
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

# Bonferroni across the full family of tests
prop_test_df$p_adj <- p.adjust(prop_test_df$p_value, method = "bonferroni")
prop_test_df <- prop_test_df[order(prop_test_df$p_adj), ]
write.csv(prop_test_df,
  file = paste0(fig_dir, "cell_prop_wilcox_bonferroni.csv"), row.names = FALSE
)
message("Significant cell-type composition changes (Bonferroni p_adj < 0.05):")
print(subset(prop_test_df, p_adj < 0.05))

#####################################################################
# OP sensitivity: samples sharing recorded age + sex may be repeats of one person.
# Average proportions within age x sex x exposure-group blocks, rerun the OP tests.
# Samples without age/sex stay as their own block. Bonferroni over the same family as above.
#####################################################################
repo_dir   <- "/scratch/icbb/igunduz/irem_github/exposure_atlas_manuscript"
suppTables <- file.path(repo_dir, "sample_annots/All_Supplementary_Tables_updated.xlsx")
if (!file.exists(suppTables)) {
  suppTables <- file.path(repo_dir, "sample_annots/All_Supplementary_Tables.xlsx")
}
stopifnot(file.exists(suppTables))

op_meta <- readxl::read_excel(suppTables, sheet = "Table S1") %>%
  dplyr::filter(exposure_type == "OP") %>%
  dplyr::transmute(Sample = arrow_name, Age = as.numeric(Age), Sex = as.character(Sex))

op_df <- df %>%
  dplyr::filter(grepl("^OP_", as.character(sample_exposure_group))) %>%
  dplyr::mutate(Sample = as.character(Sample)) %>%
  dplyr::left_join(op_meta, by = "Sample")
stopifnot(all(unique(op_df$Sample) %in% op_meta$Sample))

op_df <- op_df %>%
  dplyr::mutate(block = ifelse(is.na(Age) | is.na(Sex), Sample,
                               paste(Age, Sex, sample_exposure_group, sep = "_")))

# age/sex blocks that span exposure groups cannot be averaged into one group
cross_grp <- op_df %>%
  dplyr::filter(!is.na(Age), !is.na(Sex)) %>%
  dplyr::distinct(Sample, Age, Sex, sample_exposure_group) %>%
  dplyr::group_by(Age, Sex) %>%
  dplyr::filter(dplyr::n_distinct(sample_exposure_group) > 1) %>%
  dplyr::summarise(groups = paste(sort(unique(as.character(sample_exposure_group))), collapse = ","),
                   n_samples = dplyr::n(), .groups = "drop")
if (nrow(cross_grp) > 0) {
  message("Age/sex blocks spanning OP groups (kept apart by group):")
  print(as.data.frame(cross_grp))
}

op_blocks <- op_df %>%
  dplyr::group_by(ClusterCellTypes, sample_exposure_group, block) %>%
  dplyr::summarise(freq = mean(freq), n_in_block = dplyr::n(), .groups = "drop")

n_samp  <- dplyr::n_distinct(op_df$Sample)
n_block <- dplyr::n_distinct(op_blocks$block)
message(sprintf("OP: %d samples -> %d age/sex/group blocks", n_samp, n_block))
print(table(droplevels(dplyr::distinct(op_blocks, block, sample_exposure_group)$sample_exposure_group)))

op_comparisons <- Filter(function(cmp) all(grepl("^OP_", cmp)), comparisons)
op_block_test <- do.call(rbind, lapply(op_comparisons, function(cmp) {
  do.call(rbind, lapply(unique(as.character(op_blocks$ClusterCellTypes)), function(ct) {
    sub <- op_blocks[as.character(op_blocks$ClusterCellTypes) == ct &
      op_blocks$sample_exposure_group %in% cmp, ]
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

# same Bonferroni family size as the primary table: OP rows swapped for the block-level ones
stopifnot(nrow(op_block_test) == sum(grepl("^OP_", prop_test_df$group1)))
op_block_test$p_adj <- pmin(1, op_block_test$p_value * nrow(prop_test_df))

op_compare <- prop_test_df %>%
  dplyr::filter(grepl("^OP_", group1)) %>%
  dplyr::select(cell_type, group1, group2, n1, n2, p_value, p_adj) %>%
  dplyr::inner_join(op_block_test, by = c("cell_type", "group1", "group2"),
                    suffix = c("_samples", "_blocks")) %>%
  dplyr::mutate(sig_samples = p_adj_samples < 0.05,
                sig_blocks  = p_adj_blocks < 0.05,
                changed     = sig_samples != sig_blocks) %>%
  dplyr::arrange(p_adj_blocks)
stopifnot(nrow(op_compare) == nrow(op_block_test))

write.csv(op_compare,
  file = file.path(fig_dir, "cell_prop_wilcox_bonferroni_OP_agesex_blocks.csv"), row.names = FALSE
)
message(sprintf("OP tests significant (Bonferroni < 0.05): %d per sample, %d per age/sex block; %d changed",
                sum(op_compare$sig_samples), sum(op_compare$sig_blocks), sum(op_compare$changed)))
print(subset(op_compare, sig_samples | sig_blocks), digits = 3)

#####################################################################
# Writes the raw statistics CSV; table assembly is in 12_assemble_supplementary_tables.R
#####################################################################

#####################################################################
