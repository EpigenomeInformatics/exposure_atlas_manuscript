#!/usr/bin/env Rscript

#####################################################################
# 10_zdiff.R
# created on 07-05-2025 by Irem Gunduz
# Z-score difference between T-cell subsets using mTFR and cVAR
#####################################################################

get_groupname <- function(x) {
  return(unlist(strsplit(x, split = "_"))[1])
}
source("/icbb/projects/igunduz/exposure_atlas_manuscript/src/utils/mtfr_plots_helpers.R")

cell_type_colors <- c(
  # "Other-cell" = "#CCCCCC",
  "Monocyte" = "#CC4C02",
  "Th-Naive" = "#C7E9B4",
  "Tc-Naive" = "#888FB5",
  "Th-Mem" = "#41B6C4",
  "Tc-Mem" = "#4292C6",
  "NK-cell" = "#A65628",
  "B-cell" = "#AE017E"
)

## Load Libraries
set.seed(42)
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(methylTFR)
  library(muLogR)
  library(gplots)
  library(muRtools)
  library(ComplexHeatmap)
  library(factoextra)
  library(ggfortify)
  library(ArchR)
  library(chromVAR)
})

# Set the paths
plot_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/figures/"
# plot_dir <- "/icbb/projects/igunduz/Figure_4_040425"
if (!dir.exists(plot_dir)) dir.create(plot_dir)
motifset <- "jaspar2020_distal"
sannot_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/sample_annots/"

#####################################################################################
# Differential analysis for T-cell subsets using full-pseudobulks
#####################################################################################
mtfr_dir <- "/icbb/projects/igunduz/mTFR_bias_fix_v3/all_pseudobulks_310824/jaspar2020_distal/"
result <- list.files(mtfr_dir, pattern = "deviations.RDS", full.names = TRUE)

# Loop through combinations and read the data
mtfr_devs <- list()
for (x in seq_along(result)) {
  path <- result[x]
  cur_dev <- readRDS(path)
  cur_dev <- methylTFR::deviations(cur_dev)
  mtfr_devs[[x]] <- as.matrix(cur_dev)
}
mtfr_devs <- do.call(base::cbind, mtfr_devs)
# tdf <- as.data.frame(t(mtfr_devs))
full_groups <- sub(".bedGraph", "", unlist(lapply(FUN = get_groupname, X = colnames(mtfr_devs))))
# tdf$cell_type <- groups

# Subset the data for Th-naive and Th-mem groups
deviations_matrix <- mtfr_devs[, full_groups %in% c("Th-Naive", "Th-Mem")]
groups <- sub(".bedGraph", "", unlist(lapply(FUN = get_groupname, X = colnames(deviations_matrix))))

# Run the differential deviation test
th_result <- differential_deviation_test(
  deviations = deviations_matrix,
  groups = groups,
  alternative = "two.sided",
  parametric = TRUE,
  padjMethod = "BH"
)
saveRDS(th_result, file.path(sannot_dir, "diff_res_TH_two.sided.rds"))

# Do the same for Tc-Naive and Tc-Mem groups
deviations_matrix <- mtfr_devs[, full_groups %in% c("Tc-Naive", "Tc-Mem")]
groups <- sub(".bedGraph", "", unlist(lapply(FUN = get_groupname, X = colnames(deviations_matrix))))

# Run the differential deviation test
tc_result <- differential_deviation_test(
  deviations = deviations_matrix,
  groups = groups,
  alternative = "two.sided",
  parametric = TRUE,
  padjMethod = "BH"
)
saveRDS(tc_result, file.path(sannot_dir, "diff_res_TC_two.sided.rds"))

######################################################################################
# Compute z-diff for Th-Naive and Th-Mem, and Tc-Naive and Tc-Mem
######################################################################################

mtfr_dir <- "/icbb/projects/igunduz/mTFR_bias_fix_v3/mTFR_310824/jaspar2020_distal/"
result <- list.files(mtfr_dir, pattern = "deviations.RDS", full.names = TRUE)
# Loop through combinations and read the data
mtfr_devs <- list()
for (x in seq_along(result)) {
  path <- result[x]
  cur_dev <- readRDS(path)
  cur_dev <- methylTFR::deviations(cur_dev)
  mtfr_devs[[x]] <- as.matrix(cur_dev)
}
mtfr_devs <- do.call(base::cbind, mtfr_devs)
mtfr_devs <- methylTFR:::computeRowZScore(mtfr_devs)
tdf <- as.data.frame(t(mtfr_devs))

# Get the group names
groups <- unlist(lapply(FUN = get_groupname, X = colnames(mtfr_devs)))
groups <- sub(".bedGraph", "", groups)
tdf$cell_type <- groups

# Select Th-Naive and Th-Mem cells
tdf <- tdf %>% filter(cell_type %in% c("Th-Naive", "Th-Mem"))

# Convert tdf to row as motifs, columns as samples
tdf_transposed <- t(as.matrix(tdf[, -ncol(tdf)]))
tdf_transposed <- as.data.frame(tdf_transposed)

# Get the cell type groups
groups <- ifelse(grepl("Th-Naive", colnames(tdf_transposed)), "Th-Naive",
  ifelse(grepl("Th-Mem", colnames(tdf_transposed)), "Th-Mem", NA)
)


diff_res <- readRDS(file.path(sannot_dir, "diff_res_TH_two.sided.rds"))
diff_res$isDiff <- ifelse(diff_res$p_value_adjusted < 0.05, "Diff", "NotDiff")
# saveRDS(diff_res, file.path(plot_dir, "diff_res_TH_two.sided.rds"))

group_diff <- tdf_transposed[, groups == "Th-Naive"] - tdf_transposed[, groups == "Th-Mem"]
diff_df <- data.frame(diff = group_diff)
rownames(diff_df) <- rownames(tdf_transposed)
saveRDS(diff_df, file.path(sannot_dir, "group_diff_thmem_mtfr.rds"))

# Do the same for Tc-Naive and Tc-Mem cells
tdf <- as.data.frame(t(mtfr_devs))
groups <- sub(".bedGraph", "", unlist(lapply(FUN = get_groupname, X = colnames(mtfr_devs))))
tdf$cell_type <- groups
tdf <- tdf %>% filter(cell_type %in% c("Tc-Naive", "Tc-Mem"))

# Convert tdf to row as motifs, columns as samples
tdf_transposed <- t(as.matrix(tdf[, -ncol(tdf)]))
tdf_transposed <- as.data.frame(tdf_transposed)

# Get the cell type groups
groups <- ifelse(grepl("Tc-Naive", colnames(tdf_transposed)), "Tc-Naive",
  ifelse(grepl("Tc-Mem", colnames(tdf_transposed)), "Tc-Mem", NA)
)

diff_res <- readRDS(file.path(sannot_dir, "diff_res_TC_two.sided.rds"))
diff_res$isDiff <- ifelse(diff_res$p_value_adjusted < 0.05, "Diff", "NotDiff")
group_diff <- tdf_transposed[, groups == "Tc-Naive"] - tdf_transposed[, groups == "Tc-Mem"]
diff_df <- data.frame(diff = group_diff)
rownames(diff_df) <- rownames(tdf_transposed)
saveRDS(diff_df, file.path(sannot_dir, "group_diff_tcmem_mtfr.rds"))

#####################################################################
# Compute chromVAR Z-scores
#####################################################################

# Load archr project
outputDir <- "/icbb/projects/igunduz/archr_300824/icbb/projects/igunduz/archr_project_011023"
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)
# annot_acc <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/cell_matching/annot_acc_cellType.rds")
annot_acc <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/cell_matching/annot_acc_nearest_meth_for_acc.rds")

# Subset the cells
idxSample <- BiocGenerics::which(project$cellNames %in% rownames(annot_acc))
cellsSample <- project$cellNames[idxSample]
project <- project[cellsSample, ]

# Subset T-cells
idxSample <- BiocGenerics::which(project$ClusterCellTypes %in% c("T_mem_CD4", "T_mem_CD8", "T_naive"))
cellsSample <- project$cellNames[idxSample]
project <- project[cellsSample, ]

# Get the chromVAR matrix
pseudo_chrom <- getMatrixFromProject(project, "jaspar2020Matrix")
z_scores <- assays(pseudo_chrom)$z
cell_types <- colData(pseudo_chrom)$cell_types
unique_cell_types <- unique(cell_types)
group_summary <- sapply(unique_cell_types, function(ct) {
  mean_z_scores <- rowMeans(z_scores[, cell_types == ct], na.rm = TRUE)
  return(mean_z_scores)
})
group_summary_df <- as.data.frame(group_summary)


# Compute the rowMean Z-score per group difference for T_mem_CD4 and T_naive
naive_means <- group_summary_df[rownames(diff_df), "T_naive"]
mem_means <- group_summary_df[rownames(diff_df), "T_mem_CD4"]
mem8_means <- group_summary_df[rownames(diff_df), "T_mem_CD8"]

# Group difference for T_mem_CD4 and T_naive
group_diff <- naive_means - mem_means
group_diff <- as.data.frame(group_diff)
rownames(group_diff) <- rownames(diff_df)
saveRDS(group_diff, file.path(sannot_dir, "group_diff_cd4tmem_cvar.rds"))

# Group difference for T_mem_CD8 and T_naive
group_diff <- naive_means - mem8_means
group_diff <- as.data.frame(group_diff)
rownames(group_diff) <- rownames(diff_df)
saveRDS(group_diff, file.path(sannot_dir, "group_diff_cd8tmem_cvar.rds"))

# Convert the pseudo_chrom summarized experiment object into chromvar object
dev <- assays(pseudo_chrom)$deviations
rownames(pseudo_chrom) <- rownames(z_scores) <- rownames(dev) <- NULL

# Convert the zmat to a SummarizedExperiment
se <- SummarizedExperiment(
  assays = list(z = z_scores, deviations = dev),
  colData = colData(pseudo_chrom), rowData = rowData(pseudo_chrom)
)

chromvar_devs <- as(se, "chromVARDeviations")
rownames(chromvar_devs) <- rowData(pseudo_chrom)$name

# Subset CD4 and naive cells
idxSample <- BiocGenerics::which(colData(chromvar_devs)$cell_types %in% c("T_mem_CD4", "T_naive"))
cellsSample <- chromvar_devs[, idxSample]

# Run differential analysis
diff_acc <- differentialDeviations(cellsSample, "cell_types", "two.sided")
diff_acc$isDiff <- ifelse(diff_acc$p_value_adjusted < 0.05, "Diff", "NotDiff")
saveRDS(diff_acc, file.path(sannot_dir, "diff_res_cd4tmem_cvar.rds"))

# Subset CD8 and naive cells
idxSample <- BiocGenerics::which(colData(chromvar_devs)$cell_types %in% c("T_mem_CD8", "T_naive"))
cellsSample <- chromvar_devs[, idxSample]

# Run differential analysis
diff_acc <- differentialDeviations(cellsSample, "cell_types", "two.sided")
diff_acc$isDiff <- ifelse(diff_acc$p_value_adjusted < 0.05, "Diff", "NotDiff")
saveRDS(diff_acc, file.path(sannot_dir, "diff_res_cd8tmem_cvar.rds"))

#####################################################################
# Plot z-difference combined with differential analysis for CD4+ cells
#####################################################################

# Load the data
group_diff_cd4 <- readRDS(file.path(sannot_dir, "group_diff_cd4tmem_cvar.rds"))
group_diff_thmem <- readRDS(file.path(sannot_dir, "group_diff_thmem_mtfr.rds"))
diff_res_cd4 <- readRDS(file.path(sannot_dir, "diff_res_cd4tmem_cvar.rds"))
diff_res_thmem <- readRDS(file.path(sannot_dir, "diff_res_TH_two.sided.rds"))
rownames(diff_res_thmem) <- diff_res_thmem$motifs

# Combine the data for CD4+ cells
combined_cd4 <- data.frame(
  cVAR_diff = group_diff_cd4$group_diff,
  mTFR_diff = group_diff_thmem$diff,
  chromVAR_padj = diff_res_cd4[rownames(group_diff_cd4), ]$p_value_adjusted,
  mTFR_padj = diff_res_thmem[rownames(group_diff_thmem), ]$p_value_adjusted,
  motifs = rownames(group_diff_cd4)
)


combined_cd4 <- combined_cd4 %>%
  mutate(color = case_when(
    abs(mTFR_diff) > 0.5 & mTFR_padj < 0.05 & abs(cVAR_diff) > 0.5 & chromVAR_padj < 0.05 ~ "BothDiff",
    abs(mTFR_diff) > 0.5 & mTFR_padj < 0.05 ~ "mTFRDiff",
    abs(cVAR_diff) > 0.5 & chromVAR_padj < 0.05 ~ "cVARDiff",
    TRUE ~ "NotDiff"
  ))

# Calculate the correlation between cVAR_diff and mTFR_diff
correlation_score <- cor(combined_cd4$cVAR_diff, combined_cd4$mTFR_diff, method = "pearson")
# -0.2813568

# Add a new column that contains the names of motifs to be annotated based on significance
combined_cd4 <- combined_cd4 %>%
  mutate(motif_label = if_else(
    (abs(mTFR_diff) > 0.5 & mTFR_padj < 0.01) |
      (abs(cVAR_diff) > 0.5 & chromVAR_padj < 0.01),
    true = motifs,
    false = NA_character_
  ))

# Plot CD4+ cells with dashed lines at x=0 and y=0, and add correlation score
p_cd4 <- ggplot(combined_cd4, aes(x = cVAR_diff, y = mTFR_diff, color = color)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c(
    "BothDiff" = "red", "mTFRDiff" = "blue",
    "cVARDiff" = "purple", "NotDiff" = "skyblue"
  )) +
  labs(
    x = "cVAR Z-score difference", y = "mTFR Z-score difference",
    title = "CD4 T Memory Cells (T_mem_CD4 vs Th-Mem)"
  ) +
  theme_classic(base_size = 16) + # Increase base text size
  annotate("text",
    x = 1.8, y = 1.8, # Adjusted to the top-right corner
    label = paste("Correlation: ", round(correlation_score, 2)),
    size = 6, hjust = 1
  ) + # Increase annotation text size
  xlim(-2, 2) + # Set x-axis limits
  ylim(-2, 2) + # Set y-axis limits
  geom_text(aes(label = motif_label),
    size = 5, vjust = 1.5, hjust = 0.5, check_overlap = TRUE
  ) # Increase motif label text size

# Save the plot as PDF
pdf_file <- file.path(plot_dir, "scatter_plot_cd4tmem_vs_thmem_diffs_with_motifs.pdf")
ggsave(pdf_file, plot = p_cd4, width = 8, height = 6)

#####################################################################
# Do the same for CD8+ cells
#####################################################################

# Load the data
group_diff_cd8 <- readRDS(file.path(sannot_dir, "group_diff_cd8tmem_cvar.rds"))
group_diff_tc <- readRDS(file.path(sannot_dir, "group_diff_tcmem_mtfr.rds"))
diff_res_cd8 <- readRDS(file.path(sannot_dir, "diff_res_cd8tmem_cvar.rds"))
diff_res_tc <- readRDS(file.path(sannot_dir, "diff_res_TC_two.sided.rds"))
rownames(diff_res_tc) <- diff_res_tc$motifs

# Combine the data for CD8+ cells
combined_cd8 <- data.frame(
  cVAR_diff = group_diff_cd8$group_diff,
  mTFR_diff = group_diff_tc$diff,
  chromVAR_padj = diff_res_cd8[rownames(group_diff_cd8), ]$p_value_adjusted,
  mTFR_padj = diff_res_tc[rownames(group_diff_tc), ]$p_value_adjusted,
  motifs = rownames(group_diff_cd8)
)

combined_cd8 <- combined_cd8 %>%
  mutate(color = case_when(
    abs(mTFR_diff) > 0.5 & mTFR_padj < 0.05 & abs(cVAR_diff) > 0.5 & chromVAR_padj < 0.05 ~ "BothDiff",
    abs(mTFR_diff) > 0.5 & mTFR_padj < 0.05 ~ "mTFRDiff",
    abs(cVAR_diff) > 0.5 & chromVAR_padj < 0.05 ~ "cVARDiff",
    TRUE ~ "NotDiff"
  ))

# Calculate the correlation between cVAR_diff and mTFR_diff
correlation_score <- cor(combined_cd8$cVAR_diff, combined_cd8$mTFR_diff, method = "pearson")
#-0.242243

# Add a new column that contains the names of motifs to be annotated based on significance
combined_cd8 <- combined_cd8 %>%
  mutate(motif_label = if_else(
    (abs(mTFR_diff) > 0.5 & mTFR_padj < 0.05) |
      (abs(cVAR_diff) > 0.5 & chromVAR_padj < 0.05),
    true = motifs,
    false = NA_character_
  ))

# Plot CD8+ cells with dashed lines at x=0 and y=0, and add correlation score
p_cd8 <- ggplot(combined_cd8, aes(x = cVAR_diff, y = mTFR_diff, color = color)) +
  geom_point(size = 2, alpha = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  scale_color_manual(values = c(
    "BothDiff" = "red", "mTFRDiff" = "blue",
    "cVARDiff" = "purple", "NotDiff" = "skyblue"
  )) +
  labs(
    x = "cVAR Z-score difference", y = "mTFR Z-score difference",
    title = "CD8 T Memory Cells (T_mem_CD8 vs Tc-Mem)"
  ) +
  theme_classic(base_size = 16) + # Increase base text size
  annotate("text",
    x = 1.8, y = 1.8, # Adjusted to the top-right corner
    label = paste("Correlation: ", round(correlation_score, 2)),
    size = 6, hjust = 1
  ) + # Increase annotation text size
  xlim(-2, 2) + # Set x-axis limits
  ylim(-2, 2) + # Set y-axis limits
  geom_text(aes(label = motif_label),
    size = 5, vjust = 1.5, hjust = 0.5, check_overlap = TRUE
  ) # Increase motif label text sizee

# Save the plot as PDF
pdf_file <- file.path(plot_dir, "scatter_plot_cd8tmem_vs_tcmem_diffs_with_motifs.pdf")
ggsave(pdf_file, plot = p_cd8, width = 8, height = 6)

#####################################################################
