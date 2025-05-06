#!/usr/bin/env Rscript

##########################################################################################################################################
# sub_hm_090425.R
# Created by Irem Gunduz on 27-04-2025
# This script generates heatmaps for the methylation data per monocyte subtypes
##########################################################################################################################################

suppressPackageStartupMessages({
  library(methylTFR)
  library(circlize)
  library(ComplexHeatmap)
})
set.seed(42)
mtfr_dir <- "/icbb/projects/igunduz/Figure_5_040425/covid_mtfr_070425/jaspar2020_distal/"
motifset <- "jaspar2020"
id <- "C19_mono2_vs_sev_090425"
condition <- c("C19_ctrl", "C19_sev")
plot_dir <- "/icbb/projects/igunduz/exposure_atlas_manuscript/src/meth_processing/Fig_5_sub_080425/plots/"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}
ds_dir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/Mono_CD14/data/dsATAC_filtered/"

##################################################################################
# Plot PCA
##################################################################################
# Read methylTFR deviations
mtfr_devs <- list.files(mtfr_dir, pattern = "Monocyte", full.names = TRUE)
mtfr_devs <- rlist::list.cbind(lapply(condition, function(x) {
  methylTFR::deviations(readRDS(mtfr_devs[grepl(x, mtfr_devs)]))
}))
mtfr_devs <- as.data.frame(mtfr_devs)

# Remove outliers
mtfr_devs <- mtfr_devs[, !colnames(mtfr_devs) %in% outliers]
mtfr_devs_z <- methylTFR:::computeRowZScore(as.matrix(mtfr_devs))

# Get the groups for mtfr_devs_z
groups <- unlist(lapply(colnames(mtfr_devs_z), get_group))

# Remove .bedGraph.bed suffix from colnames
colnames(mtfr_devs_z) <- gsub(".bedGraph.bed", "", colnames(mtfr_devs_z))
# Remove Monocyte_ prefix from colnames
colnames(mtfr_devs_z) <- gsub("Monocyte_", "", colnames(mtfr_devs_z))


# Apply PCA
tdf <- as.data.frame(t(mtfr_devs_z))
pca_result <- prcomp(tdf, center = TRUE, scale. = TRUE)
tdf$groups <- groups
sample_names <- rownames(tdf)

fig_path <- paste0(plot_dir, "PCA_mTFR_allmotifs_no_outliers.pdf")
pdf(fig_path, width = 10, height = 10)
# Plot PCA
autoplot(pca_result,
  data = tdf,
  colour = "groups",
  main = "PCA",
  size = 5
) +
  geom_text_repel(aes(label = sample_names), size = 3) + # Add readable labels
  theme_classic() +
  # scale_color_manual(values = cell_type_colors) +
  theme(legend.position = "bottom")
dev.off()

# Keep only C19_sev_mono1 and C19_ctrl
mtfr_devs_z <- mtfr_devs_z[, groups %in% c("C19_ctrl", "C19_sev_mono1")]
groups <- groups[groups %in% c("C19_ctrl", "C19_sev_mono1")]

# Apply PCA
tdf <- as.data.frame(t(mtfr_devs_z))
pca_result <- prcomp(tdf, center = TRUE, scale. = TRUE)
tdf$groups <- groups
sample_names <- rownames(tdf)
fig_path <- paste0(plot_dir, "PCA_mTFR_allmotifs_no_outliers_mono1.pdf")
pdf(fig_path, width = 10, height = 10)
# Plot PCA
autoplot(pca_result,
  data = tdf,
  colour = "groups",
  main = "PCA",
  size = 5
) +
  geom_text_repel(aes(label = sample_names), size = 3) + # Add readable labels
  theme_classic() +
  # scale_color_manual(values = cell_type_colors) +
  theme(legend.position = "bottom")
dev.off()
##################################################################################
# Plot heatmap for all group
##################################################################################

outliers <- c(
  "Monocyte_CoV_S_S8_D1.bedGraph.bed",
  "Monocyte_CoV_S_S15_D7.bedGraph.bed",
  "Monocyte_Ctrl_1_F_White_45yo.bedGraph.bed",
  "Monocyte_Ctrl_12_F_White_32yo.bedGraph.bed",
  "Monocyte_Ctrl_11_M_White_57yo.bedGraph.bed"
)[5]

mono_1 <- c(
  "Monocyte_CoV_S_S15_D1.bedGraph.bed",
  "Monocyte_CoV_S_S11_D3.bedGraph.bed",
  "Monocyte_CoV_S_S7_D1.bedGraph.bed",
  "Monocyte_CoV_S_S11_D1.bedGraph.bed"
)

# combined <- c(outliers, mono_1)

# Read methylTFR deviations
mtfr_devs <- list.files(mtfr_dir, pattern = "Monocyte", full.names = TRUE)
mtfr_devs <- rlist::list.cbind(lapply(condition, function(x) {
  methylTFR::deviations(readRDS(mtfr_devs[grepl(x, mtfr_devs)]))
}))
mtfr_devs <- as.data.frame(mtfr_devs)

# Remove outliers
mtfr_devs <- mtfr_devs[, !colnames(mtfr_devs) %in% outliers]
mtfr_devs_z <- methylTFR:::computeRowZScore(as.matrix(mtfr_devs))
motifs <- readRDS(paste0(ds_dir, "cvar_motifs_diff_080425.rds"))
mtfr_devs_z <- mtfr_devs_z[motifs, ]

# Remove .bedGraph.bed suffix from colnames
# colnames(mtfr_devs_z) <- gsub(".bedGraph.bed", "", colnames(mtfr_devs_z))

# Remove Monocyte_ prefix from colnames
# colnames(mtfr_devs_z) <- gsub("Monocyte_", "", colnames(mtfr_devs_z))

# Define the grouping logic
get_group <- function(filename) {
  if (filename %in% mono_1) {
    return("C19_sev_mono1")
  } else if (grepl("CoV", filename)) {
    return("C19_sev_mono2")
  } else if (grepl("Ctrl", filename)) {
    return("C19_ctrl")
  } else {
    return("Unknown")
  }
}

# Get the groups for mtfr_devs_z
groups <- unlist(lapply(colnames(mtfr_devs_z), get_group))

# Subset the data to include only the desired groups
valid_groups <- c("C19_ctrl", "C19_sev_mono1", "C19_sev_mono2")
mtfr_devs_z <- mtfr_devs_z[, groups %in% valid_groups]
groups <- groups[groups %in% valid_groups]

# Create the annotation data frame
ann <- data.frame(Exposure = groups)
rownames(ann) <- colnames(mtfr_devs_z)

# Define the group colors
group_colors <- c("C19_ctrl" = "#0072B2", "C19_sev_mono1" = "#D55E00", "C19_sev_mono2" = "#009E73")

# Create the heatmap annotation
column_ha <- HeatmapAnnotation(
  df = ann,
  col = list(Exposure = group_colors),
  annotation_legend_param = list(title = "Exposure")
)
c <- muRtools::colpal.cont(nrow(ann), "cptcity.arendal_temperature") # Daha fazla renk değeri

# Generate the heatmap
heatmap_obj <- Heatmap(
  mtfr_devs_z,
  cluster_rows = TRUE,
  cluster_columns = TRUE, # Enable clustering of columns
  column_split = ann$Exposure, # Split columns by Exposure groups
  top_annotation = column_ha, # Add the annotation to the top
  col = c, # Color palette
  heatmap_legend_param = list(
    title = "methylTFR"
  ),
  show_column_names = FALSE
)

# Save the heatmap to a PDF
pdf(paste0(plot_dir, "mTFR_allmotifs_no_outliers_mono2_CLUSTERED.pdf"), width = 20, height = 20)
draw(heatmap_obj)
dev.off()

##################################################################################
# Plot only C19_sev_mono1 and control
##################################################################################

# Read methylTFR deviations
mtfr_devs <- list.files(mtfr_dir, pattern = "Monocyte", full.names = TRUE)
mtfr_devs <- rlist::list.cbind(lapply(condition, function(x) {
  methylTFR::deviations(readRDS(mtfr_devs[grepl(x, mtfr_devs)]))
}))
mtfr_devs <- as.data.frame(mtfr_devs)

outliers <- c(
  "Monocyte_Ctrl_1_F_White_45yo.bedGraph.bed",
  "Monocyte_CCtrl_12_F_White_32yo.bedGraph.bed",
  "Monocyte_Ctrl_11_M_White_57yo.bedGraph.bed"
)
# Remove outliers
mtfr_devs <- mtfr_devs[, !colnames(mtfr_devs) %in% outliers]
mtfr_devs <- mtfr_devs[, groups %in% c("C19_ctrl", "C19_sev_mono1")]
mtfr_devs_z <- methylTFR:::computeRowZScore(as.matrix(mtfr_devs))
mtfr_devs_z <- mtfr_devs_z[motifs, ]

# Get the groups for mtfr_devs_z
groups <- unlist(lapply(colnames(mtfr_devs_z), get_group))

# Subset the data to include only the desired groups
valid_groups <- c("C19_ctrl", "C19_sev_mono1")
mtfr_devs_z <- mtfr_devs_z[, groups %in% valid_groups]
groups <- groups[groups %in% valid_groups]

# Remove .bedGraph.bed suffix from colnames
colnames(mtfr_devs_z) <- gsub(".bedGraph.bed", "", colnames(mtfr_devs_z))

# Remove Monocyte_ prefix from colnames
colnames(mtfr_devs_z) <- gsub("Monocyte_", "", colnames(mtfr_devs_z))

# Create the annotation data frame
ann <- data.frame(Exposure = groups)
rownames(ann) <- colnames(mtfr_devs_z)

# Define the group colors
group_colors <- c("C19_ctrl" = "#0072B2", "C19_sev_mono1" = "#D55E00")
# Create the heatmap annotation
column_ha <- HeatmapAnnotation(
  df = ann,
  col = list(Exposure = group_colors),
  annotation_legend_param = list(title = "Exposure")
)
c <- muRtools::colpal.cont(nrow(ann), "cptcity.arendal_temperature") # Daha fazla renk değeri
# Generate the heatmap
heatmap_obj <- Heatmap(
  mtfr_devs_z,
  cluster_rows = TRUE,
  cluster_columns = TRUE, # Enable clustering of columns
  column_split = ann$Exposure, # Split columns by Exposure groups
  top_annotation = column_ha, # Add the annotation to the top
  col = c, # Color palette
  heatmap_legend_param = list(
    title = "methylTFR"
  ),
  show_column_names = TRUE
)

# Save the heatmap to a PDF
pdf(paste0(plot_dir, "mTFR_allmotifs_no_outliers_mono1_CLUSTERED.pdf"), width = 20, height = 20)
draw(heatmap_obj)
dev.off()

# Get the row order
row_order <- row_order(heatmap_obj)
motifs <- rownames(mtfr_devs_z)[row_order]
saveRDS(motifs, paste0(plot_dir, "motifs_ordered.rds"))


##################################################################################
# Plot cVAR heatmap with mTFR order (Row Z-score based on batch corrected)
##################################################################################

# Plot heatmap
logger.start("Plotting heatmap for chromVAR deviations")
chromvar_mat <- readRDS(paste0(ds_dir, "chromvar_jaspar2020_101224_corrected_zscores.rds"))

# chromvar_mat <- readRDS(paste0(ds_dir, "chromvar_jaspar2020_080425_corrected.rds"))
rownames(chromvar_mat) <- sub(".*_", "", rownames(chromvar_mat))
# zmat <- methylTFR:::computeRowZScore(chromvar_mat)
zmat <- chromvar_mat

# Arrange max and min as 5 and -5
zmat[zmat > 5] <- 5
zmat[zmat < -5] <- -5

ann <- data.table(Exposure = ifelse(grepl("ctrl", colnames(zmat)), "C19_ctrl", "C19_sev"))
rownames(ann) <- colnames(zmat)
zmat <- zmat[motifs, ]


# Define the color scheme for the annotations
column_ha <- HeatmapAnnotation(
  Exposure = ann$Exposure, # Use Exposure groups for the annotation
  col = list(Exposure = c("C19_ctrl" = "blue", "C19_sev" = "red")), # Colors for groups
  annotation_legend_param = list(title = "Exposure") # Title for annotation legend
)

# Set the color scheme for the heatmap
colors.cv <- ChrAccR::getConfigElement("colorSchemesCont")
colors.cv <- colors.cv[[".default.div"]]
c <- grDevices::colorRampPalette(colors.cv)(nrow(ann))

# Create the heatmap object
heatmap_obj <- Heatmap(zmat,
  cluster_rows = FALSE,
  cluster_columns = TRUE, # Enable clustering of columns
  column_split = ann$Exposure, # Split columns by Exposure groups
  top_annotation = column_ha, # Add the annotation to the top
  col = c,
  heatmap_legend_param = list(
    title = "chromVAR"
  ),
  show_column_names = FALSE
)

logger.completed()

logger.info("Plotting the heatmap for chromVAR differentials")
pdf(paste0(plot_dir, "chromVAR_motifs_080425.pdf"), width = 20, height = 20)
draw(heatmap_obj)
dev.off()
##################################################################################
