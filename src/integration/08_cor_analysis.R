#!/usr/bin/env Rscript

#####################################################################
# 08_cor_analysis.R
# created on 17-04-2025 by Irem Gunduz
# Correlation analysis for the motifs
#####################################################################

cell_type_colors_atac <- c(
  "B_mem" = "#AE017E",
  "B_naive" = "#F768A1",
  "DC" = "#67000D",
  "Mono_CD14" = "#FE9929",
  "Mono_CD16" = "#CC4C02",
  "NK_CD16" = "#A65628",
  "Plasma" = "#A106BD",
  "T_mait" = "#41B6C4",
  "T_mem_CD4" = "#4292c6",
  "T_mem_CD8" = "#0074cc",
  "T_mix" = "#888FB5",
  "T_naive" = "#C7E9B4"
)


cell_type_colors <- c(
  # "Other-cell" = "#CCCCCC",
  "B-cell" = "#AE017E",
  "Monocyte" = "#CC4C02",
  "NK-cell" = "#A65628",
  "Th-Mem" = "#41B6C4",
  "Tc-Mem" = "#4292C6",
  "Tc-Naive" = "#888FB5",
  "Th-Naive" = "#C7E9B4"
)

get_groupname <- function(x) {
  return(unlist(strsplit(x, split = "_"))[1])
}
# load libraries
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(muLogR)
  library(GenomicRanges)
  library(methylTFR)
  library(ggplot2)
  library(circlize)
  library(ComplexHeatmap)
})

set.seed(12) # set seed
mtfr_dir <- "/icbb/projects/igunduz/methylTFR_081124/all_pseudobulks_121124/jaspar2020_distal/"

# Load the project
outputDir <- "/icbb/projects/igunduz/archr_300824/icbb/projects/igunduz/archr_project_011023/"
plot_dir <- "/icbb/projects/igunduz/methylTFR_081124/Figure_4_heatmaps_051224/"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)
annot_acc <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/cell_matching/annot_acc_cellType.rds")
row_order <- readRDS("/icbb/projects/igunduz/Figure_4_270824/motifs.rds")


cell_annot_atac <- read.delim("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/sample_annot_atac.tsv")
rownames(cell_annot_atac) <- cell_annot_atac$cellId_archr
annot_acc <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/cell_matching/annot_acc_cellType.rds")
annot_acc$archr_id <- rownames(annot_acc)
new_annot <- merge(cell_annot_atac, annot_acc, by = "row.names")
new_annot$atac_pb_id <- paste0(new_annot$class, "_", new_annot$sample_sampleId_cminid)

# Subset the cells
idxSample <- BiocGenerics::which(project$cellNames %in% new_annot$archr_id)
cellsSample <- project$cellNames[idxSample]
project <- project[cellsSample, ]

# Add mapped cell-type information
project <- addCellColData(
  ArchRProj = project, data = new_annot$atac_pb_id,
  name = "atac_pb_id", cells = new_annot$archr_id
)

# Get the chromVAR matrix
pseudo_chrom <- ArchR::getGroupSE(project, "jaspar2020Matrix", groupBy = "atac_pb_id", divideN = TRUE)
seZ <- pseudo_chrom[rowData(pseudo_chrom)$seqnames == "z", ]
zmat <- assay(seZ)
zmat <- as.matrix(zmat)
rownames(zmat) <- rowData(seZ)$name
# zmat <- zmat[row_order, ]
saveRDS(zmat, "/icbb/projects/igunduz/Figure_4_270824/cVAR_zscore.rds")

zmat <- methylTFR:::computeRowZScore(zmat)

zmat[round(zmat, 2) > 5] <- 5
zmat[round(zmat, 2) < -5] <- -5

zmat <- zmat[row_order, ]

# Filter the motifs
# zmat[abs(round(zmat, 2)) < 0.01] <- NA
# zmat[abs(round(zmat, 2)) > 5] <- NA
# zmat <- zmat[complete.cases(zmat), ]
# dim(zmat)
# max(zmat)
# min(zmat)


# Create the heatmap for chromVAR
logger.start("Create a heatmap for atac...")

# Order the cells
# Order the cells
new_order <- c("B_naive", "B_mem", "Plasma", "Mono_CD14", "Mono_CD16", "DC", "NK_CD16", "T_mait", "T_mem_CD4", "T_mem_CD8", "T_mix", "T_naive")

# Find the positions of the new_order items in the zmat column names
ordered_indices <- sapply(new_order, function(x) grep(paste0("^", x), colnames(zmat)))
ordered_indices <- unlist(ordered_indices)
zmat <- zmat[, ordered_indices]

# Extract the base cell types from colnames(zmat) by matching them to new_order
base_cell_types <- sapply(colnames(zmat), function(x) {
  match <- new_order[sapply(new_order, function(y) grepl(paste0("^", y), x))]
  if (length(match) > 0) {
    return(match[1])
  } else {
    return(NA)
  }
})

# Create the annotation data table
ann <- data.table(Cell = base_cell_types)
rownames(ann) <- colnames(zmat)
ann$Cell <- factor(ann$Cell, levels = new_order) # Match levels to new_order

# Create the column annotation object
column_ha <- HeatmapAnnotation(
  df = ann, # Pass the updated annotation data
  col = list(Cell = cell_type_colors_atac) # Provide the colors for each cell type
)

# Set the color scheme
colors.cv <- ChrAccR::getConfigElement("colorSchemesCont")
colors.cv <- colors.cv[[".default.div"]]
c <- grDevices::colorRampPalette(colors.cv)(nrow(ann))

# Order columns by new_order and cluster within each group
zmat <- zmat[, order(factor(base_cell_types, levels = new_order))]
column_split_factor <- factor(base_cell_types, levels = new_order)

# Create the heatmap object
heatmap_obj <- Heatmap(
  zmat,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  top_annotation = column_ha,
  column_split = column_split_factor, # Split columns by cell types
  cluster_column_slices = FALSE, # Disable clustering between groups
  col = c,
  heatmap_legend_param = list(
    title = paste0("chromVAR")
  ),
  show_column_names = FALSE
)

logger.completed()

# Plot the heatmap
logger.info("Plotting the heatmap for chromVAR")
pdf(paste0(plot_dir, "chromVAR_motifs.pdf"), width = 20, height = 20)
draw(heatmap_obj)
dev.off()

row_order <- rownames(zmat[row_order(heatmap_obj), ])
saveRDS(row_order, "/icbb/projects/igunduz/Figure_4_270824/motifs.rds")


#######################################################################
# methylTFR Heatmap
#######################################################################

mtfr_dir <- "/icbb/projects/igunduz/mTFR_bias_fix_v3/all_pseudobulks_310824/jaspar2020_distal/"
result <- list.files(mtfr_dir, pattern = "deviations.RDS", full.names = TRUE)

# Get methylTFR matrix
methylTFR <- list()
for (x in seq_along(result)) {
  path <- result[x]
  cur_dev <- readRDS(path)
  cur_dev <- methylTFR::deviations(cur_dev)
  methylTFR[[x]] <- as.matrix(cur_dev)
}
methylTFR <- do.call(base::cbind, methylTFR)
methylTFR <- methylTFR:::computeRowZScore(methylTFR)

methylTFR[round(methylTFR, 2) > 5] <- 5
methylTFR[round(methylTFR, 2) < -5] <- -5

# Remove ".bedGraph.bed" from column names in methylTFR
colnames(methylTFR) <- sub("\\.bedGraph\\.bed$", "", colnames(methylTFR))
saveRDS(methylTFR, "/icbb/projects/igunduz/Figure_4_270824/mTFR_zscore.rds")

mtfr_devs <- methylTFR[row_order, ]

# Extract cell types from column names
get_celltype <- function(x) strsplit(x, "_")[[1]][1]
groups <- sapply(colnames(mtfr_devs), get_celltype)

logger.start("Create a data frame for the samples' conditions")
ann <- data.frame(Cell = groups)
rownames(ann) <- colnames(mtfr_devs)

# Define the desired order of cell types
new_order <- c("B-cell", "Monocyte", "NK-cell", "Th-Naive", "Tc-Naive", "Tc-Mem", "Th-Mem")

# Ensure the annotation matches the new order
ann$Cell <- factor(ann$Cell, levels = new_order)

# Reorder columns of mtfr_devs and annotation data frame to match new_order
ordered_indices <- order(factor(groups, levels = new_order)) # Order based on factor levels
mtfr_devs <- mtfr_devs[, ordered_indices]
ann <- ann[ordered_indices, , drop = FALSE]

logger.info("Creating the heatmap object")
column_ha <- HeatmapAnnotation(df = ann, col = list(Cell = cell_type_colors))

# Define a color palette
c <- muRtools::colpal.cont(100, "cptcity.arendal_temperature")

# Create a factor for column splitting based on the new order
column_split_factor <- factor(ann$Cell, levels = new_order)

# Create the heatmap object with clustering within each cell type
heatmap_obj <- Heatmap(as.matrix(mtfr_devs),
  cluster_rows = FALSE, # Disable row clustering
  cluster_columns = TRUE, # Enable column clustering
  top_annotation = column_ha,
  col = c,
  column_split = column_split_factor, # Split columns by the desired order
  cluster_column_slices = FALSE, # Prevent clustering of the groups themselves
  heatmap_legend_param = list(
    title = paste0("methylTFR")
  ),
  show_column_names = FALSE
)

# Save the heatmap to a PDF
pdf(paste0(plot_dir, "methylTFR_motifs2.pdf"), width = 20, height = 20)
draw(heatmap_obj)
dev.off()

#######################################################################
# Add SELEX annotation
#######################################################################
selex <- read.csv("/icbb/projects/igunduz/exposure_atlas_manuscript/sample_annots/Selex_data.csv", skip = 20, header = 21, sep = ";")[, c(1, 2, 3, 4, 6)]
motifs <- readRDS("/icbb/projects/igunduz/Figure_4_270824/motifs.rds")

table(motifs %in% selex$TF.name)
#> table(row_order %in% selex$TF.name)

# FALSE  TRUE
#   47    15


# Add a new annotation column based on selex$Call
# Ensure the order of `motifs` matches the heatmap rows
selex_annotation <- data.frame(
  Call = ifelse(motifs %in% selex$TF.name,
    selex$Call[match(motifs, selex$TF.name)],
    "Inconclusive"
  ) # Replace NA with "Inconclusive"
)

# Replace NA values with "Inconclusive"
rownames(selex_annotation) <- motifs

# Create a heatmap annotation for the new column
row_ha <- rowAnnotation(
  Call = selex_annotation$Call,
  col = list(Call = c(
    "MethylPlus" = "blue",
    "MethylMinus" = "red",
    "Little effect and MethylPlus" = "green",
    "Inconclusive" = "gray"
  ))
)


selex_annotation$Call[is.na(selex_annotation$Call)] <- "Inconclusive"

# Create a single-column heatmap annotation
row_ha <- Heatmap(selex_annotation$Call,
  name = "Call",
  col = c(
    "MethylPlus" = "blue",
    "MethylMinus" = "red",
    "Little effect and MethylPlus" = "green",
    "Inconclusive" = "gray"
  ),
  cluster_rows = FALSE,
  show_row_names = TRUE,
  width = unit(1, "cm")
)

# Save the heatmap to a PDF
pdf(paste0(plot_dir, "selex_annotation_heatmap.pdf"), width = 5, height = 20)
draw(row_ha)
dev.off()



#######################################################################
# Correlation heatmap (based on FACS label - Pseudobulk level)
#######################################################################
new_annot$atac_pb_id_mapped <- paste0(new_annot$mapped_cell_class, "_", new_annot$sample_sampleId_cminid)
# Add mapped cell-type information
project <- addCellColData(
  ArchRProj = project, data = new_annot$atac_pb_id_mapped,
  name = "atac_pb_id_mapped", cells = new_annot$archr_id,
  force = TRUE
)
row_order <- readRDS("/icbb/projects/igunduz/Figure_4_270824/motifs.rds")

# Get the chromVAR matrix
pseudo_chrom <- ArchR::getGroupSE(project, "jaspar2020Matrix", groupBy = "atac_pb_id_mapped", divideN = TRUE)
seZ <- pseudo_chrom[rowData(pseudo_chrom)$seqnames == "z", ]
zmat <- assay(seZ)
zmat <- as.matrix(zmat)
rownames(zmat) <- rowData(seZ)$name
# zmat <- zmat[row_order, ]
# saveRDS(zmat, "/icbb/projects/igunduz/Figure_4_270824/cVAR_zscore.rds")

zmat <- methylTFR:::computeRowZScore(zmat)
zmat[round(zmat, 2) > 5] <- 5
zmat[round(zmat, 2) < -5] <- -5
zmat <- zmat[row_order, ]

# Get methylTFR matrix
methylTFR <- list()
for (x in seq_along(result)) {
  path <- result[x]
  cur_dev <- readRDS(path)
  cur_dev <- methylTFR::deviations(cur_dev)
  methylTFR[[x]] <- as.matrix(cur_dev)
}
methylTFR <- do.call(base::cbind, methylTFR)
methylTFR <- methylTFR:::computeRowZScore(methylTFR)
colnames(methylTFR) <- sub("\\.bedGraph\\.bed$", "", colnames(methylTFR))

methylTFR[round(methylTFR, 2) > 5] <- 5
methylTFR[round(methylTFR, 2) < -5] <- -5

# Sort to match the order of methylTFR
zmat <- zmat[, colnames(methylTFR)]
mtfr_devs <- methylTFR[row_order, ]

# Calculate row-wise correlation
row_correlation <- sapply(seq_len(nrow(zmat)), function(i) {
  cor(zmat[i, ], mtfr_devs[i, ]) # ), use = "complete.obs")  # Only use rows with complete observations
})

# Convert to a data frame for better readability
row_correlation_df <- data.frame(
  RowName = rownames(zmat),
  Correlation = row_correlation
)

# Convert row-wise correlation vector to a matrix for heatmap plotting
correlation_matrix <- as.matrix(row_correlation)
rownames(correlation_matrix) <- rownames(zmat) # Add row names for better interpretation

# Create the heatmap
cm <- Heatmap(
  correlation_matrix,
  name = "Correlation", # Name for the heatmap legend
  cluster_rows = FALSE, # Cluster rows
  show_row_names = TRUE, # Show row names
  cluster_columns = FALSE, # No clustering for columns since it's a single column
  # col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),  # Color scale
  heatmap_legend_param = list(title = "Row Correlation") # Legend settings
)

pdf(paste0(plot_dir, "correlation_heatmap.pdf"), width = 20, height = 20)
draw(cm)
dev.off()

#######################################################################
# Correlation plots for SELEX annotation
#######################################################################
# SELEX data
selex <- read.csv("/icbb/projects/igunduz/exposure_atlas_manuscript/sample_annots/Selex_data.csv", skip = 20, header = 21, sep = ";")[, c(1, 2, 3, 4, 6)]
selex <- selex[, c("TF.name", "methyl.SELEX.call")]

# cVAR matrix for all motifs
seZ <- pseudo_chrom[rowData(pseudo_chrom)$seqnames == "z", ]
zmat <- assay(seZ)
zmat <- as.matrix(zmat)
rownames(zmat) <- rowData(seZ)$name
zmat <- methylTFR:::computeRowZScore(zmat)
zmat[round(zmat, 2) > 5] <- 5
zmat[round(zmat, 2) < -5] <- -5

# mTFR matrix for all motifs
methylTFR <- list()
for (x in seq_along(result)) {
  path <- result[x]
  cur_dev <- readRDS(path)
  cur_dev <- methylTFR::deviations(cur_dev)
  methylTFR[[x]] <- as.matrix(cur_dev)
}
methylTFR <- do.call(base::cbind, methylTFR)
methylTFR <- methylTFR:::computeRowZScore(methylTFR)
colnames(methylTFR) <- sub("\\.bedGraph\\.bed$", "", colnames(methylTFR))
methylTFR[round(methylTFR, 2) > 5] <- 5
methylTFR[round(methylTFR, 2) < -5] <- -5

# Sort to match the order of methylTFR
zmat <- zmat[, colnames(methylTFR)]
methylTFR <- methylTFR[rownames(zmat), ]

# Calculate row-wise correlation
row_correlation <- sapply(seq_len(nrow(zmat)), function(i) {
  cor(zmat[i, ], methylTFR[i, ]) # ), use = "complete.obs")  # Only use rows with complete observations
})

# Convert to a data frame for better readability
row_correlation_df <- data.frame(
  TF.name = rownames(zmat),
  Correlation = row_correlation
)

# Add SELEX annotation and filter methylPlus and methylMinus
row_correlation_df <- merge(row_correlation_df, selex, by = "TF.name")
row_correlation_df <- row_correlation_df[row_correlation_df$methyl.SELEX.call %in% c("MethylPlus", "MethylMinus"), ]

# Ensure row_correlation_df is a proper data frame
row_correlation_df <- data.frame(row_correlation_df)

# Count the number of motifs in each class using base R
motif_counts <- as.data.frame(table(row_correlation_df$methyl.SELEX.call))
colnames(motif_counts) <- c("methyl.SELEX.call", "count")


# Create the boxplot with jittered points and annotations
p <- ggplot(row_correlation_df, aes(x = methyl.SELEX.call, y = Correlation, fill = methyl.SELEX.call)) +
  geom_boxplot(color = "black", outlier.shape = NA, width = 0.6) + # Boxplot without outliers
  geom_jitter(aes(color = methyl.SELEX.call), width = 0.2, size = 1.5, alpha = 0.7) + # Add jittered points
  scale_fill_manual(values = c("MethylPlus" = "darkgreen", "MethylMinus" = "darkred")) + # Custom colors for boxplot
  scale_color_manual(values = c("MethylPlus" = "darkgreen", "MethylMinus" = "darkred")) + # Custom colors for points
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) + # Set y-axis limits and breaks
  # Add whisker caps
  stat_boxplot(geom = "errorbar", width = 0.3, size = 0.8) + # Horizontal caps at whisker ends
  labs(
    x = "SELEX Group",
    y = "Correlation"
  ) +
  theme_classic(base_size = 14) + # Use a larger base font size for readability
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Center and bold the title
    axis.title = element_text(face = "bold"), # Bold axis titles
    legend.position = "none" # Remove legend if unnecessary
  ) +
  # Annotate the number of motifs in each class
  geom_text(
    data = motif_counts,
    aes(x = methyl.SELEX.call, y = 1, label = count), # Position annotations above the boxplots
    vjust = -0.5,
    size = 4
  )

# Save the plot
ggsave(
  filename = paste0(plot_dir, "correlation_boxplot_with_caps.pdf"),
  plot = p,
  width = 8,
  height = 6
)

#######################################################################
