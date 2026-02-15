
#!/usr/bin/env Rscript

#####################################################################
# one_vs_other.R
# Created on 23-04-2025 by Irem Gunduz
# Compare one cell type vs other cell types
#####################################################################


# Load required libraries

suppressPackageStartupMessages({
  library(dplyr)
  library(ggrepel)
  library(ggplot2)
  library(ArchR)
  library(muLogR)
  library(circlize)
  library(jmuOutlier)
  library(ComplexHeatmap)
  library(ggplot2)
  library(ggpubr) 
  library(VennDiagram)
  library(RColorBrewer)
})
set.seed(12)

fig_dir <- "/icbb/projects/igunduz/Figure_4_230425/"
if (!dir.exists(fig_dir)) {dir.create(fig_dir)}

# Load the project
outputDir <- "/icbb/projects/igunduz/archr_300824/icbb/projects/igunduz/archr_project_011023/"
plot_dir <- "/icbb/projects/igunduz/methylTFR_081124/Figure_4_heatmaps_051224/"
if (!dir.exists(plot_dir)) {dir.create(plot_dir)}
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)

cell_annot_atac <- read.delim("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/sample_annot_atac.tsv")
rownames(cell_annot_atac) <- cell_annot_atac$cellId_archr
annot_acc <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/cell_matching/annot_acc_nearest_meth_for_acc.rds")
annot_acc$archr_id <- rownames(annot_acc)
new_annot <- merge(cell_annot_atac, annot_acc, by = "row.names")
new_annot$atac_pb_id <- paste0(new_annot$mapped_cell_class, "_", new_annot$sample_sampleId_cminid)

# Subset the cells
idxSample <- BiocGenerics::which(project$cellNames %in% new_annot$archr_id)
cellsSample <- project$cellNames[idxSample]
project <- project[cellsSample, ]

# Add mapped cell-type information
project <- addCellColData(
  ArchRProj = project, data = new_annot$mapped_cell_class,
  name = "mapped_cell_class", cells = new_annot$archr_id
)

# Add mapped cell-type information
project <- addCellColData(
  ArchRProj = project, data = new_annot$atac_pb_id,
  name = "atac_pb_id", cells = new_annot$archr_id
)

# Remove Other-cells from the analysis
idxSample <- BiocGenerics::which(project$mapped_cell_class != "Other-cell")
cellsSample <- project$cellNames[idxSample]
project <- project[cellsSample, ]



#####################################################################
# One vs Other wilcox test for chromVAR
#####################################################################
cell_order <- c("Monocyte", "NK-cell", "B-cell", "Th-Naive", "Tc-Naive", "Th-Mem", "Tc-Mem")

# Get the chromVAR matrix
pseudo_chrom <- ArchR::getGroupSE(project, "jaspar2020Matrix", groupBy = "atac_pb_id", divideN = TRUE)
seZ <- pseudo_chrom[rowData(pseudo_chrom)$seqnames == "z", ]
zmat <- assay(seZ)
zmat <- as.matrix(zmat)
rownames(zmat) <- rowData(seZ)$name
#zmat <- zmat[row_order, ]
zmat <- methylTFR:::computeRowZScore(zmat)
saveRDS(zmat, paste0(fig_dir, "zmat_pb.rds"))

# Define the wilcoxon_helper function
wilcoxon_helper <- function(x, groups, alternative) {
  splitx <- split(x, groups)
  return(wilcox.test(splitx[[1]], splitx[[2]],
    alternative = alternative,
    paired = FALSE
  )$p.value)
}

# Initialize an empty list to store results
all_results <- list()

# Loop through each cell type
for (cell_type in cell_order) {
  # Get the cell type
  idx_cell <- grepl(cell_type, colnames(zmat))
  groups <- ifelse(idx_cell, "cell", "other")
  
  # Perform Wilcoxon test using wilcoxon_helper
  wilcox_test <- apply(zmat, 1, wilcoxon_helper, groups = groups, alternative = "two.sided")
  
  # Adjust p-values
  wilcox_test_adj <- p.adjust(wilcox_test, method = "fdr")
  
  # Combine results into a data.frame
  cell_results <- data.frame(
    TF = rownames(zmat),
    p_value = wilcox_test,
    p_value_adj = wilcox_test_adj,
    cell_type = cell_type#,isDiff = wilcox_test_adj < 0.05
  )
  
  # Append to the list
  all_results[[cell_type]] <- cell_results
}

# Combine all results into a single data.frame
final_results <- do.call(rbind, all_results)
rownames(final_results) <- NULL

# Save the combined results
saveRDS(final_results, file.path(fig_dir, "cVAR_combined_wilcox_results.rds"))


#####################################################################
# methylTFR one vs other
#####################################################################
cell_order <- c("Monocyte", "NK-cell", "B-cell", "Th-Naive", "Tc-Naive", "Th-Mem", "Tc-Mem")
mtfr_dir <- "/icbb/projects/igunduz/mTFR_bias_fix_v3/all_pseudobulks_310824/jaspar2020_distal/"
result <- list.files(mtfr_dir, pattern = "deviations.RDS", full.names = TRUE)

# Get methylTFR matrix
methylTFR <-  list()
for (x in seq_along(result)) {
    path <- result[x]
    cur_dev <- readRDS(path)
    cur_dev <- methylTFR::deviations(cur_dev) 
    methylTFR[[x]] <- as.matrix(cur_dev)
  }
methylTFR <- do.call(base::cbind, methylTFR)
methylTFR <- methylTFR:::computeRowZScore(methylTFR)
colnames(methylTFR) <- sub("\\.bedGraph\\.bed$", "", colnames(methylTFR))


# Initialize an empty list to store results
all_results_mtfr <- list()

# Loop through each cell type
for (cell_type in cell_order) {
  # Get the cell type
  idx_cell <- grepl(cell_type, colnames(methylTFR))
  groups <- ifelse(idx_cell, "cell", "other")
  
  # Perform Wilcoxon test using wilcoxon_helper
  wilcox_test <- apply(methylTFR, 1, wilcoxon_helper, groups = groups, alternative = "two.sided")
  
  # Adjust p-values
  wilcox_test_adj <- p.adjust(wilcox_test, method = "fdr")
  
  # Combine results into a data.frame
  cell_results <- data.frame(
    TF = rownames(methylTFR),
    p_value = wilcox_test,
    p_value_adj = wilcox_test_adj,
    cell_type = cell_type#,isDiff = wilcox_test_adj < 0.05
  )
  
  # Append to the list
  all_results_mtfr[[cell_type]] <- cell_results
}
# Combine all results into a single data.frame
final_results_mtfr <- do.call(rbind, all_results_mtfr)
rownames(final_results_mtfr) <- NULL

# Save the combined results
saveRDS(final_results_mtfr, file.path(fig_dir, "methylTFR_combined_wilcox_results.rds"))

#####################################################################
# Plot cell-specific mTFR and cVAR results with differentials
#####################################################################
cell_order <- c("Monocyte", "NK-cell", "B-cell", "Th-Naive", "Tc-Naive", "Th-Mem", "Tc-Mem")
final_results_mtfr <- readRDS(file.path(fig_dir, "methylTFR_combined_wilcox_results.rds"))
final_results <- readRDS(file.path(fig_dir, "cVAR_combined_wilcox_results.rds"))

# Get the chromVAR matrix
pseudo_chrom <- ArchR::getGroupSE(project, "jaspar2020Matrix", groupBy = "mapped_cell_class", divideN = TRUE)
seZ <- pseudo_chrom[rowData(pseudo_chrom)$seqnames == "z", ]
zmat <- assay(seZ)
zmat <- as.matrix(zmat)
rownames(zmat) <- rowData(seZ)$name
#zmat <- methylTFR:::computeRowZScore(zmat)

# Get the methylTFR matrix
mtfr_zmat <- readRDS("/icbb/projects/igunduz/mTFR_bias_fix_v3/mTFR_310824/jaspar2020_distal/deviations.RDS")
mtfr_zmat <- methylTFR:::deviations(mtfr_zmat)
mtfr_zmat <- methylTFR:::computeRowZScore(mtfr_zmat)
colnames(mtfr_zmat) <- sub("\\.bedGraph", "", colnames(mtfr_zmat))

plot_list <- list()  # Store all plots

# Loop through each cell type
for(cell_type in cell_order){
  
  # Subset the data for the current cell type
  idx <- grepl(cell_type, colnames(mtfr_zmat))
  idx_cvar <- grepl(cell_type, colnames(zmat))
  methylTFR_sub <- as.data.frame(mtfr_zmat[, idx])
  zmat_sub <- as.data.frame(zmat[, idx_cvar])
  colnames(zmat_sub) <- colnames(methylTFR_sub) <- "zscore"
  res_mtfr_sub <- final_results_mtfr[final_results_mtfr$cell_type == cell_type, ]
  res_cvar_sub <- final_results[final_results$cell_type == cell_type, ]

  # Merge the results
  merged_results <- merge(res_mtfr_sub, res_cvar_sub, by = "TF", suffixes = c("_mtfr", "_cvar"))
  
  # Merge cVAR and methylTFR values
  zscores <- merge(methylTFR_sub, zmat_sub, by = "row.names", suffixes = c("_mtfr", "_cvar"))
  colnames(zscores) <- c("TF", "zscore_mtfr", "zscore_cvar")

 # Merge with the results
  final_sub <- merge(merged_results, zscores, by = "TF", suffixes = c("_mtfr", "_cvar"))

# Add isDiff_mtfr column
final_sub$isDiff_mtfr <- with(final_sub, p_value_mtfr < 0.05 & abs(zscore_mtfr) > 1)

# Add isDiff_cvar column
final_sub$isDiff_cvar <- with(final_sub, p_value_cvar < 0.05 & abs(zscore_cvar) > 1)

# Add isDiff column
final_sub$isDiff <- with(final_sub, ifelse(
  isDiff_mtfr & isDiff_cvar, "diff_both",
  ifelse(isDiff_mtfr, "diff_mtfr",
    ifelse(isDiff_cvar, "diff_cvar", "non_sig")
  )
))

# Remove NAs
final_sub <- final_sub[!is.na(final_sub$isDiff), ]

# Filter the data for differentials only
#differentials <- final_sub[final_sub$isDiff != "non_sig", ]

# Calculate the correlation for differentials only
correlation <- cor(final_sub$zscore_mtfr, final_sub$zscore_cvar, method = "spearman")
# Calculate axis limits based on the maximum absolute values
x_limit <- max(abs(final_sub$zscore_mtfr), na.rm = TRUE)
y_limit <- max(abs(final_sub$zscore_cvar), na.rm = TRUE)

# Create the scatter plot with correlation score and dynamic axis limits
p <- ggplot(final_sub, aes(x = zscore_mtfr, y = zscore_cvar)) +
        geom_point(aes(color = isDiff), alpha = 0.5) +
        scale_color_manual(values = c(
          "non_sig" = "grey",
          "diff_both" = "darkred",
          "diff_cvar" = "darkgreen",
          "diff_mtfr" = "darkblue"
        )) +
        geom_text_repel(data = final_sub[final_sub$isDiff != "non_sig", ], aes(label = TF), size = 3) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") + # Add dashed line for y-axis
        geom_vline(xintercept = 0, linetype = "dashed", color = "black") + # Add dashed line for x-axis
        #geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "solid") + # Add correlation line
        annotate("text", x = Inf, y = Inf, label = paste0("Correlation: ", round(correlation, 2)), 
                 hjust = 1.1, vjust = 1.5, size = 5, color = "black") + # Add correlation score
        labs(x = "mTFR z-score", y = "chromVAR z-score") +
        coord_cartesian(xlim = c(-x_limit, x_limit), ylim = c(-y_limit, y_limit)) + # Set axis limits
        theme_classic() +
        theme(legend.position = "bottom")

# Save the plot
ggsave(filename = paste0(fig_dir, cell_type, "_mtfr_vs_cvar.png"), plot = p, width = 8, height = 6)
ggsave(filename = paste0(fig_dir, cell_type, "_mtfr_vs_cvar.pdf"), plot = p, width = 8, height = 6)

    # Store the plot in the list
    plot_list[[cell_type]] <- p
}

# Combine all plots into a single PDF using ggarrange
pdf(file.path(fig_dir, "mtfr_vs_cvar_combined.pdf"), width = 15, height = ceiling(length(plot_list) / 4) * 5)

ggarrange(plotlist = plot_list, 
          ncol = 4, 
          nrow = ceiling(length(plot_list) / 4), 
          labels = names(plot_list),  # Cell-type labels ONCE, on top
          font.label = list(size = 14, face = "bold"))

dev.off()

#####################################################################
# Heatmap of cVAR with FACS sorted labels
#####################################################################
cell_type_colors <- c(
  #"Other-cell" = "#CCCCCC",
  "B-cell" = "#AE017E",
  "Monocyte" = "#CC4C02",
  "NK-cell" = "#A65628",
  "Th-Mem" = "#41B6C4",
  "Tc-Mem" = "#4292C6",
  "Tc-Naive" = "#888FB5",
  "Th-Naive" = "#C7E9B4"
)
row_order <- readRDS("/icbb/projects/igunduz/Figure_4_270824/motifs.rds")

# Get the chromVAR matrix
pseudo_chrom <- ArchR::getGroupSE(project, "jaspar2020Matrix", groupBy = "mapped_cell_class", divideN = TRUE)
seZ <- pseudo_chrom[rowData(pseudo_chrom)$seqnames == "z", ]
zmat <- assay(seZ)
zmat <- as.matrix(zmat)
rownames(zmat) <- rowData(seZ)$name
zmat <- methylTFR:::computeRowZScore(zmat)
zmat <- zmat[row_order, ]


logger.start("Create a data frame for the samples' conditions")
ann <- data.frame(Cell = colnames(zmat))
rownames(ann) <- colnames(zmat)

# Define the desired order of cell types
new_order <- c("B-cell", "Monocyte", "NK-cell", "Th-Naive", "Tc-Naive", "Tc-Mem", "Th-Mem")


# Create the column annotation object
column_ha <- HeatmapAnnotation(
  df = ann,  # Pass the updated annotation data
  col = list(Cell = cell_type_colors)  # Provide the colors for each cell type
)

# Set the color scheme
colors.cv <- ChrAccR::getConfigElement("colorSchemesCont")
colors.cv <- colors.cv[[".default.div"]]
c <- grDevices::colorRampPalette(colors.cv)(nrow(ann))

# Order columns by new_order and cluster within each group
zmat <- zmat[, order(factor(colnames(zmat), levels = new_order))]
column_split_factor <- factor(colnames(zmat), levels = new_order)

# Create the heatmap object
heatmap_obj <- Heatmap(
  zmat,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  top_annotation = column_ha,
  column_split = column_split_factor,  # Split columns by cell types
  cluster_column_slices = FALSE,  # Disable clustering between groups
  col = c,
  heatmap_legend_param = list(
    title = paste0("chromVAR")
  ),
  show_column_names = FALSE
)

logger.completed()

# Plot the heatmap
logger.info("Plotting the heatmap for chromVAR")
pdf(paste0(fig_dir, "chromVAR_motifs.pdf"), width = 20, height = 20)
draw(heatmap_obj)
dev.off()

