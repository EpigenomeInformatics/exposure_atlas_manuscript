#!/usr/bin/env Rscript

#####################################################################
# 08_1_chraccr_plots.R
# created on 2023-10-02 by Irem Gunduz
# ChrAccR plots for COVID samples
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ChrAccR)
  library(ggplot2)
  library(data.table)
  library(stringr)
  library(muRtools)
  library(ggplot2)
  library(circlize)
  library(ComplexHeatmap)
  library(ggrepel)
  library(GenomicRanges)
  library(muLogR)
  library(clue)
})
set.seed(12)
source("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/utils/lola.R")
source("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/utils/chraccr_plots.R")
source("/icbb/projects/igunduz/sc_epigenome_exp/utils/helpers.R")

color_mapping <- c(
  "C19_ctrl" = "#6CA6CD",
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

diffCompNames <- c(
  "HIV_ctrl vs HIV_chr",
  "HIV_ctrl vs HIV_acu",
  "Influenza_d3 vs Influenza_ctrl",
  "Influenza_d30 vs Influenza_ctrl",
  "Influenza_d6 vs Influenza_ctrl",
  "OP_high vs OP_low",
  "OP_high vs OP_med",
  "OP_med vs OP_low"
)
cells <- c("B_mem", "B_naive", "Mono_CD14", "Mono_CD16", "NK_CD16", "T_mem_CD8", "T_mem_CD4", "T_naive", "T_mix", "T_naive", "T_mait")
outputDir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_2023-10-02/"
plotDir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/plots/"
if (!dir.exists(plotDir)) {
  dir.create(plotDir)
}
for (cell in cells) {
  # Iterate over each comparison in the comparison table
  comparisonTable <- readRDS(paste0(outputDir, cell, "/reports/differential_data/comparisonTable.rds"))
  comparisonTable$comps <- paste0(comparisonTable$grp1, " vs ", comparisonTable$grp2)
  not_present <- diffCompNames[which(!diffCompNames %in% comparisonTable$comps)]

  for (i in 1:nrow(comparisonTable)) {
    # Get the relevant information for the current comparison
    grp1 <- comparisonTable$grp1[i]
    grp2 <- comparisonTable$grp2[i]
    comp <- comparisonTable$comps[i]
    pp <- plotMAwithChrAccR(cell, outputDir, i, "archrPeaks")
    filename <- paste0(plotDir, "MA_plot_", cell, "_", grp1, " vs ", grp2, ".pdf")
    ggsave(filename = filename, plot = pp, width = 20, height = 20)
  }
}

outputDir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/"
diffCompNames <- c(
  "C19_mild vs C19_ctrl",
  "C19_mod vs C19_ctrl",
  "C19_sev vs C19_ctrl"
)
cells <- c("Mono_CD14", "NK_CD16", "T_mem_CD8", "T_mem_CD4", "T_naive")
plotDir <- "/icbb/projects/igunduz/finalize_echo_050824/MA_plots_251124/"
if (!dir.exists(plotDir)) {
  dir.create(plotDir)
}
for (cell in cells) {
  # Iterate over each comparison in the comparison table
  comparisonTable <- readRDS(paste0(outputDir, cell, "/reports/differential_data/comparisonTable.rds"))
  comparisonTable$comps <- paste0(comparisonTable$grp1, " vs ", comparisonTable$grp2)
  not_present <- diffCompNames[which(!diffCompNames %in% comparisonTable$comps)]

  for (i in 1:nrow(comparisonTable)) {
    # Get the relevant information for the current comparison
    grp1 <- comparisonTable$grp1[i]
    grp2 <- comparisonTable$grp2[i]
    comp <- comparisonTable$comps[i]
    pp <- plotMAwithChrAccR(cell, outputDir, i, "archrPeaks")
    filename <- paste0(plotDir, "MA_plot_", cell, "_", grp1, " vs ", grp2, ".pdf")
    ggsave(filename = filename, plot = pp, width = 20, height = 20)
  }
}

###############################################################################
# l2fc plots for COVID-MONO

# archr_peaks dataframe to the annotations
global_peaklist <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/regionsGR.rds") %>%
  data.table::as.data.table() %>%
  dplyr::select(!width & !strand) %>%
  as.data.frame()

outputDir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/"
cell <- "Mono_CD14"
plot_dir <- "/icbb/projects/igunduz/finalize_echo_050824/C19/"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}
cor_list <- list()
# ATAC get the differentially accessble regions
atac <- prepareDARforPlot(cell, outputDir, 3, "archrPeaks") # get the C19_severe
for (comp in diffCompNames[1:2]) {
  if (!file.exists(paste0(plot_dir, "L2FC_", cell, "_", comp, ".pdf"))) {
    j <- ifelse(comp == "C19_mild vs C19_ctrl", 1, 2)

    # Methylation, get the differntiall methylated regions
    atac2 <- prepareDARforPlot(cell, outputDir, j, "archrPeaks")
    # dmt <-  prepareDMRforPlot(cell,comp,path,global_peaklist,"sites")

    # filter methylation by overlapped regions of atac
    datatable <- prepareScatterMethAtac(atac, atac2, 4)
    cor_list[[comp]] <- cor.test(datatable$log2FC_1, datatable$log2FC_2)
    # plot the scatter plot
    scplot <- plotScatterL2FC(datatable, y_lab = paste0("Covid-19 ", j), x_lab = "Covid-19 Severe", comb = paste0("Covid-19 "), group1 = "log2FC_1", group2 = "log2FC_2")

    # save the plot
    ggsave(plot = scplot, paste0(plot_dir, "L2FC_", cell, "_", comp, ".pdf"), width = 20, height = 20)
    ggsave(plot = scplot, paste0(plot_dir, "L2FC_", cell, "_", comp, ".png"), width = 20, height = 20)
  }
}
saveRDS(cor_list, file = paste0(plot_dir, "cor_list.rds"))

###############################################################################

# load loladb database
lolaDb <- "/icbb/projects/share/annotations/lolaDB/hg38/"
lolaDb <- RnBeads::loadLolaDbs(lolaDb)
# lola_path <- "/icbb/projects/igunduz/DARPA_analysis/LOLA/"
lola_path <- "/icbb/projects/igunduz/finalize_echo_050824/C19/"
if (!dir.exists(lola_path)) {
  dir.create(lola_path)
}

# volcano plot for atac
p <- lolaVolcanoPlotC19("Mono_CD14", lolaDb, outputDir, pValCut = 1.5, region = "archrPeaks", database = "TF_motif_clusters", signifCol = "qValue", cnames = c("loss", "gain"), n = 3, colorpanel = c())$plot
ggsave(plot = p, paste0(lola_path, "Mono_atac_TF_motif_clusters.pdf"), width = 20, height = 20)

p <- lolaVolcanoPlotC19("Mono_CD14", lolaDb, outputDir, pValCut = 1.5, region = "archrPeaks", database = "TF_motifs", signifCol = "qValue", cnames = c("loss", "gain"), n = 3, colorpanel = c())$plot
ggsave(plot = p, paste0(lola_path, "Mono_atac_TF_motifs.pdf"), width = 20, height = 20)

p <- lolaVolcanoPlotC19("Mono_CD14", lolaDb, outputDir, pValCut = 1.5, region = "archrPeaks", database = "TF_motifs_vert", signifCol = "qValue", cnames = c("loss", "gain"), n = 3, colorpanel = c())$plot
ggsave(plot = p, paste0(lola_path, "Mono_atac_TF_motifs_vert.pdf"), width = 20, height = 20)

###############################################################################
# differential heatmap for COVID-MONO
logger.start("Starting differential heatmap for COVID-MONO")
diff <- prepareDARforPlot(cell, outputDir, 3, "archrPeaks") # , padj = 0.05) # get the C19_severe
# diff <- diff[isDiff == TRUE,] #filter the differentially accessible regions
diff$Start <- diff$Start + 1
diff$id <- paste0(diff$Chromosome, ":", diff$Start, "-", diff$End)
diffm <- diff
diff <- dplyr::select(diff, c("id", "isDiff"))

ds <- ChrAccR::loadDsAcc(paste0(outputDir, cell, "/data/dsATAC_processed/"))
logger.info("Getting fragment GR")
counts <- getCounts(ds, "archr_peaks")
# Get the counts matrix
coords <- data.table::as.data.table(getCoord(ds, "archr_peaks")) %>%
  dplyr::mutate(id = paste0(seqnames, ":", start, "-", end)) %>%
  dplyr::select(id)
# Apply k-means clustering
counts <- data.table::as.data.table(computeZScore(counts))

counts <- cbind(counts, coords)
counts <- merge(counts, diff)
counts <- counts[isDiff == TRUE, ]
counts <- dplyr::select(counts, !id & !isDiff)
assignments <- consensus_kmeans(counts, 5, 1)
diffm[diffm$isDiff == TRUE, "cluster"] <- assignments
diffm <- diffm[diffm$isDiff == TRUE, ]

# Reorder columns in the counts matrix and column annotations
counts_ordered <- setcolorder(counts, sort(names(counts)))

# Create column annotations based on the sample names
col_ann <- data.frame(Exposure = gsub("^([^_]*_[^_]*).*", "\\1", colnames(counts_ordered)))
rownames(col_ann) <- colnames(counts)
colnames(col_ann) <- "Exposure"

# Create row annotations based on cluster assignments
row_ann <- data.frame(Cluster = as.factor(assignments))
row_ann_ordered <- row_ann[order(row_ann$Cluster), ]

# counts_ordered$cluster_assignments <- assignments # add cluster info
setDT(counts_ordered)
rownames(counts_ordered) <- diffm$id

col_ha <- HeatmapAnnotation(
  Exposure = col_ann$Exposure,
  col = list(Exposure = c(
    "C19_ctrl" = "#4F609C",
    "C19_mild" = "#C43E96",
    "C19_mod" = "#06948E",
    "C19_sev" = "#C03830"
  ))
)

# Define the order of the columns based on their groups
group_order <- c("C19_ctrl", "C19_mild", "C19_mod", "C19_sev")
column_order <- colnames(counts_ordered)
column_order <- column_order[order(match(column_order, group_order))]

dend_cols <- cluster_within_group(counts_ordered, col_ann$Exposure)

# Create heatmap
ht <- Heatmap(
  counts_ordered,
  column_order = column_order,
  cluster_columns = FALSE,
  col = ArchR::paletteContinuous("solarExtra"),
  name = "Accessibility Z-score",
  cluster_rows = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  heatmap_legend_param = list(at = c(-max(round(counts_ordered)), 0, max(round(counts_ordered)))),
  column_title_side = "top",
  top_annotation = col_ha
)


pdf(paste0("/icbb/projects/igunduz/pheatmap_", cell, "_samplewise_C19_diff_all.pdf"), width = 20, height = 20)
draw(ht)
dev.off()


###############################################################################
# Add label to peaks
###############################################################################

ctrl_columns <- grep("C19_ctrl", colnames(counts_ordered), value = TRUE)
sev_columns <- grep("C19_sev", colnames(counts_ordered), value = TRUE)
ctrl_means <- rowMeans(counts_ordered[, ..ctrl_columns])
sev_means <- rowMeans(counts_ordered[, ..sev_columns])
labels <- ifelse(ctrl_means > sev_means, "control", "severe")

labels <- as.data.frame(labels)
rownames(labels) <- rownames(counts_ordered)
saveRDS(labels, file = "/icbb/projects/igunduz/mono_projection_070824/labels.rds")

###############################################################################
