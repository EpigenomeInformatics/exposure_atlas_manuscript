#!/usr/bin/env Rscript

#####################################################################
# 03_RnBeads_plots.R
# created on 2023-08-24 by Irem Gunduz
# Generate plots for RnBeads results
#####################################################################

set.seed(42)
suppressPackageStartupMessages({
  library(RnBeads)
  library(dplyr)
  library(ChrAccR)
  library(ggplot2)
  library(data.table)
})

# Load the functions
source("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/utils/plot_utils.R")
source("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/utils/chraccr_plots.R")

cells <- data.table::fread("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/sample_annots/allc_sample_annot_final.csv") %>%
  dplyr::select(!V1) %>%
  dplyr::filter(!cell_type %in% c("Other-cell", "Tc-Eff", "Th-Eff"))

# get the unique cell types
cells <- names(table(cells$cell_type))

diffCompNames <- c(
  "C19_mild_vs_Ctrl", "HIV_acu_vs_Ctrl", "Influenza_ctrl_vs_d30",
  "HIV_chr_vs_Ctrl", "OP_high_vs_low", "OP_high_vs_med", "OP_low_vs_med"
)

plot.dir <- "/icbb/projects/igunduz/DARPA_analysis/RnBeads_0111023/plots_140225/"
if (!dir.exists(plot.dir)) {
  dir.create(plot.dir)
}
path <- "/icbb/projects/igunduz/DARPA_analysis/RnBeads_0111023/"

for (cell in cells) {
  for (comp in diffCompNames) {
    pp <- rnbeadsDensityScatter(cell, comp, path, "archr_peaks")
    ggsave(plot = pp, paste0(plot.dir, "scatter_", cell, "_", comp, ".pdf"), width = 10, height = 10)
  }
}


#####################################################################
# C19 plots
#####################################################################


# archr_peaks dataframe to the annotations
global_peaklist <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/regionsGR.rds") %>%
  data.table::as.data.table() %>%
  dplyr::select(!width & !strand) %>%
  as.data.frame()

outputDir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/" # atac directory
analysis.dir <- "/icbb/projects/igunduz/DARPA_analysis/rnbeads_C19_monosub_080425_V2/Monocyte/C19_sev_vs_Ctrl/"
report.dir <- file.path(analysis.dir, "reports")
plot.dir <- "/icbb/projects/igunduz/exposure_atlas_manuscript/src/meth_processing/plots_270425/"
if (!dir.exists(plot.dir)) {
  dir.create(plot.dir)
}


logger.info("Loading RnBeads objects")
diffMeth <- load.rnb.diffmeth(paste0(analysis.dir, "/reports/differential_methylation_data/differential_rnbDiffMeth/"))
rnb_set <- load.rnb.set(paste0(analysis.dir, "/reports/data_import_data/rnb.set_preprocessed/"))

# ATAC get the differentially accessble regions
atac <- prepareDARforPlot("Mono_CD14", outputDir, 3, "archrPeaks")

# Get differentially methylated regions
dmrs <- prepareDMRforPlot(rnb_set, diffMeth, global_peaklist, "archr_peaks", changeMethod = "meandiff")

for (i in 1:length(dmrs)) {
  # filter methylation by overlapped regions of atac
  datatable <- prepareScatterMethAtac(atac, dmrs[[i]], threshold = 10)
  cor <- cor(datatable$log2FC_1, datatable$log2FC_2) # -0.07179112
  print(paste0("Correlation in ", names(dmrs)[i], " is ", cor))
  plot.dir <- "/icbb/projects/igunduz/exposure_atlas_manuscript/src/meth_processing/Fig_5_sub_080425/plots/"
  scplot <- plotScatterL2FC(datatable, y_lab = "mean-diff (METH)", x_lab = "log2FC (ATAC)", comb = "(ATAC vs METH)", group1 = "log2FC_1", group2 = "log2FC_2")
  ggsave(plot = scplot, paste0(plot.dir, "L2FC_Mono_CD14_meth_atac_", names(dmrs)[i], ".pdf"), width = 20, height = 20)
}


# Plot the density scatter plot for comparasions
plot_path <- "/icbb/projects/igunduz/exposure_atlas_manuscript/src/meth_processing/Fig_5_sub_080425/plots/"
rnbeadsDensityScatter_sub(diffMeth, "archr_peaks", plot_path)

# load loladb database
lolaDb <- "/icbb/projects/share/annotations/lolaDB/hg38/"
lolaDb <- RnBeads::loadLolaDbs(lolaDb)

# LOLA for scMETH
outputDir <- "/icbb/projects/igunduz/DARPA_analysis/rnbeads_C19_monosub_080425_V2/Monocyte/C19_sev_vs_Ctrl/reports/differential_methylation_data/differential_rnbDiffMeth/lola_results.rds"
p <- lolaVolcanoPlotC19("Mono_CD14", lolaDb, outputDir, pValCut = 1.5, region = "archrPeaks", database = "TF_motif_clusters", signifCol = "qValue", cnames = c("loss", "gain"), n = 3, colorpanel = c())$plot
path <- paste0(plot.dir, "Mono_CD14_lolaVolcanoPlot_clustered.pdf")
ggsave(plot = p, filename = path, width = 20, height = 20)

p <- lolaVolcanoPlotC19("Mono_CD14", lolaDb, outputDir, pValCut = 1.5, region = "archrPeaks", database = "TF_motifs_vert", signifCol = "qValue", cnames = c("loss", "gain"), n = 3, colorpanel = c())$plot
path <- paste0(plot.dir, "Mono_CD14_lolaVolcanoPlot_motifs.pdf")
ggsave(plot = p, filename = path, width = 20, height = 20)

#####################################################################
