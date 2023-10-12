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
source("/icbb/projects/igunduz/sc_epigenome_pathogen_exposure/utils/plots_utils.R")

cells <- data.table::fread("/icbb/projects/igunduz/DARPA/allc_sample_annot_final.csv") %>%
  dplyr::select(!V1) %>%
  dplyr::filter(!cell_type %in% c("Other-cell","Tc-Eff","Th-Eff"))

# get the unique cell types
cells <- names(table(cells$cell_type))

diffCompNames <- c(
  "C19_mild_vs_Ctrl", "C19_sev_vs_Ctrl", "HIV_acu_vs_Ctrl", "Influenza_ctrl_vs_d30",
  "HIV_chr_vs_Ctrl", "OP_high_vs_low", "OP_high_vs_med", "OP_low_vs_med"
)

plot.dir <- "/icbb/projects/igunduz/DARPA_analysis/rnbeads_run_061023/plots/"
if(!dir.exists(plot.dir)){dir.create(plot.dir)}
path <- "/icbb/projects/igunduz/DARPA_analysis/rnbeads_run_061023/"

for (cell in cells) {
  for (comp in diffCompNames) {
    pp <- rnbeadsDensityScatter(cell, comp, path, "archr_peaks")
    ggsave(plot = pp, paste0(plot.dir,"scatter_", cell, "_", comp, ".pdf"), width = 10, height = 10)
  }
}



#############################################################################################################################################################
## ATAC-METH Plots FOR COVID

# get the unique cell types
unique_cells <- as.data.frame(c("Monocyte", "NK-cell", "Tc-Mem", "Th-Mem", "B-cell"))
unique_cells$atac <- c("Mono_CD14", "NK_CD16", "T_mem_CD8", "T_mem_CD4", "B_comb")
colnames(unique_cells) <- c("meth", "atac")
unique_cells <- unique_cells[1, ]
diffCompNames <- diffCompNames[1:2]

# archr_peaks dataframe to the annotations
global_peaklist <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/regionsGR.rds") %>%
  data.table::as.data.table() %>%
  dplyr::select(!width & !strand) %>%
  as.data.frame()


for (i in 1:nrow(unique_cells)) {
  cell <- unique_cells[i, 1] # for methylation
  cells <- unique_cells[i, 2] # for atac
  outputDir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/" # atac directory

  for (comp in diffCompNames) {
    if (!file.exists(paste0(plot.dir,"L2FC_", cell, "_meth_atac_", comp, ".pdf"))) {
      j <- ifelse(comp == "C19_mild_vs_Ctrl", 1, 3)
      # ATAC get the differentially accessble regions
      atac <- prepareDARforPlot(cells, outputDir, j, "archrPeaks") 

      # Methylation, get the differntiall methylated regions
      dmr <- prepareDMRforPlot(cell, comp, path, global_peaklist, "archr_peaks")
      # dmt <-  prepareDMRforPlot(cell,comp,path,global_peaklist,"sites")

      # filter methylation by overlapped regions of atac
      datatable <- prepareScatterMethAtac(atac, dmr, threshold = 10)
      # datatable2 <- prepareScatterMethAtac(atac,dmt,threshold=5)
      cor <- cor(datatable$log2FC_1,datatable$log2FC_2) # -0.09375142
      # plot the scatter plot
      scplot <- plotScatterL2FC(datatable, y_lab = "mean-diff (METH)", x_lab = "log2FC (ATAC)", comb = "(ATAC vs METH)", group1 = "log2FC_1", group2 = "log2FC_2")
      # scplot2 <- plotScatterL2FC(datatable2,y_lab = "mean-diff (METH)",x_lab = "log2FC (ATAC)",comb ="(ATAC vs METH)",group1="log2FC_1",group2 ="log2FC_2")

      # save the plot
      ggsave(plot = scplot, paste0(plot.dir,"L2FC_", cell, "_meth_atac_", comp, ".pdf"), width = 20, height = 20)
      ggsave(plot = scplot, paste0(plot.dir,"L2FC_", cell, "_meth_atac_", comp, ".png"), width = 20, height = 20)
      # ggsave(plot=scplot2,paste0("/icbb/projects/igunduz/DARPA/Plots/L2FC_",cell,"_meth_atac_",comp,"sites.pdf"), width = 20, height = 20)
      # ggsave(plot=scplot2,paste0("/icbb/projects/igunduz/DARPA/Plots/L2FC_",cell,"_meth_atac_",comp,"sites.png"), width = 20, height = 20)
    }
  }
}

# lola meth
outputDir <- "/icbb/projects/igunduz/DARPA_analysis/rnbeads_run_061023/Monocyte/C19_sev_vs_Ctrl/reports/differential_methylation_data/differential_rnbDiffMeth/enrichment_lola.rds"
p <- lolaVolcanoPlotC19("Mono_CD14", lolaDb, outputDir, pValCut = 1.5, region = "archrPeaks", database = "TF_motif_clusters", signifCol = "qValue", cnames = c("loss", "gain"), n = 3, colorpanel = c())$plot
ggsave(plot = p, "/icbb/projects/igunduz/sc_epigenome_pathogen_exposure/figures/LOLA/TF_motif_clusters_METH.pdf", width = 20, height = 20)

p <- lolaVolcanoPlotC19("Mono_CD14", lolaDb, outputDir, pValCut = 1.5, region = "archrPeaks", database = "TF_motifs_vert", signifCol = "qValue", cnames = c("loss", "gain"), n = 3, colorpanel = c())$plot
ggsave(plot = p, "/icbb/projects/igunduz/sc_epigenome_pathogen_exposure/figures/LOLA/TF_motif_vert_METH.pdf", width = 20, height = 20)

#####################################################################
