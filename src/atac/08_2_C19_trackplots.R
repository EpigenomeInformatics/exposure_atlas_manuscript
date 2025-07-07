#!/usr/bin/env Rscript

#####################################################################
# 08_2_C19_trackplots.R
# created on 2025-05-05
# Track plots for C19 samples
#####################################################################

## Load Libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(Seurat)
  library(ArchR)
  library(chromVARmotifs)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(muLogR)
})
set.seed(12) # set seed
source("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/utils/helpers.R")

# Load the data
addArchRThreads(threads = 30) # set the cores
project <- ArchR::loadArchRProject("/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/", showLogo = FALSE)
bedfile <- paste0("/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/Mono_CD14/reports/differential_data/diffTab_3_archrPeaks.tsv")

# subset monocytes
idxSample <- BiocGenerics::which(project$ClusterCellTypes %in% c("Mono_CD14", "Mono_CD16"))
cellsSample <- project$cellNames[idxSample]
project <- project[cellsSample, ]

# subset C19 samples
idxSample <- BiocGenerics::which(project$sample_exposure_type %in% c("C19"))
cellsSample <- project$cellNames[idxSample]
project <- project[cellsSample, ]

DiffPeaks <- data.table::fread(bedfile)
DiffPeaks$isDiff <- cutL2FCpadj(DiffPeaks, lfc = 1, padj = 0.01)


# Filter for differential peaks and resize to extend upstream and downstream by 100000
DiffPeaks <- DiffPeaks %>%
  dplyr::filter(isDiff == TRUE) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# Get peaks from the ArchR project
peaks <- as.data.table(getPeakSet(project))
peaks <- as(peaks, "GRanges")


# Create a list of peaks and differential peaks
peaks_matched <- list(
  Peaks = peaks,
  DiffPeaks = DiffPeaks
)

region <- GenomicRanges::resize(DiffPeaks, width = 200000, fix = "center")

# Plot track with gene names, peaks, and differential peaks
plot <- plotBrowserTrack(
  ArchRProj = project,
  features = peaks_matched,
  minCells = 50,
  tileSize = 500,
  region = region, # Plot the top 5 regions (adjust as needed)
  ylim = c(0.001, 0.999),
  upstream = 50000,
  downstream = 50000,
  groupBy = "sample_exposure_group",
  useGroups = c("C19_ctrl", "C19_mild", "C19_mod", "C19_sev"),
  pal = c("#4F619D", "#8FBC8F", "#2C948F", "#006400")
)


# Save the track plot as PDF
plotPDF(
  plotList = plot,
  name = "Plot-Differential-Markers-Track-plots-log2FC1.pdf",
  ArchRProj = project,
  addDOC = FALSE, width = 5, height = 5
)


# set the markers
markerGenes <- c(
  "MZF1", "RUNX3", "SPIB", "REL", "RELB", "IGKC", "ATF3", "IL1B", "DOCK4",
  "RBM47", "IFI27", "PLSCR1", "STAT1", "OAS1", "DDX60", "PARP9", "DDX58", "TNFSF13B", "APOBEC3A",
  "SIGLEC1", "MX2", "OASL", "EIF2AK2", "OAS3", "IFIH1", "IRF7", "ISG15", "HERC5", "RSAD2", "IFIT1",
  "IFIT3", "IFIT2"
)

allDiffPeaks <- data.table::fread(bedfile)
allDiffPeaks$isDiff <- cutL2FCpadj(allDiffPeaks, lfc = 0.5, padj = 0.05)

# Filter for differential peaks and resize to extend upstream and downstream by 100000
allDiffPeaks <- allDiffPeaks %>%
  dplyr::filter(isDiff == TRUE) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)


# Create a list of peaks and differential peaks
peaks_matched <- list(
  Peaks = peaks,
  DiffPeaks = allDiffPeaks
)

# Plot track with gene names, peaks, and differential peaks
plot <- plotBrowserTrack(
  ArchRProj = project,
  features = peaks_matched,
  minCells = 50,
  tileSize = 500,
  geneSymbol = markerGenes,
  ylim = c(0.001, 0.999),
  upstream = 75000,
  downstream = 75000,
  groupBy = "sample_exposure_group",
  useGroups = c("C19_ctrl", "C19_mild", "C19_mod", "C19_sev"),
  pal = c("#4F619D", "#8FBC8F", "#2C948F", "#006400")
)

plotPDF(
  plotList = plot,
  name = "Plot-Differential-Markers-Track-plots-Mono.pdf",
  ArchRProj = project,
  addDOC = FALSE, width = 5, height = 5
)

##########################################################################
# genes from Wilk 2021
##########################################################################
# set the markers
genes <- c("FCGR1A",
  "HLA-E","CCL2", "CD4", "TNF", "IL6", "IL1B", "TXNIP", "CD74", "HLA-DRA",
   "HLA-DPA1", "HLA-DRB1",
  "HLA-DRB5", "HLA-DQB1", "EIF1", "ZFP36", "S100A8", "S100A9", "XAF1", "IFI6",
  "OAS1", "OAS2", "IFITM3", "LGALS1", "PLBD1", "S100A12", "FCGR3A", "CCL3",
  "CCL4", "HLA-DPB1", "HLA-DMA", "S100A9", "IRF7", "CXCL2", "ACTB", "CLU",
  "PTMA", "ITM3", "DYF", "MX1", "RNY1", "LGALS1", "PPIA"
)

allDiffPeaks <- data.table::fread(bedfile)
allDiffPeaks$isDiff <- cutL2FCpadj(allDiffPeaks, lfc = 0.5, padj = 0.05)

# Filter for differential peaks and resize to extend upstream and downstream by 100000
allDiffPeaks <- allDiffPeaks %>%
  dplyr::filter(isDiff == TRUE) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)


# Create a list of peaks and differential peaks
peaks_matched <- list(
  Peaks = peaks,
  DiffPeaks = allDiffPeaks
)

# Plot track with gene names, peaks, and differential peaks
plot <- plotBrowserTrack(
  ArchRProj = project,
  features = peaks_matched,
  minCells = 50,
  tileSize = 500,
  geneSymbol = genes,
  ylim = c(0.001, 0.999),
  upstream = 75000,
  downstream = 75000,
  groupBy = "sample_exposure_group",
  useGroups = c("C19_ctrl", "C19_mild", "C19_mod", "C19_sev"),
  pal = c("#4F619D", "#8FBC8F", "#2C948F", "#006400")
)

plotPDF(
  plotList = plot,
  name = "Plot-Wilk2021-Markers-Track-plots-Mono-HLAE.pdf",
  ArchRProj = project,
  addDOC = FALSE, width = 5, height = 5
)
##########################################################################

genes <- c(
  "HLA-E", "CCL2", "CD4", "TNF", "IL6", "IL1B", "TXNIP", "CD74", "HLA-DRA", "HLA-DPA1", "HLA-DRB1",
  "HLA-DRB5", "HLA-DQB1", "EIF1", "ZFP36", "S100A8", "S100A9", "XAF1", "IFI6",
  "OAS1", "OAS2", "IFITM3", "LGALS1", "PLBD1", "S100A12", "FCGR3A", "CCL3",
  "CCL4", "HLA-DPB1", "HLA-DMA", "S100A9", "IRF7", "CXCL2", "ACTB", "CLU",
  "PTMA", "ITM3", "DYF", "MX1", "RNY1", "LGALS1", "PPIA"
)
wimmers_genes <- c(
  "IGKC", "ATF3", "IL1B", "DOCK4", "RBM47", "IFI27", "PLSCR1", "STAT1",
  "OAS1", "DDX60", "PARP9", "DDX58", "TNFSF13B", "APOBEC3A", "SIGLEC1",
  "MX2", "OASL", "EIF2AK2", "OAS3", "IFIH1", "IRF7", "ISG15", "HERC5",
  "RSAD2", "IFIT1", "IFIT3", "IFIT2"
)

data("human_pwms_v2")
motifs <- human_pwms_v2
motif_name <- c("JUN", "FOS", "FOSL2", "NFKB1", "NFKB2", "RELA", "RELB","SPIB","CEBPD","FOSL2::JUN"
,"CREB1","BATF","BATF::JUN","SPIC")
# Combine multiple motif names into a single pattern
pattern <- paste(motif_name, collapse = "|")

# Find matching names in motifs
interesting <- names(motifs)[grepl(pattern, names(motifs), ignore.case = TRUE)]

motif_positions <- motifmatchr::matchMotifs(
  pwms = motifs[interesting],
  subject = getPeakSet(project),
  genome = eval(parse(text = getGenomeAnnotation(project)$genome)),
  out = "matches"
  #p.cutoff=1e-04
)

genes <- union(genes, wimmers_genes)
allDiffPeaks <- data.table::fread(bedfile)
allDiffPeaks$isDiff <- cutL2FCpadj(allDiffPeaks, lfc = 0.5, padj = 0.05)

# Filter for differential peaks and resize to extend upstream and downstream by 100000
allDiffPeaks <- allDiffPeaks %>%
  dplyr::filter(isDiff == TRUE) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)



for(motif in colnames(motif_positions)){
if(!file.exists(paste0("Plot-Motif-", motif, "-Track-plots-Mono.pdf"))){
# Create a list of peaks and differential peaks
peaks_matched <- list(
  Peaks = peaks,
  DiffPeaks = allDiffPeaks,
  Motif = getPeakSet(project)[assay(motif_positions)[,motif]]
)

# Plot track with gene names, peaks, and differential peaks
plot <- plotBrowserTrack(
  ArchRProj = project,
  features = peaks_matched,
  minCells = 50,
  tileSize = 500,
  geneSymbol = genes,
  ylim = c(0.001, 0.999),
  upstream = 75000,
  downstream = 75000,
  groupBy = "sample_exposure_group",
  useGroups = c("C19_ctrl", "C19_mild", "C19_mod", "C19_sev"),
  pal = c("#4F619D", "#8FBC8F", "#2C948F", "#006400")
)

plotPDF(
  plotList = plot,
  name = paste0("Plot-Motif-", motif, "-Track-plots-Mono.pdf"),
  ArchRProj = project,
  addDOC = FALSE, width = 5, height = 5
)
}
}
