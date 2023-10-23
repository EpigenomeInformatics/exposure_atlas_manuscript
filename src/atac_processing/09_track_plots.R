suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(data.table)
  library(GenomicRanges)
  library(chromVARmotifs)
  library(motifmatchr)
})
set.seed(12) # set seed
outputDir <- "/icbb/projects/igunduz/ATAC_processed_final/"
data("human_pwms_v2")


addArchRThreads(threads = 20) # set the cores
addArchRGenome("hg38") # set the reference genome
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)

# subset monocytes
idxSample <- BiocGenerics::which(project$ClusterCellTypes %in% "Mono_CD14")
cellsSample <- project$cellNames[idxSample]
project <- project[cellsSample, ]
bedfile <- paste0("/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/Mono_CD14/reports/differential_data/diffTab_3_archrPeaks.tsv")


diffs <- NULL
# extract the first differential peak
peaks <- data.table::fread(bedfile)
peaks$isDiff <- rescueR::cutL2FCpadj(peaks, padj = 0.01)
peaks <- peaks %>%
  dplyr::filter(isDiff == TRUE) %>%
  dplyr::arrange(padj) %>%
  GenomicRanges::makeGRangesFromDataFrame()
peaks <- GenomicRanges::resize(peaks, width = 50000, fix = "center")
diffs <- c(peaks, diffs)


# Convert GRanges to data.table
project_dt <- as.data.table(getPeakSet(project))

subset_dt <- project_dt[grepl("^Mono_CD14", GroupReplicate)]

# Convert back to GRanges object
subset_granges <- as(subset_dt, "GRanges")

motif_name <- c(
  "ENSG00000077150_LINE3188_NFKB2_D_N1", "ENSG00000109320_LINE3202_NFKB1_D_N3", "ENSG00000125347_LINE2729_IRF1_D_N2",
  "ENSG00000104856_LINE3201_RELB_D", "ENSG00000162924_LINE3217_REL_D_N3", "ENSG00000173039_LINE3227_RELA_D_N13"
)

motif_positions <- motifmatchr::matchMotifs(
  pwms = human_pwms_v2[motif_name],
  subject = subset_granges,
  genome = "hg38",
  out = "matches"
)
motif_names <- colnames(assay(motif_positions))
for (motif in motif_names) {
  peaks_matched <- list(
    `Peaks` = subset_granges,
    `Motif` = subset_granges[assay(motif_positions)[, motif]]
  )
  names(peaks_matched) <- recode(names(peaks_matched), `TRUE` = motif, `FALSE` = "no motif")
  plot <- plotBrowserTrack(project,
    features = peaks_matched,
    minCells = 50,
    tileSize = 500,
    region = diffs,
    ylim = c(0.001, 0.999),
    groupBy = "sample_exposure_group",
    useGroups = c("C19_ctrl", "C19_mild", "C19_mod", "C19_sev"),
    pal = c("#4F609C", "#C43E96", "#06948E", "#C03830")
  )
  # save them as pdf
  plotPDF(
    plotList = plot,
    name = paste0("Plot-differential-mono-track-plots", motif, ".pdf"),
    ArchRProj = project,
    addDOC = FALSE, width = 5, height = 5
  )
}

# set the markers
markerGenes <- c("MZF1", "RUNX3", "SPIB", "REL", "RELB")

# plot tracks
p <- plotBrowserTrack(
  ArchRProj = project,
  groupBy = "sample_exposure_group",
  geneSymbol = markerGenes,
  plotSummary = c("bulkTrack", "featureTrack", "loopTrack", "geneTrack"),
  minCells = 50,
  tileSize = 500,
  pal = c("#4F609C", "#C43E96", "#06948E", "#C03830"),
  useGroups = c("C19_ctrl", "C19_mild", "C19_mod", "C19_sev"),
  upstream = 70000,
  downstream = 70000
)

# save them as pdf
plotPDF(
  plotList = p,
  name = "Plot-motif_track_plots.pdf",
  ArchRProj = project,
  addDOC = FALSE, width = 5, height = 5
)
