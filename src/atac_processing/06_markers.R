#!/usr/bin/env Rscript

#####################################################################
# 05_markers.R
# created on 2023-08-24 by Irem Gunduz
# Add motif annotations and plot the UMAPs
#####################################################################

set.seed(12) # set seed
outputDir <- "/icbb/projects/igunduz/archr_project_011023/"

# load libraries
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(ChrAccR)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(JASPAR2018)
  library(chromVARmotifs)
  library(muLogR)
})
addArchRThreads(threads = 30) # set the cores
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)


# markers per cell type
markerGenes <- c(
  "CD34", # Early Progenitor
  "GATA1", # Erythroid
  "PAX5", "MS4A1", "MME", # B-Cell
  "CD14", "MPO", # Monocytes
  "IRF8",
  "CD3D", "CD8A", "CD4", # TCells
  "TBX21", # NK
  "CD3G", "NCAM1", "FCGR3A",
  "FOXP3", "GATA3", "RORC", "PDCD1",
  "HLA-DRA", "CD28", "IL2RA", "CD69", "CD44"
)

cbPalette <- c(
  "#ae017e", "#f768a1", "#67000d",
  "#fe9929", "#cc4c02", "#B5651D",
  "#a106bd", "#41b6c4", "#4292c6",
  "#888fb5", "#c7e9b4", "#0074cc"
)

# add impute weights
project <- addImputeWeights(project)

# plot GeneScoreMatrix UMAPs
p2 <- plotEmbedding(
  ArchRProj = project,
  colorBy = "GeneScoreMatrix",
  continuousSet = "horizonExtra",
  name = markerGenes,
  embedding = "UMAPHarmony"
)

# save the UMAP
plotPDF(
  plotList = p2,
  name = "Plot-UMAP-GeneScoreMatrix_final.pdf",
  ArchRProj = project,
  addDOC = FALSE, width = 5, height = 5
)

if ("Motif" %ni% names(project@peakAnnotation)) {
  project <- addMotifAnnotations(ArchRProj = project, motifSet = "cisbp", name = "Motif")
}

# add dev matrix for motif
project <- addDeviationsMatrix(
  ArchRProj = project,
  peakAnnotation = "Motif",
  force = TRUE
)

motifs <- c("GATA1", "GATA3", "RORC", "CEBPA", "EBF1", "IRF4", "TBX21", "PAX5", "FOXP3")
markerMotifs <- getFeatures(project, select = paste(motifs, collapse = "|"), useMatrix = "MotifMatrix")

# plot the motif activity
p1 <- plotEmbedding(
  ArchRProj = project,
  colorBy = "MotifMatrix",
  name = markerMotifs,
  embedding = "UMAPHarmony"
)

# save the UMAP
plotPDF(
  plotList = p1,
  name = "Plot-UMAP-chromVAR_final.pdf",
  ArchRProj = project,
  addDOC = FALSE, width = 5, height = 5
)

data("human_pwms_v2") # load the motifs
motifNames <- grep("PAX5|GATA3|GATA1|IRF4|TBX21|CD14", names(human_pwms_v2), value = TRUE, ignore.case = TRUE) # search NFKB motifs
markersConcatenated <- paste(markerGenes, collapse = "|")
motifNames <- grep(markersConcatenated, names(human_pwms_v2), value = TRUE, ignore.case = TRUE) # search NFKB motifs

motif_positions <- motifmatchr::matchMotifs(
  pwms = human_pwms_v2[motifNames],
  subject = getPeakSet(project),
  genome = "hg38",
  out = "matches"
)

for (motif in motifNames) {
  clean_motif <- sub("^.*_(\\w+)_D_N\\d+$", "\\1", motif)
  peaks_matched <- list(
    `All Peaks` = getPeakSet(project),
    `Motif` = getPeakSet(project)[assay(motif_positions)[, motif]]
  )
  names(peaks_matched) <- recode(names(peaks_matched), `TRUE` = motif, `FALSE` = "no motif")

  # plot track plots for markers
  p <- plotBrowserTrack(
    ArchRProj = project,
    groupBy = "ClusterCellTypes",
    geneSymbol = clean_motif,
    features = peaks_matched,
    upstream = 50000,
    downstream = 50000,
    pal = cbPalette
  )

  # save the track plots for markers
  plotPDF(
    plotList = p,
    name = paste0(clean_motif, "-plot-track.pdf"),
    ArchRProj = project,
    addDOC = FALSE, width = 5, height = 5
  )
}

saveArchRProject(project, outputDirectory = outputDir, load = FALSE)
