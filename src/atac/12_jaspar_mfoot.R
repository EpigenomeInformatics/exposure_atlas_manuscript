#!/usr/bin/env Rscript

#####################################################################
# archr_footprint.R
# created on 02-04-24 by Irem Gunduz
# Plot footprints for the marker motifs
#####################################################################


# load libraries
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(muLogR)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

set.seed(12) # set seed
outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023"
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)
addArchRThreads(threads = 30) # set the cores
#annot_acc <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/cell_matching/annot_acc_cellType.rds")

# Subset the Monocytes
idxSample <- BiocGenerics::which(project$ClusterCellTypes %in% c("Mono_CD14", "Mono_CD16"))
cellsSample <- project$cellNames[idxSample]
project <- project[cellsSample, ]

# Subset the covid
idxSample <- BiocGenerics::which(project$sample_exposure_group %in% c("C19_sev","C19_ctrl"))
cellsSample <- project$cellNames[idxSample]
project <- project[cellsSample, ]


# get motif positions
motifPositions <- getPositions(project,name="jaspar2020")

# Get the marker motifs
#motifs <- c("FOS", "RUNX3", "FOSL2", "BATF", "TBX21", "NFATC1", "IRF4", "EGR2", "STAT","REL","TCF4","SPI","FOXP1")
#markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
#markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs <- names(motifPositions)
nfkb_family <- c("NFKB1", "NFKB2", "RELA", "RELB", "REL")

# add group coverages
project <- addGroupCoverages(ArchRProj = project, groupBy = "sample_exposure_group",force = TRUE)

# get footprints
seFoot <- getFootprints(
  ArchRProj = project, 
  positions = motifPositions[nfkb_family], 
  groupBy = "sample_exposure_group"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = project, 
  normMethod = "divide",
  plotName = "C19-Footprints-Divide-Bias-nfkb_family",
  addDOC = FALSE
)
#saveArchRProject(project, outputDir)
