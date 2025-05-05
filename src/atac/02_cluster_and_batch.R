#!/usr/bin/env Rscript

#####################################################################
# 02_cluster_and_batch.R
# created on 2023-08-24 by Irem Gunduz
# Add iterative LSI, batch correction and cluster the cells
#####################################################################

## Load Libraries
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
})
set.seed(12) # set seed

addArchRThreads(threads = 30) # set the cores
addArchRGenome("hg38") # set the reference genome

## Load and subset the ArchR project
outputDir <- "/icbb/projects/igunduz/archr_project_011023/"

# load the project and cell annotation table
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)

# add iterative LSI
project <- addIterativeLSI(
  ArchRProj = project,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 2,
  clusterParams = list(
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:30,
  saveIterations = TRUE,
  nPlot = 100000
)

# add clusters
project <- addClusters(
  input = project,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  force = TRUE,
  name = "Clusters",
  resolution = 0.8
)

# add umap
project <- addUMAP(
  ArchRProj = project,
  reducedDims = "IterativeLSI",
  name = "UMAP",
  force = TRUE,
  nNeighbors = 40,
  minDist = 0.5,
  metric = "cosine"
)


## Apply batch correction using Harmony
# perform harmony
project <- addHarmony(
  ArchRProj = project,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "sample_exposure_type"
)

# add clusters
project <- addClusters(
  input = project,
  reducedDims = "Harmony",
  method = "Seurat",
  name = "Clusters_0.8",
  resolution = 0.8
)

# add umap
project <- addUMAP(
  ArchRProj = project,
  reducedDims = "Harmony",
  name = "UMAPHarmony",
  nNeighbors = 40,
  minDist = 0.5,
  metric = "cosine"
)

# plot UMAP based on sample and clusters
har_umaps <- plotEmbedding(
  ArchRProj = project,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAPHarmony"
)

har_umapc <- plotEmbedding(
  ArchRProj = project,
  colorBy = "cellColData",
  name = "Clusters_0.8",
  embedding = "UMAPHarmony"
)

har_umapsg <- plotEmbedding(
  ArchRProj = project,
  colorBy = "cellColData",
  name = "sample_exposure_type",
  embedding = "UMAPHarmony"
)

# save the UMAP plots
plotPDF(har_umapc, har_umapc, har_umapsg, har_umaps,
  name = "Harmony-Plots-Samples-Clusters-Tile-Matrix.pdf",
  ArchRProj = project, addDOC = FALSE, width = 5, height = 5
)


u1 <- plotEmbedding(
  ArchRProj = project,
  colorBy = "cellColData",
  name = "Clusters",
  embedding = "UMAP"
)

u2 <- plotEmbedding(
  ArchRProj = project,
  colorBy = "cellColData",
  name = "sample_exposure_type",
  embedding = "UMAP"
)

# save the UMAP plots
plotPDF(u1, u2,
  name = "Harmony-Plots-Samples-Clusters-Tile-Matrix-NO-HARMONY.pdf",
  ArchRProj = project, addDOC = FALSE, width = 5, height = 5
)

# save the project to the subdirectory
saveArchRProject(project, outputDirectory = outputDir, load = FALSE)

#####################################################################
