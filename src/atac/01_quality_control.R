#!/usr/bin/env Rscript

#####################################################################
# 02_quality_control.R
# created on 2023-08-24 by Irem Gunduz
# Subset ArchR project and perform quality control
#####################################################################

## Load Libraries
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(ggplot2)
  library(BSgenome.Hsapiens.UCSC.hg38)
})
set.seed(12) # set seed

source("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/utils/archr_utils.R")

addArchRThreads(threads = 30) # set the cores
addArchRGenome("hg38") # set the reference genome

## Load and subset the ArchR project
outputDir <- "/icbb/projects/igunduz/DARPA/ArchRProject_5x/"

# load the project and cell annotation table
project <- ArchR::loadArchRProject(outputDir, force = T)

# exclude the BA samples
ba_cells <- data.table::fread("/icbb/projects/igunduz/DARPA/ArchRProject_5x/cellAnnot_archr.tsv") %>%
  dplyr::filter(!sample_exposure_type == "BA") %>%
  dplyr::select(Sample)
idx <- BiocGenerics::which(project$Sample %in% ba_cells$Sample)

# extract the cells of the samples
cellsSample <- project$cellNames[idx]
outputDir <- "/icbb/projects/igunduz/archr_project_011023/"

# subset the project to directory
project <- ArchR::subsetArchRProject(project,
  cells = cellsSample,
  outputDirectory = outputDir
)
project <- project[cellsSample, ]
arrows <- list.files(paste0(outputDir, "ArrowFiles"), full.names = TRUE)
project@sampleColData <- DataFrame(ArrowFiles = arrows)
rownames(project@sampleColData) <- gsub(x = gsub(x = arrows, ".arrow", ""), paste0(outputDir, "ArrowFiles/"), "")

# Get the row names
row_names <- rownames(project@sampleColData)
# Create a SimpleList with the row names
simple_list <- SimpleList(vector("list", length = 92))
names(simple_list) <- row_names

project@sampleMetadata <- simple_list
project@projectMetadata$outputDirectory <- outputDir
# project@projectMetadata$GroupCoverages$ClusterCellTypes$coverageMetadata$File = gsub("DARPA/ArchRProject_5x/ATAC_processed","ATAC_processed_final",project@projectMetadata$GroupCoverages$ClusterCellTypes$coverageMetadata$File)

# save the project to the subdirectory
saveRDS(project, paste0(outputDir, "Save-ArchR-Project.rds"))

## Add exposure annotations to cells
cellann <- data.table::fread("/icbb/projects/igunduz/DARPA/ArchRProject_5x/cellAnnot_archr.tsv") %>%
  dplyr::filter(!sample_exposure_type == "BA") %>%
  dplyr::filter(cellId_archr %in% project$cellNames) %>%
  dplyr::distinct(.keep_all = TRUE)

# add sample explosure type
project <- addCellColData(
  ArchRProj = project, data = cellann$sample_exposure_type,
  name = "sample_exposure_type", cells = project$cellNames, force = T
)
# add sample explosure group
project <- addCellColData(
  ArchRProj = project, data = cellann$sample_exposure_group,
  name = "sample_exposure_group", cells = project$cellNames, force = T
)

# add annotations
project <- addArchRAnnotations(project)

# before filtering
p <- plotUniqueFragsvsTSS(project)
# save pdf
plotPDF(p, name = "TSS-vs-Frags-before-filtering.pdf", ArchRProj = project, addDOC = FALSE)

# Filter based on TSS enrichment
idxPass <- which(project$TSSEnrichment >= 8)
cellsPass <- project$cellNames[idxPass]
project <- project[cellsPass, ]

# Filter based on unique fragments
idxPass <- which(log10(project$nFrags) >= 3)
cellsPass <- project$cellNames[idxPass]
project <- project[cellsPass, ]

# after filtering
p <- plotUniqueFragsvsTSS(project)
# save pdf
plotPDF(p, name = "TSS-vs-Frags-after-filtering.pdf", ArchRProj = project, addDOC = FALSE)

# remove doublets
project <- filterDoublets(project)

# save the project to the subdirectory
saveArchRProject(project, outputDirectory = outputDir, load = FALSE)

#####################################################################
# Plotting QC metrics
#####################################################################

# Define the color mapping
color_mapping <- c(
  "C19_Ctrl" = "#4F619D",
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

# Plot TSS enrichment
p2 <- plotGroups(
  ArchRProj = project,
  groupBy = "sample_exposure_group",
  colorBy = "cellColData",
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.8,
  addBoxPlot = TRUE
) +
  theme_classic()

p2 <- p2 + scale_fill_manual(values = color_mapping) + scale_color_manual(values = color_mapping)

# Plot number of unique fragments
p4 <- plotGroups(
  ArchRProj = project,
  groupBy = "sample_exposure_group",
  colorBy = "cellColData",
  name = "log10(nFrags)",
  plotAs = "violin",
  alpha = 0.8,
  addBoxPlot = TRUE
) +
  theme_classic()

p4 <- p4 + scale_fill_manual(values = color_mapping) + scale_color_manual(values = color_mapping) + theme(legend.position = "bottom")

# Save them into the project
plotPDF(p2, name = "QC-Sample-Statistics-NFragments.pdf", ArchRProj = project, addDOC = FALSE, width = 10, height = 10)
ggsave(p2, filename = "/icbb/projects/igunduz/QC-NFragments.pdf", units = "in", dpi = 300)
plotPDF(p4, name = "QC-Sample-Statistics-TSS-Enrichment.pdf", ArchRProj = project, addDOC = FALSE, width = 10, height = 10)
ggsave(p4, filename = "/icbb/projects/igunduz/QC-TSS-Enrichment-Final.pdf", width = 10, height = 10, units = "in", dpi = 300)


#####################################################################
