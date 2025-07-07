#!/usr/bin/env Rscript

#####################################################################
# 10_tcells.R
# Created on 28-04-2025 by Irem B. Gunduz, written by Bei Wei
# Last modified on 28-04-2025
# This script is used to analyze T cells
#####################################################################

## Load Libraries
set.seed(42)
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ArchR)
  library(ggplot2)
  library(tidyr)
})

addArchRThreads(threads = 20) 

outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/"
echo_full <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)

# Get the cell data from echo
cell_data <- as.data.frame(echo_full@cellColData) %>%
  dplyr::filter(sample_exposure_type %in% c("HIV")) %>%
  dplyr::select(Sample, ClusterCellTypes) %>%
  dplyr::filter(ClusterCellTypes %in% c("T_mem_CD8"))

# Subset echo project for this cells
cellsSample <- rownames(cell_data)
project <- echo_full[cellsSample, ]

# Add IterativeLSI
projHeme2 <- addIterativeLSI(
  ArchRProj = project,
  useMatrix = "TileMatrix",
  name = "IterativeLSI_HIV",
  iterations = 3,
  sampleCellsPre = 400000,
  dimsToUse = 1:25,
  varFeatures = 50000,
  LSIMethod = 1,
  scaleDims = TRUE,
  corCutOff = 0.75,
  binarize = TRUE,
  selectionMethod = "var",
  scaleTo = 5000,
  totalFeatures = 500000,
  filterQuantile = 0.99,
  saveIterations = TRUE,
  UMAPParams = list(n_neighbors = 15, min_dist = 0.2, metric = "euclidean", verbose = FALSE, fast_sgd = TRUE),
  nPlot = 100000,
  threads = getArchRThreads(),
  seed = 1
)


projHeme2 <- addClusters(
  input = projHeme2,
  reducedDims = "IterativeLSI_HIV",
  method = "Seurat",
  name = "ClustersHIV",
  resolution = 0.8
)

projHeme2 <- addUMAP(
  ArchRProj = projHeme2,
  reducedDims = "IterativeLSI_HIV",
  name = "UMAP_HIV",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "euclidean"
)

project <- readRDS("/icbb/projects/igunduz/ArchR-Project-cd8t.rds")
fig_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/figures/"

# Create a vector for sample file names (matching order in your ArchR project)
sample_files <- c(
  "hiv6_fragments.tsv.gz", "hiv12_fragments.tsv.gz", "hiv9_fragments.tsv.gz",
  "hiv8_fragments.tsv.gz", "hiv4_fragments.tsv.gz", "hiv1_fragments.tsv.gz",
  "hiv2_fragments.tsv.gz", "hiv7_fragments.tsv.gz", "hiv3_fragments.tsv.gz",
  "hiv11_fragments.tsv.gz", "hiv10_fragments.tsv.gz", "hiv5_fragments.tsv.gz"
)

# Create a vector for the subject assignment corresponding to each sample
subject_info <- c(
  "sub1", "sub1", "sub1",
  "sub2", "sub2", "sub2",
  "sub3", "sub3", "sub3",
  "sub4", "sub4", "sub4"
)

# Create a vector for time points
time_point <- c(
  "pre", "acute", "chronic",
  "pre", "acute", "chronic",
  "pre", "acute", "chronic",
  "pre", "acute", "chronic"
)

# Add the time point information to the ArchR project
project <- addSampleColData(
  ArchRProj = project,
  data = time_point,
  name = "TimePoint",
  samples = sample_files,
  force = TRUE
)

# Add the subject information to the ArchR project
project <- addSampleColData(
  ArchRProj = project,
  data = subject_info,
  name = "Subject",
  samples = sample_files,
  force = TRUE # Set to TRUE to overwrite if the column already exists
)

df <- getEmbedding(project, embedding = "UMAP", returnDF = TRUE)
colnames(df) <- c("UMAP1", "UMAP2")

# Create a named vector where names are sample files and values are subjects
sample_to_subject <- c(
  "hiv6_fragments.tsv.gz" = "sub1", "hiv12_fragments.tsv.gz" = "sub1", "hiv9_fragments.tsv.gz" = "sub1",
  "hiv8_fragments.tsv.gz" = "sub2", "hiv4_fragments.tsv.gz" = "sub2", "hiv1_fragments.tsv.gz" = "sub2",
  "hiv2_fragments.tsv.gz" = "sub3", "hiv7_fragments.tsv.gz" = "sub3", "hiv3_fragments.tsv.gz" = "sub3",
  "hiv11_fragments.tsv.gz" = "sub4", "hiv10_fragments.tsv.gz" = "sub4", "hiv5_fragments.tsv.gz" = "sub4"
)
df$Sample <- sapply(strsplit(rownames(df), "#"), `[`, 1)
df$Subject <- sample_to_subject[df$Sample]

# Define the color palette for the subjects
colorPalette <- c(
  "sub1" = "#1f77b4", # blue for Subject 1
  "sub2" = "#ff7f0e", # orange for Subject 2
  "sub3" = "#7f7f7f", # gray for Subject 3
  "sub4" = "#ffbb00"
) # yellow for Subject 4

# Create the UMAP plot with the specified color palette
umap <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = Subject)) +
  geom_point() +
  scale_color_manual(values = colorPalette) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme_void() +
  theme(legend.position = "bottom")

ggsave(umap, file = paste0(fig_dir, "umap_tcell_subject.pdf"), width = 7, height = 7)


project <- readRDS("/icbb/projects/igunduz/ArchR-Project-cd8t.rds")
outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/"
# Reorgnaize arrow file path
arrows <- list.files(paste0(outputDir, "ArrowFiles"), full.names = TRUE)
project@sampleColData <- DataFrame(ArrowFiles = arrows)
rownames(project@sampleColData) <- arrows#gsub(x = gsub(x = arrows, ".arrow", ""), paste0(outputDir, "ArrowFile# Reorgnaize arrow file paths/"), "")

p1 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Sample", embedding = "UMAP")
p2 <- plotEmbedding(ArchRProj = project, colorBy = "cellColData", name = "Clusters", embedding = "UMAP")
ggAlignPlots(p1, p2, type = "h")
ggsave(p2, filename = paste0(fig_dir, "hiv_umaps.pdf"), width = 10, height = 5)


# Subject | Pre | Acute | Chronic
df <- data.frame(
  Subject = c("S1", "S2", "S3", "S4"),
  Pre = c(0, 0, 40, 40),
  Acute = c(1866453, 852580, 415480, 27250),
  Chronic = c(794, 110691, 244551, 12341)
)


# Reshape the data for ggplot and order the 'Stage' factor levels
df_long <- df %>%
  pivot_longer(cols = -Subject, names_to = "Stage", values_to = "Viral_Load") %>%
  mutate(Stage = factor(Stage, levels = c("Pre", "Acute", "Chronic")))


# Plotting
lp <- ggplot(df_long, aes(x = Stage, y = Viral_Load, color = Subject, group = Subject)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  labs(
    title = "Viral Load per Subject Across Stages",
    x = "Stage",
    y = "Viral Load"
  ) +
  theme_classic() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("S1" = "#1f77b4", "S2" = "#ff7f0e", "S3" = "#7f7f7f", "S4" = "#ffbb00"))
ggsave(lp, file = paste0(fig_dir, "viral_load_subject.pdf"), width = 7, height = 7)


# Exhaustion markers
markerGenes <- c(
  "HAVCR2", "CTLA4", "NCAM1", "ROBO2", "ROBO1", "TIGIT"
)

# Get the marker genes
markersGS <- getMarkerFeatures(
  ArchRProj = project,
  useMatrix = "GeneScoreMatrix",
  groupBy = "Clusters",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 1.25")

# Create the heatmap
heatmapGS <- markerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 0.01 & Log2FC >= 1.25",
  labelMarkers = markerGenes,
  transpose = TRUE
)

pdf(file = paste0(fig_dir, "exhaustion_markers_heatmap.pdf"), width = 10, height = 10)
draw(heatmapGS)
dev.off()

projHeme2 <- addImputeWeights(projHeme2)
p <- plotEmbedding(
  ArchRProj = projHeme2,
  colorBy = "GeneScoreMatrix",
  name = markerGenes,
  embedding = "UMAP_HIV",
  quantCut = c(0.01, 0.95)
)

p2 <- lapply(p, function(x) {
  x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
})
pdf(file = paste0(fig_dir, "exhaustion_markers_umap.pdf"), width = 10, height = 10)
do.call(cowplot::plot_grid, c(list(ncol = 3), p2))
dev.off()

#####################################################################
