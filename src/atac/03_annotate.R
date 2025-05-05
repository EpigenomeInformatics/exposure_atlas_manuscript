#!/usr/bin/env Rscript

#####################################################################
# 03_annotate.R
# created on 2023-08-24 by Irem Gunduz
# Annotate the cells based on the scRNAseq data
#####################################################################

set.seed(12) # set seed
outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/"
fig_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/figures/"

# load libraries
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(mclust)
  library(reshape2)
  library(plyr)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(pheatmap)
  library(ComplexHeatmap)
})

# load the functions
source("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/utils/helpers.R")
source("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/utils/archr_utils.R")

cbPalette <- c(
  "#ae017e", "#f768a1", "#67000d",
  "#fe9929", "#cc4c02", "#B5651D",
  "#a106bd", "#41b6c4", "#4292c6",
  "#0074cc", "#888fb5", "#c7e9b4"
)

addArchRThreads(threads = 30) # set the cores
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)

# load the scRNAseq data
sRNA <- readRDS("/icbb/projects/igunduz/DARPA/blood_imatlas_seurat.RDS")

# integration
project <- addGeneIntegrationMatrix(
  ArchRProj = project,
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = sRNA,
  addToArrow = FALSE,
  force = TRUE,
  groupRNA = "Manually_curated_celltype",
  nameCell = "predictedCellAtlas",
  nameGroup = "predictedGroupAtlas",
  nameScore = "predictedScoreAtlas"
)

# save the project to the subdirectory
saveArchRProject(project, outputDirectory = outputDir, load = T)

# Assuming project$predictedGroupAtlas is a factor variable
project$predictedGroupAtlas <- recode(project$predictedGroupAtlas,
  "ABCs" = "ABC",
  "Classical monocytes" = "Mono_CD14",
  "Cycling T&NK" = "Cycling_T_NK",
  "DC1" = "DC",
  "DC2" = "DC",
  "MAIT" = "T_mait",
  "Megakaryocytes" = "Megakaryocytes",
  "Memory B cells" = "B_mem",
  "Naive B cells" = "B_naive",
  "NK_CD16+" = "NK_CD16",
  "NK_CD56bright_CD16-" = "NK_CD16",
  "Nonclassical monocytes" = "Mono_CD16",
  "Plasma cells" = "Plasma",
  "Plasmablasts" = "Plasmablasts",
  "Progenitor" = "Progenitor",
  "T_CD4/CD8" = "T_mix",
  "Teffector/EM_CD4" = "T_mem_CD4",
  "Tem/emra_CD8" = "T_mem_CD8",
  "Tfh" = "T_mix",
  "Tgd_CRTAM+" = "T_mix",
  "Tnaive/CM_CD4" = "T_naive",
  "Tnaive/CM_CD8" = "T_naive",
  "Tregs" = "T_mix",
  "Trm_gut_CD8" = "T_mix",
  "Trm/em_CD8" = "T_mem_CD8"
)

# Check the updated levels
table(project$predictedGroupAtlas)

# compute ARI
rand <- mclust::adjustedRandIndex(
  project$predictedGroupAtlas,
  project$Clusters_0.8
)
# 0.4682788

# construct confusion matrix
df <- as.data.frame(as.matrix(confusionMatrix(
  project$predictedGroupAtlas,
  project$Clusters_0.8
)))

# plot Jaccard index as heatmap
jacch <- computeJaccardIndex(df, heatmap = TRUE)
pdf(file = paste0(outputDir, "Plots/jaccard_pheatmap_annotated.pdf"), width = 7, height = 7)
jacch
dev.off()

# Relabel the cell types
cM <- confusionMatrix(
  project$Clusters_0.8,
  project$predictedGroupAtlas
)
labelOld <- rownames(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew2 <- c(
  "Mono_CD14", "NK_CD16", "NK_CD16", "T_mem_CD8", "Mono_CD14", "Mono_CD16",
  "Mono_CD14", "T_mem_CD4", "T_mem_CD8", "T_mait", "T_mix", "B_naive",
  "T_mem_CD4", "T_naive", "T_mem_CD4", "T_naive", "B_mem", "T_mem_CD8",
  "DC", "Plasma", "T_mix", "B_naive"
)
project$ClusterCellTypes <- mapLabels(project$Clusters_0.8, oldLabels = labelOld, newLabels = labelNew2)


df <- getEmbedding(project, embedding = "UMAPHarmony", returnDF = TRUE)
colnames(df) <- c("UMAP1", "UMAP2")
df$Group <- project$ClusterCellTypes
umap <- plotNiceArchRumap(df, colorPalette = cbPalette)
ggsave(umap, file = paste0(fig_dir, "umap_annotated.pdf"), width = 7, height = 7)


# plot cell-type proportion plot
cell <- as.data.frame(project@cellColData) %>%
  dplyr::filter(!sample_exposure_group %in% c("BA_na", "BA_vac"))
cell <- table(cell$ClusterCellTypes, cell$sample_exposure_group)
cell <- as.matrix(cell)

stacked <- cellTypeProportionPlot(cell,
  scale = TRUE, center = FALSE,
  groupName = "Exposure", colorPalette = cbPalette, order = c(
    "C19_ctrl", "C19_mild", "C19_mod", "C19_sev",
    "HIV_ctrl", "HIV_acu", "HIV_chr", "Influenza_ctrl",
    "Influenza_d3", "Influenza_d6", "Influenza_d30",
    "OP_low", "OP_med", "OP_high"
  )
)
ggsave(stacked, file = paste0(fig_dir, "cellTypeProportion_stacked.pdf"), width = 13, height = 10)


# save the project to the subdirectory
saveArchRProject(project, outputDirectory = outputDir, load = FALSE)

#####################################################################
# Jaccard Index
#####################################################################

# construct confusion matrix
df <- as.data.frame(as.matrix(confusionMatrix(
  project$Clusters_0.8,
  project$predictedGroupAtlas
)))
cM <- computeJaccardIndex(df, heatmap = FALSE)
cM <- prettyOrderMat(t(cM), clusterCols = TRUE)$mat %>% t()

whitePurple <- c(
  "#f7fcfd", "#e0ecf4", "#bfd3e6", "#9ebcda", "#8c96c6",
  "#8c6bb1", "#88419d", "#810f7c", "#4d004b"
)

ht_opt$simple_anno_size <- unit(0.25, "cm")

# Plot the heatmap
hm <- BORHeatmap(
  cM,
  dataColorMidPoint = 0.4,
  labelCols = TRUE, labelRows = TRUE,
  dataColors = whitePurple,
  showColDendrogram = F,
  showRowDendrogram = F,
  row_names_side = "left",
  width = ncol(cM) * unit(0.5, "cm"),
  height = nrow(cM) * unit(0.5, "cm"),
  border_gp = gpar(col = "black")
)
pdf(paste0(fig_dir, "jaccard_heatmap.pdf"), width = 6, height = 6)
draw(hm)
dev.off()

#####################################################################
