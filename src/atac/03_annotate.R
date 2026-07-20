#!/usr/bin/env Rscript

#####################################################################
# 03_annotate.R
# created on 2023-08-24 by Irem Gunduz
# Annotate the cells based on the scRNAseq data
#
# REVISION (reviewer R2.3/R2.4): added export of a per-cluster
# annotation mapping table (cluster -> final label, top predicted
# scRNA label, Jaccard similarity, ambiguity flag) as a TSV for the
# supplementary tables.
#####################################################################

set.seed(12) # set seed
outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/"
fig_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/figures/"
# directory for supplementary table output
sannot_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/sample_annots/"

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
# REVISION (R2.3/R2.4): export per-cluster annotation mapping table
#
# Builds one row per graph-based cluster (Clusters_0.8) with:
#   - final assigned cell-type label (ClusterCellTypes / labelNew2)
#   - top predicted scRNA label (highest confusion-matrix count)
#   - Jaccard similarity of the cluster to its top predicted label
#   - number of cells in the cluster
#   - an ambiguity flag (top Jaccard below a threshold, or top-2 close)
#####################################################################

# 1. Confusion matrix: rows = clusters (Clusters_0.8), cols = predicted scRNA groups
cM_counts <- confusionMatrix(
  project$Clusters_0.8,
  project$predictedGroupAtlas
)
cM_counts <- as.matrix(cM_counts)

# 2. Jaccard index matrix 
jacc_mat <- as.matrix(
  computeJaccardIndex(as.data.frame(cM_counts), heatmap = FALSE)
)

# FIX: Transpose jacc_mat if the rows/cols are inverted compared to cM_counts
if (nrow(jacc_mat) == ncol(cM_counts) && ncol(jacc_mat) == nrow(cM_counts)) {
  jacc_mat <- t(jacc_mat)
}

# align jaccard matrix to the confusion matrix orientation (clusters as rows)
if (!all(rownames(jacc_mat) == rownames(cM_counts))) {
  jacc_mat <- jacc_mat[rownames(cM_counts), colnames(cM_counts), drop = FALSE]
}

# 3. Map final labels (labelNew2) to clusters, in the row order of cM_counts
final_label_map <- setNames(labelNew2, labelOld)   # labelOld = rownames(cM) earlier

# 4. Assemble per-cluster rows
clusters <- rownames(cM_counts)
mapping <- lapply(clusters, function(cl) {
  counts_row <- cM_counts[cl, ]
  jacc_row   <- jacc_mat[cl, ]

  # top predicted scRNA label by Jaccard (fall back to counts if ties)
  ord_j <- order(jacc_row, decreasing = TRUE)
  top1_label <- colnames(cM_counts)[ord_j[1]]
  top1_jacc  <- round(jacc_row[ord_j[1]], 4)
  top2_label <- colnames(cM_counts)[ord_j[2]]
  top2_jacc  <- round(jacc_row[ord_j[2]], 4)

  n_cells <- sum(counts_row)
  final   <- final_label_map[cl]

  # ambiguity flag: top Jaccard low, or top-2 Jaccard within 0.1 of each other
  ambiguous <- (top1_jacc < 0.2) | ((top1_jacc - top2_jacc) < 0.1)

  data.frame(
    Cluster              = cl,
    Final_CellType       = unname(final),
    Top_predicted_scRNA  = top1_label,
    Top_Jaccard          = top1_jacc,
    Second_predicted_scRNA = top2_label,
    Second_Jaccard       = top2_jacc,
    N_cells              = n_cells,
    Ambiguous            = ifelse(ambiguous, "yes", ""),
    stringsAsFactors     = FALSE
  )
})
mapping_df <- do.call(rbind, mapping)

# order by final cell type then cluster for readability
mapping_df <- mapping_df[order(mapping_df$Final_CellType, mapping_df$Cluster), ]

# 5. Write TSV for the supplementary tables
write.table(
  mapping_df,
  file = file.path(sannot_dir, "cluster_annotation_mapping.tsv"),
  sep = "\t", row.names = FALSE, quote = FALSE
)
message("Wrote cluster annotation mapping: ",
        file.path(sannot_dir, "cluster_annotation_mapping.tsv"))
print(mapping_df)

#####################################################################