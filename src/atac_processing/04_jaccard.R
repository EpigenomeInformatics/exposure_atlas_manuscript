# load libraries
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(mclust)
  library(ComplexHeatmap)
  library(ggplot2)
  library(rescueR)
})

source("/icbb/projects/igunduz/exposure_atlas_manuscript/src/utils/jaccard.R")

project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)
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

# ra <- HeatmapAnnotation(atac_cluster=rownames(cM),col=list(atac_cluster=atac_label_cmap), which="row", show_legend=c("atac_cluster"=FALSE))
# ta <- HeatmapAnnotation(rna_cluster=colnames(cM),col=list(rna_cluster=rna_label_cmap), show_legend=c("rna_cluster"=FALSE))
hm <- BORHeatmap(
  cM,
  # limits=c(0,1),
  dataColorMidPoint = 0.4,
  # clusterCols=TRUE, clusterRows=TRUE,
  labelCols = TRUE, labelRows = TRUE,
  dataColors = whitePurple,
  # left_annotation = ra,
  # top_annotation = ta,
  showColDendrogram = F, # Should the column dendrogram be shown
  showRowDendrogram = F,
  row_names_side = "left",
  width = ncol(cM) * unit(0.5, "cm"),
  height = nrow(cM) * unit(0.5, "cm"),
  border_gp = gpar(col = "black") # Add a black border to entire heatmap
)
pdf(paste0("/icbb/projects/igunduz/jaccard_heatmap.pdf"), width = 6, height = 6)
draw(hm)
dev.off()
