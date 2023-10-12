#!/usr/bin/env Rscript

#####################################################################
# 05_plot_cca.R
# created on 2023-08-26 written by Fabian Mueller adapted by Irem Gunduz
# Plot CCA results
#####################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexHeatmap)
  library(muLogR)
  library(muRtools)
  library(GenomicRanges)
  library(SummarizedExperiment)
  library(ChrAccR)
  library(parallel)
  library(dplyr)
  library(Seurat)
})
set.seed(12) # set seed
k <- "sub11kpc50"
oDir <- file.path("/icbb/projects/igunduz/DARPA_analysis/artemis_031023", "cell_matching")
fn_anch <- file.path(oDir, paste0("ccaCellMatches_", "anchors", ".rds"))
fn_nns <- file.path(oDir, paste0("ccaCellMatches_", "nns", ".rds"))

logger.start("Load anchor data..")
matchRes_anchors <- readRDS(fn_anch)
matchRes_nns <- readRDS(fn_nns)
logger.completed()

logger.start("Load annotation")
cellAnnot_atac <- read.delim("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/cellAnnot_atac.tsv")
rownames(cellAnnot_atac) <- cellAnnot_atac[, "cellId_archr"]
cellAnnot_atac$Cluster_atac <- cellAnnot_atac$Clusters_0.8
cellAnnot_atac$exposure <- cellAnnot_atac$sample_exposure_group
cellAnnot_atac$exposure <- ifelse(cellAnnot_atac$sample_exposure_group =="HIV_acu", "HIV_acute",cellAnnot_atac$sample_exposure_group)
cellAnnot_atac$exposure <- ifelse(cellAnnot_atac$sample_exposure_group =="HIV_chr", "HIV_chronic",cellAnnot_atac$sample_exposure_group)
cellAnnot_atac$exposure <- ifelse(cellAnnot_atac$sample_exposure_group =="HIV_ctr", "HIV_pre",cellAnnot_atac$sample_exposure_group)
cellAnnot_atac$exposure <- ifelse(cellAnnot_atac$sample_exposure_group =="OP_med", "OP_medium",cellAnnot_atac$sample_exposure_group)

cellAnnot_meth <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/cellAnnot_meth.rds")
rownames(cellAnnot_meth) <- cellAnnot_meth[, "Cell_UID"]
cellAnnot_meth$exposure <- cellAnnot_meth$condition
logger.completed()

logger.start("Plotting")
# clean-up contingency tables by omitting entries lower than a certain percentage
contTab_truncPerc <- function(ct, perc = 0.05, by = "row", removeZeroCols = TRUE) {
  if (by == "row") {
    percM <- ct / rowSums(ct)
  } else if (by == "column") {
    percM <- t(t(ct) / colSums(ct))
  } else {
    logger.error("Unknown by")
  }
  ct[percM < perc] <- 0L
  if (removeZeroCols) ct <- ct[, colSums(ct, na.rm = TRUE) > 0]
  return(ct)
}
colorSchemes <- list(
  ClusterCellTypes = c(
    B_mem = "#AE017E",
    B_naive = "#F768A1",
    DC =  "#67000D",
    Mono_CD14 = "#FE9929",
    Mono_CD16 = "#CC4C02",
    NK_CD16 = "#A65628",
    Plasma = "#A106BD",
    T_mait = "#41B6C4",
    T_mem_CD4 = "#4292c6",
    T_mem_CD8 = "#0074cc",
    T_mix = "#888FB5",
    T_naive = "#C7E9B4"
  ),
  cell_type = c(
    "B-cell" = "#AE017E",
    "Monocyte" = "#CC4C02",
    "NK-cell" = "#A65628",
    "Th-Mem" = "#41B6C4",
    "Tc-Mem" = "#4292C6",
    "Tc-Naive" = "#888FB5",
    "Th-Naive" = "#C7E9B4",
    "Other-cell" = "#CCCCCC"
  ),
  Clusters = c(
    c0 = "#FF0000",
    c1 = "#0000FF",
    c2 = "#008000",
    c3 = "#FFA500",
    c4 = "#800080",
    c5 = "#00FFFF",
    c6 = "#FFC0CB",
    c7 = "#A52A2A",
    c8 = "#808080",
    c9 = "#FF00FF",
    c10 = "#FFFF00",
    c11 = "#006400",
    c12 = "#00008B"
  ),
  exposure = c(
    HIV_acute = "#DC143C",
    HIV_chronic = "#800080",
    HIV_pre =  "#FFC0CB",
    OP_high = "#A0522D",
    OP_low = "#5E3D23"
  ),
  Cluster_atac = c(
    C0 = "#FF0000",
    C1 = "#0000FF",
    C2 = "#008000",
    C3 = "#FFA500",
    C4 = "#800080",
    C5 = "#00FFFF",
    C6 = "#FFC0CB",
    C7 = "#A52A2A",
    C8 = "#808080",
    C9 = "#FF00FF",
    C10 = "#FFFF00",
    C11 = "#006400",
    C12 = "#00008B",
    C13 = "#00FF00",
    C14 = "#FFD700",
    C15 = "#008080",
    C16 = "#EE82EE",
    C17 = "#4B0082",
    C18 = "#000080",
    C19 = "#FA8072",
    C20 = "#40E0D0",
    C21 = "#E6E6FA"
  )
)

annotPairs <- data.frame(
  annotName = c("cellType", "cluster", "exposure"),
  annot_atac = c("ClusterCellTypes", "Cluster_atac", "exposure"),
  annot_meth = c("cell_type", "Clusters", "exposure"),
  stringsAsFactors = FALSE
)
rownames(annotPairs) <- annotPairs[, "annotName"]

clName <- "cellType"
cs_cl_atac <- colorSchemes[[annotPairs[clName, "annot_atac"]]]
cs_cl_meth <- colorSchemes[[annotPairs[clName, "annot_meth"]]]


logger.start("Nearest neighbor matching")
matchRes <- matchRes_nns[matchRes_nns[, "mapping"] == "nearest_meth_for_acc", ]
annot_acc <- data.frame(class = cellAnnot_atac[, annotPairs[clName, "annot_atac"]])
rownames(annot_acc) <- rownames(cellAnnot_atac)
annot_acc[, "mapped_cell"] <- sapply(rownames(annot_acc), FUN = function(cid) {
  res <- NA
  X <- matchRes[matchRes[, "cell_acc"] == cid, , drop = FALSE]
  if (nrow(X) > 0) {
    res <- X[which.min(X[, "dist"]), "cell_meth"]
  }
  return(res)
})

# remove missing ones
annot_acc <- annot_acc[complete.cases(annot_acc[, "mapped_cell"]), ]

# map cells
annot_acc[, "mapped_cell_class"] <- cellAnnot_meth[annot_acc[, "mapped_cell"], annotPairs[clName, "annot_meth"]]

logger.info("Plotting chord diagrams")
contTab_a2m <- table(annot_acc[, "class"], annot_acc[, "mapped_cell_class"])
names(dimnames(contTab_a2m)) <- c("A", "M")
fn <- file.path(oDir, paste0("chord_nn_match", "_nearest_meth_for_acc", ".pdf"))
pdf(fn, width = 10, height = 10)
par(cex = 0.45) # Adjust the font size
chordDiagramFromContingencyTable(contTab_truncPerc(contTab_a2m), chordColorByCol = TRUE, cs_row = cs_cl_atac, cs_column = cs_cl_meth)
dev.off()

matchRes <- matchRes_nns[matchRes_nns[, "mapping"] == "nearest_acc_for_meth", ]
annot_meth <- data.frame(class = cellAnnot_meth[, annotPairs[clName, "annot_meth"]])
rownames(annot_meth) <- rownames(cellAnnot_meth)
annot_meth[, "mapped_cell"] <- sapply(rownames(annot_meth), FUN = function(cid) {
  res <- NA
  X <- matchRes[matchRes[, "cell_meth"] == cid, , drop = FALSE]
  if (nrow(X) > 0) {
    res <- X[which.min(X[, "dist"]), "cell_acc"]
  }
  return(res)
})
# remove missing ones
annot_meth <- annot_meth[complete.cases(annot_meth[, "mapped_cell"]), ]

annot_meth[, "mapped_cell_class"] <- cellAnnot_atac[annot_meth[, "mapped_cell"], annotPairs[clName, "annot_atac"]]
# chord diagram
contTab_m2a <- table(annot_meth[, "class"], annot_meth[, "mapped_cell_class"])
names(dimnames(contTab_m2a)) <- c("M", "A")
fn <- file.path(oDir, paste0("chord_nn_match", "_nearest_acc_for_meth", ".pdf"))
pdf(fn, width = 10, height = 10)
par(cex = 0.45) # Adjust the font size
chordDiagramFromContingencyTable(contTab_truncPerc(contTab_m2a,by="column"), chordColorByCol = TRUE, cs_row = cs_cl_meth, cs_column = cs_cl_atac)
dev.off()

logger.start("Contingency Tables")
# row-scale tables
contTab_a2ms <- contTab_a2m / rowSums(contTab_a2m)
class(contTab_a2ms) <- "matrix"
contTab_m2as <- t(contTab_m2a / rowSums(contTab_m2a))
class(contTab_m2as) <- "matrix"

# Remove the "Plasma" row from contTab_a2ms
#contTab_a2ms <- contTab_a2ms[rownames(contTab_a2ms) %in% rownames(contTab_m2as), ]

rowAnnot <- rowAnnotation(
  class_atac = rownames(contTab_a2ms),
  col = list(
    class_atac = cs_cl_atac
  ),
  show_annotation_name = FALSE
)
colAnnot <- HeatmapAnnotation(
  class_meth = colnames(contTab_a2ms),
  col = list(
    class_meth = cs_cl_meth
  ),
  show_annotation_name = FALSE
)

cs <- colpal.cont(9, "cb.Reds")
cs <- circlize::colorRamp2(seq(0, 1, length.out = length(cs)), cs)
# cs <- NULL
hm <- diagDivCellHeatmap(
  contTab_a2ms, contTab_m2as,
  col.l = cs, col.r = cs,
  name.l = "ATAC <- METH", name.r = "METH <- ATAC",
  left_annotation = rowAnnot,
  top_annotation = colAnnot,
  cluster_rows = getClusteringDendrogram(t(contTab_a2ms), distMethod = "euclidean", linkMethod = "ward.D", corMethod = "pearson"),
  cluster_columns = getClusteringDendrogram(contTab_m2as, distMethod = "euclidean", linkMethod = "ward.D", corMethod = "pearson")
)
fn <- file.path(oDir, paste0("contTab_nn_match", "_both", ".pdf"))
pdf(fn, width = 10, height = 10, onefile = FALSE)
draw(hm)
dev.off()
logger.completed()
logger.completed()


logger.start("Preapaing UMAP data")
umap_atac <- cellAnnot_atac[rownames(cellAnnot_atac) %in% rownames(annot_acc), paste0("UMAP_", 1:2)]
umap_meth <- readRDS(paste0("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/itLSI_res/itLSI_res_",k,".rds"))$umapCoord
umap_meth <- umap_meth[rownames(umap_meth) %in% rownames(annot_meth), ]
logger.completed()

logger.start("Side-by-side dimred plot (crosscluster)")
pp_acc <- getDimRedPlot(umap_atac, annot = annot_acc, colorCol = "mapped_cell_class", shapeCol = FALSE, colScheme = cs_cl_meth, ptSize = 0.25, addLabels = FALSE, addDensity = FALSE, annot.text = NULL) + coord_fixed() + guides(colour = guide_legend(override.aes = list(size = 5))) + theme_void()+ theme(legend.position = "bottom")
pp_meth <- getDimRedPlot(umap_meth, annot = annot_meth, colorCol = "mapped_cell_class", shapeCol = FALSE, colScheme = cs_cl_atac, ptSize = 0.25, addLabels = FALSE, addDensity = FALSE, annot.text = NULL) + coord_fixed() + guides(colour = guide_legend(override.aes = list(size = 5))) + theme_void() + theme(legend.position = "bottom")
pp <- cowplot::plot_grid(pp_acc, pp_meth)
ggsave4doc(file.path(oDir, paste0("umap_nearestNeighbor_crossmap_", "crosscluster", ".pdf")), pp, width = 384, height = 128)
logger.completed()

logger.start("Plotting other UMAPs")
logger.info("Arranging annotation")
annot_acc$sample_exposure_type <- cellAnnot_atac[rownames(annot_acc), "sample_exposure_group"]
annot_acc$sample_exposure_type <- cellAnnot_atac[rownames(annot_acc), "sample_exposure_group"]
annot_acc$ClusterCellTypes <- cellAnnot_atac[rownames(annot_acc), "ClusterCellTypes"]
annot_meth$sample_exposure_type <- cellAnnot_meth[rownames(cellAnnot_meth) %in% rownames(annot_meth), "condition"]
clustAss <- readRDS(paste0("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/itLSI_res/itLSI_res_",k,".rds"))$clustAss
clustAss <-  as.data.frame(clustAss)
annot_meth$Clusters <- clustAss[rownames(clustAss) %in% rownames(annot_meth),"clustAss"]
annot_acc$Cluster_atac <- cellAnnot_atac[rownames(annot_acc), "Clusters_0.8"]
annot_meth$cell_type <- cellAnnot_meth[rownames(cellAnnot_meth) %in% rownames(annot_meth), "cell_type"]
name_mapping <- c(
  "HIV_acu" = "HIV_acute",
  "HIV_chr" = "HIV_chronic",
  "HIV_ctrl" = "HIV_pre",
  "OP_high" = "OP_high",
  "OP_low" = "OP_low"
)

# Rename the categories using the mapping
annot_acc$sample_exposure_type <- name_mapping[annot_acc$sample_exposure_type]
annot_acc$exposure <- annot_acc$sample_exposure_type
annot_meth$exposure <- annot_meth$sample_exposure_type
umapDir <- file.path("/icbb/projects/igunduz/DARPA_analysis/artemis_031023", "umap")
if (!dir.exists(umapDir)) dir.create(umapDir)

	for (i in 1:nrow(annotPairs)){
		cs <- NULL
		if (is.element(annotPairs[i, "annot_atac"], names(colorSchemes))) cs <- colorSchemes[[annotPairs[i, "annot_atac"]]]
		pp_atac <- getDimRedPlot(umap_atac, annot=annot_acc, colorCol=annotPairs[i, "annot_atac"], shapeCol=FALSE, colScheme=cs, ptSize=0.25, addLabels=FALSE, addDensity=FALSE, annot.text=NULL) + coord_fixed() + theme_void()
		if (!is.numeric(annot_acc[,annotPairs[i, "annot_atac"]])){
			pp_atac <- pp_atac + guides(colour=guide_legend(override.aes=list(size=5)))+ theme(legend.position = "bottom")
		}
		cs <- NULL
		if (is.element(annotPairs[i, "annot_meth"], names(colorSchemes))) cs <- colorSchemes[[annotPairs[i, "annot_meth"]]]
		pp_meth <- getDimRedPlot(umap_meth, annot=annot_meth, colorCol=annotPairs[i, "annot_meth"], shapeCol=FALSE, colScheme=cs, ptSize=0.25, addLabels=FALSE, addDensity=FALSE, annot.text=NULL) + coord_fixed() + theme_void()
		if (!is.numeric(annot_meth[,annotPairs[i, "annot_meth"]])){
			pp_meth <- pp_meth + guides(colour=guide_legend(override.aes=list(size=5)))+ theme(legend.position = "bottom")
		}
		pp <- cowplot::plot_grid(pp_atac, pp_meth)
		ggsave4doc(file.path(umapDir, paste0("dimRed_umap_", annotPairs[i, "annotName"], ".pdf")), pp, width=384, height=128)
	}
logger.completed()


