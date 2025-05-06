#!/usr/bin/env Rscript

#####################################################################
# 05_cca.R
# created on 2023-08-25 written by Fabian Mueller adapted by Irem Gunduz
# CCA analysis of scATAC + scMethylation data of the common samples
#####################################################################
# !/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(muLogR)
  library(muRtools)
  library(GenomicRanges)
  library(SummarizedExperiment)
  library(ChrAccR)
  library(Matrix)
  library(DelayedArray)
  library(HDF5Array)
  library(ArchR)
  library(parallel)
  library(Seurat)
  library(dplyr)
})
set.seed(12) # set seed
# outputDir <- "/icbb/projects/igunduz/archr_project_011023/"

logger.start("Loading ATAC data")
project <- ArchR::loadArchRProject("/icbb/projects/igunduz/archr_project_011023", showLogo = FALSE)
# project <- ArchR::loadArchRProject(outputDir, force = T,showLogo=FALSE)
methDir <- "/icbb/projects/igunduz/DARPA_analysis/artemis_031023/"

logger.start("Preparing sample annotation")
sampleAnnotAll <- readTab("/icbb/projects/igunduz/DARPA/samples_all_20221026.tsv")
rownames(sampleAnnotAll) <- sampleAnnotAll[, "CommonMinID"]
logger.completed()


# get the colData of the project
ca_ha <- as.data.frame(getCellColData(project))
logger.status("Annotation")
cellAnnot_atac <- ca_global <- readr::read_tsv("/icbb/projects/igunduz/DARPA/ArchRProject_5x/cellAnnot_archr.tsv")

# cell annotation
df <- project@embeddings$UMAPHarmony$df
colnames(df) <- c("UMAP_1", "UMAP_2")
df$cellId_archr <- rownames(df)
cellAnnot_atac <- cellAnnot_atac[, !colnames(cellAnnot_atac) %in% c("UMAP_1", "UMAP_2", "Clusters")]
cellAnnot_atac <- merge(cellAnnot_atac, df, by = "cellId_archr")

# cluster annotation
cluster <- dplyr::select(ca_ha, Clusters_0.8)
cluster$cellId_archr <- rownames(cluster)
cellAnnot_atac <- merge(cellAnnot_atac, cluster, by = "cellId_archr")
cellAnnot_atac$Cluster <- cellAnnot_atac$Clusters_0.8

# cell type annotation
cell <- dplyr::select(ca_ha, ClusterCellTypes)
cell$cellId_archr <- rownames(cell)
cellAnnot_atac <- merge(cellAnnot_atac, cell, by = "cellId_archr")

# save sample annotation
if (!file.exists(paste0(methDir, "rawData/cellAnnot_atac.tsv"))) {
  write.table(cellAnnot_atac,
    file = paste0(methDir, "rawData/cellAnnot_atac.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}
rownames(cellAnnot_atac) <- cellAnnot_atac$cellId_archr
cellIds_atac <- rownames(cellAnnot_atac)
for (cn in c(grep("^sample", colnames(ca_global), value = TRUE))) {
  cellAnnot_atac[, cn] <- ca_global[ca_global$cellId_archr %in% cellIds_atac, cn]
}
logger.completed()

logger.start("Loading methylation cell annotation")
# cellAnnot_meth <- as.data.frame(data.table::fread("/icbb/projects/igunduz/DARPA/allc_sample_annot_final.csv"))
cellAnnot_meth <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/cellAnnot_meth.rds")
rownames(cellAnnot_meth) <- cellAnnot_meth[, "cellId"]

colname_samplId_atac <- "sample_sampleId_cminid"
colname_samplId_meth <- "CommonMinID"

logger.start("Finding common samples")
sampleIds_atac <- sort(unique(cellAnnot_atac[, colname_samplId_atac]))
sampleIds_meth <- sort(unique(cellAnnot_meth[, colname_samplId_meth]))
summarizeSetOverlap(sampleIds_atac, sampleIds_meth, set1name = "ATAC", set2name = "BS", doVenn = FALSE)
sampleIds_int <- intersect(sampleIds_atac, sampleIds_meth)

nCells_perSample_atac <- sort(table(cellAnnot_atac[cellAnnot_atac[, colname_samplId_atac] %in% sampleIds_int, colname_samplId_atac]))
nCells_perSample_meth <- sort(table(cellAnnot_meth[cellAnnot_meth[, colname_samplId_meth] %in% sampleIds_int, colname_samplId_meth]))
# exclude samples with too few cells
idx <- nCells_perSample_atac[sampleIds_int] >= 200 & nCells_perSample_meth[sampleIds_int] >= 50
logger.info(paste0("excluding ", sum(!idx), " of ", length(idx), " (", round(100 * sum(!idx) / length(idx), 2), "%) samples because they do not contain enough cells"))
sampleIds_int <- sampleIds_int[idx]

sannot_int <- sampleAnnotAll[sampleIds_int, ]
idx <- cellAnnot_atac[, colname_samplId_atac] %in% sampleIds_int
cellIds_atac <- rownames(cellAnnot_atac)[idx]
logger.info(paste0("Retained ", sum(idx), " of ", nrow(cellAnnot_atac), " (", round(100 * sum(idx) / nrow(cellAnnot_atac), 2), "%) ATAC cells for shared samples"))
idx <- cellAnnot_meth[, colname_samplId_meth] %in% sampleIds_int
cellIds_meth <- rownames(cellAnnot_meth)[idx]
logger.info(paste0("Retained ", sum(idx), " of ", nrow(cellAnnot_meth), " (", round(100 * sum(idx) / nrow(cellAnnot_meth), 2), "%) METH cells for shared samples"))

cellAnnot_atac <- cellAnnot_atac[cellIds_atac, ]
cellAnnot_meth <- cellAnnot_meth[cellIds_meth, ]
logger.completed()


logger.start("Preparing methylation data")
se <- readRDS(paste0(methDir, "rawData/methSe_filtered.rds"))
rt <- "archr_peaks"
logger.info(c("Region type:", rt))
methSeL <- assay(se, "mc")
logger.info("Converting to matrix...")
cm_meth <- as(methSeL, "matrix")
logger.completed()

logger.start("Organizing methylation data")
# cm_meth <- methSeL
logger.info("Organizing methylation regions")
rr <- readRDS(paste0("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/regionsGR_filtered.rds"))
rr <- data.table::as.data.table(rr)
rownames(cm_meth) <- paste0("peak", "_", rr$seqnames, "_", rr$start, "to", rr$end)
logger.completed()

################################################################################
# ATAC - methylation cell matching
# using TF-IDF on binary matrices
################################################################################

oDir <- "/icbb/projects/igunduz/DARPA_analysis/artemis_031023"
if (!dir.exists(oDir)) dir.create(oDir)
oDir <- file.path(oDir, "cell_matching")
if (!dir.exists(oDir)) dir.create(oDir)

# tfIdf function
tfIdf <- function(X) {
  tf <- t(t(X) / colSums(X)) # term frequency
  idf <- tf * log(1 + ncol(X) / rowSums(X)) # inverse document frequency
  na <- sum(is.na(x = idf))
  logger.info(paste0("Number of NaNs in idf:", na))
  if (na > 0) {
    idf[is.na(idf)] <- 0
  }
  return(idf)
}
logger.start("Preparing matrices")
regGr <- ArchR::getPeakSet(project)
names(regGr) <- paste0("peak", "_", seqnames(regGr), "_", start(regGr), "to", end(regGr))
# process per chromosome, otherwise the matrices will get too big (cholmod error for sparse matrices)
chroms <- as.character(sort(unique(seqnames(regGr))))
cmL <- lapply(chroms, FUN = function(cc) {
  sec <- getMatrixFromProject(project, "PeakMatrix", useSeqnames = cc, logFile = tempfile())
  sec <- sec[, cellIds_atac]
  rr <- rowRanges(sec)
  names(rowRanges(sec)) <- paste0("peak", "_", seqnames(rr), "_", start(rr), "to", end(rr))
  assay(sec, "PeakMatrix")
})
cm_atac <- do.call(rbind, cmL)
cm_atac <- cm_atac[names(regGr), ]
rm(cmL)



# binarize: is accessible / is unmethylated
logger.start("Binarizing matrices")
cm_atac <- cm_atac > 0
cm_meth <- !is.na(cm_meth) & cm_meth < 0.4
logger.completed()

# logger.start("Filtering low quality regions")
# filter regions with low quality
# keep <- rowSums(cm_meth) > 0
# logger.info(paste0("Retained ", sum(keep), " of ", nrow(cm_meth), " (", round(100 * sum(keep) / nrow(cm_meth), 2), "%) regions"))
# cm_meth <- cm_meth[keep, ]
# logger.completed()

# use the region set of peaks from ATAC intersecting with methylation measurements
logger.info("Region set overlap between ATAC and METH:")
summarizeSetOverlap(rownames(cm_atac), rownames(cm_meth), set1name = "ATAC", set2name = "METH", doVenn = FALSE)
regIds <- intersect(rownames(cm_atac), rownames(cm_meth))
cm_atac <- cm_atac[regIds, ]
cm_meth <- cm_meth[regIds, ]

logger.start("Calculating TF-IDF on scMETH")
tfidf_meth <- tfIdf(cm_meth)
logger.completed()
rm(cm_meth)
# logger.info("Converting to matrix...")
# tfidf_meth <- as(tfidf_meth, "matrix")
ChrAccR::cleanMem()
logger.completed()

seu_meth <- CreateSeuratObject(counts = tfidf_meth, project = "BS_seq", assay = "Bisulfite", meta.data = cellAnnot_meth)
logger.info("Saving Seurat object for meth")
saveRDS(seu_meth, file.path(oDir, "seuratDs_tfidf_meth.rds"))
ChrAccR::cleanMem()

logger.start("Calculating TF-IDF on scATAC")
tfidf_atac <- tfIdf(cm_atac)
logger.completed()
rm(cm_atac)

seu_atac <- CreateSeuratObject(counts = tfidf_atac, project = "ATAC_seq", assay = "ATAC", meta.data = cellAnnot_atac)
logger.info("Saving Seurat object for atac")
saveRDS(seu_atac, file.path(oDir, "seuratDs_tfidf_acc.rds"))
ChrAccR::cleanMem()


sampleIds_int <- unique(seu_atac@meta.data[, colname_samplId_atac])
ChrAccR::cleanMem()
logger.start("CCA matching")
fn_anch <- file.path(oDir, paste0("ccaCellMatches_", "anchors", ".rds"))
fn_nns <- file.path(oDir, paste0("ccaCellMatches_", "nns", ".rds"))

# Define the anchor and nearest neighbor result directories
anchorsDir <- file.path(oDir, "anchors")
nnsDir <- file.path(oDir, "nns")

# Create the anchor and nearest neighbor result directories if they don't exist
if (!file.exists(anchorsDir)) {
  dir.create(anchorsDir)
}
if (!file.exists(nnsDir)) {
  dir.create(nnsDir)
}

matchResL <- lapply(sampleIds_int, FUN = function(sid) {
  anchorFile <- file.path(anchorsDir, paste0("anchors_", sid, ".rds"))
  nnsFile <- file.path(nnsDir, paste0("nns_", sid, ".rds"))

  if (file.exists(anchorFile) && file.exists(nnsFile)) {
    logger.info(paste("Loading existing anchor and nearest neighbor results for", sid))
    anchRes <- readRDS(anchorFile)
    ccaNnRes <- readRDS(nnsFile)
  } else {
    logger.start(c("Processing sample:", sid))
    sa <- seu_atac[, seu_atac@meta.data[, colname_samplId_atac] == sid]
    sm <- seu_meth[, seu_meth@meta.data[, colname_samplId_meth] == sid]
    anch <- FindTransferAnchors(
      reference = sa, query = sm,
      features = rownames(seu_atac),
      reference.assay = "ATAC",
      query.assay = "Bisulfite",
      reduction = "cca",
      k.anchor = 10
    )


    anchRes <- data.frame(
      cell_acc = colnames(sa)[as.integer(anch@anchors[, "cell1"])],
      cell_meth = colnames(sm)[as.integer(anch@anchors[, "cell2"])],
      score = anch@anchors[, "score"],
      sampleId = sid,
      stringsAsFactors = FALSE
    )

    # Nearest neighbors in CCA embeddings
    seu <- anch@object.list[[1]]
    layer <- gsub("^(.+_)(reference|query)$", "\\2", Cells(seu))
    seu@meta.data[, "layer"] <- ifelse(layer == "reference", "accessibility", ifelse(layer == "query", "methylation", NA))
    dimRedM <- Embeddings(seu, reduction = "cca.l2")
    idx_acc <- which(seu@meta.data[rownames(dimRedM), "layer"] == "accessibility")
    idx_meth <- which(seu@meta.data[rownames(dimRedM), "layer"] == "methylation")
    drm_acc <- dimRedM[idx_acc, ]
    rownames(drm_acc) <- gsub("_reference$", "", rownames(drm_acc))
    drm_meth <- dimRedM[idx_meth, ]
    rownames(drm_meth) <- gsub("_query$", "", rownames(drm_meth))
    knnRes_expr <- FNN::get.knnx(drm_meth, drm_acc, k = 1)
    knnRes_acc <- FNN::get.knnx(drm_acc, drm_meth, k = 1)
    ccaNnRes <- data.frame(
      cell_acc = c(rownames(drm_acc), rownames(drm_acc)[knnRes_acc$nn.index[, 1]]),
      cell_meth = c(rownames(drm_meth)[knnRes_expr$nn.index[, 1]], rownames(drm_meth)),
      dist = c(knnRes_expr$nn.dist[, 1], knnRes_acc$nn.dist[, 1]),
      sampleId = sid,
      mapping = rep(c("nearest_meth_for_acc", "nearest_acc_for_meth"), times = c(length(idx_acc), length(idx_meth))),
      stringsAsFactors = FALSE
    )
    if (!file.exists(file.path(oDir, paste0("ccaSeuDataset_", sid, ".rds")))) {
      logger.status("Saving CCA dataset")
      saveRDS(seu, file.path(oDir, paste0("ccaSeuDataset_", sid, ".rds")))
      logger.completed()
    }
    # Save anchRes and ccaNnRes
    logger.start(paste("Saving anchor and nearest neighbor results for", sid))
    saveRDS(anchRes, anchorFile)
    saveRDS(ccaNnRes, nnsFile)
    logger.completed()
    ChrAccR::cleanMem()
  }
  # Return the lists
  return(list(anchors = anchRes, nns = ccaNnRes))
})


matchRes_anchors <- do.call("rbind", lapply(matchResL, FUN = function(x) {
  x[["anchors"]]
}))
matchRes_nns <- do.call("rbind", lapply(matchResL, FUN = function(x) {
  x[["nns"]]
}))

logger.status("Saving ...")
saveRDS(matchRes_anchors, fn_anch)
saveRDS(matchRes_nns, fn_nns)
logger.completed()

#####################################################################
