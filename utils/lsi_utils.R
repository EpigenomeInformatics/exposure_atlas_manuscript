itLsi <- function(methSe, # methylation SummarizedExperiment
                  doImpute = TRUE,
                  doTfIdf = FALSE,
                  it0mostVarRegs = 50000L,
                  it1pcs = 1:30,
                  it1clusterResolution = 0.8,
                  it1mostVarRegs = 50000L,
                  it2pcs = 1:30,
                  it2clusterResolution = 0.8,
                  normPcs = FALSE,
                  umapParams = list(
                    distMethod = "euclidean",
                    min_dist = 0.5,
                    n_neighbors = 25
                  )) {
  cellIds <- colnames(methSe)
  regCoord <- SummarizedExperiment::rowData(methSe)
  dummyMat <- matrix(11.0, ncol = length(cellIds), nrow = 11)
  colnames(dummyMat) <- cellIds
  rownames(dummyMat) <- paste0("df", 1:nrow(dummyMat))
  logger.start("Retrieving methylation matrix")
  # X <- rr
  # X <- Xall <- getSparseMethLvls(methSe, asDense=TRUE)
  X <- Xall <- as.matrix(assay(methSe, "mc"))
  logger.completed()

  imputeMode <- "default"
  if (!is.logical(doImpute)) {
    if (doImpute == "pcimpute") {
      imputeMode <- "pcimpute"
      doImpute <- FALSE
    } else {
      logger.error("Invalid imputation method")
    }
  }

  # logger.start("Checking duplicates...")
  # unique_indices <- !kit::fduplicated(X)
  # logger.info(c("Removing", sum(!unique_indices), "duplicated regions"))
  # X <- Xall <- X[unique_indices,]
  # regCoord <- regCoord[unique_indices, ]
  # logger.completed()

  if (!is.null(it0mostVarRegs) && it0mostVarRegs < nrow(X)) {
    logger.start(c("Identifying variable regions"))
    regVar <- matrixStats::rowVars(X, na.rm = TRUE)
    if (it0mostVarRegs < length(regVar)) {
      idx2rem <- rank(-regVar, na.last = "keep", ties.method = "min") > it0mostVarRegs
      idx2rem <- idx2rem[complete.cases(idx2rem)] # added this due to error
      logger.info(c("Retaining the", sum(!idx2rem), "most variable regions"))
      X <- X[!idx2rem, ]
      regCoord <- regCoord[!idx2rem, ] # filter out the regions as well
    }
    logger.completed()
  }
  if (doImpute && any(is.na(X))) {
    logger.start("Imputing missing values using 5-nearest-neighbor imputation")
    X <- handleMissingData(X, method = "impute.knn5")$data
    logger.completed()
  }

  logger.start("Iteration 1")
  if (doTfIdf) {
    # binarize methylation calls
    X <- X < 0.4 # is unmethylated
    logger.start("TF-IDF normalization")
    X <- tfIdf(X)
    logger.completed()
  }
  if (imputeMode == "pcimpute") {
    pcaCoord_it1 <- muRtools::getDimRedCoords.pca(t(X), components = 1:max(it1pcs), method = "prcomp_iter_impute_scbs")
  } else {
    pcaCoord_it1 <- muRtools::getDimRedCoords.pca(t(X), components = 1:max(it1pcs), method = "irlba_svd")
  }
  pcM <- pcaCoord_it1
  if (normPcs) {
    logger.info("Scaling SVDs")
    pcM <- ChrAccR::rowZscores(pcM, na.rm = TRUE)
  }
  pcM <- pcM[, it1pcs, drop = FALSE]
  logger.start(c("Clustering"))
  sObj <- Seurat::CreateSeuratObject(dummyMat, project = "scMeth", min.cells = 0, min.features = 0, assay = "meth")
  sObj[["pca"]] <- Seurat::CreateDimReducObject(embeddings = pcM, key = "PC_", assay = "meth")
  sObj <- Seurat::FindNeighbors(sObj, reduction = "pca", assay = "meth", dims = 1:ncol(pcM), k.param = 30)
  clustRes <- Seurat::FindClusters(sObj, algorithm = 1, n.start = 100, n.iter = 10, resolution = it1clusterResolution, verbose = FALSE)
  clustAss_it1 <- factor(paste0("c", clustRes@active.ident), levels = paste0("c", levels(clustRes@active.ident)))
  names(clustAss_it1) <- names(clustRes@active.ident)
  logger.info(c("Number of clusters found:", nlevels(clustAss_it1)))
  logger.completed()

  X <- Xall
  rGr <- regCoord

  if (!is.null(it1mostVarRegs) && it1mostVarRegs < nrow(X)) {
    logger.start(c("Identifying cluster-variable regions"))
    logger.start("Creating cluster pseudo-bulk samples")
    Xc <- do.call(cbind, tapply(1:ncol(X), clustAss_it1, FUN = function(iis) {
      matrixStats::rowMeans2(X[, iis], na.rm = TRUE)
    }))
    logger.completed()
    logger.start("Identifying cluster variable regions")
    regVar <- matrixStats::rowVars(Xc, na.rm = TRUE)
    if (it1mostVarRegs < length(regVar)) {
      idx2rem <- rank(-regVar, na.last = "keep", ties.method = "min") > it1mostVarRegs
      idx2rem <- idx2rem[complete.cases(idx2rem)] # added this due to error
      logger.info(c("Retaining the", sum(!idx2rem), "most variable regions across clusters"))
      Xc <- Xc[!idx2rem, ] # filter Xc instead X
      X <- X[rownames(X) %in% rownames(Xc), ] # then filter X
    }
    logger.completed()
    logger.completed()
  }
  logger.completed()
  logger.start("Iteration 2")
  if (doImpute && any(is.na(X))) {
    logger.start("Imputing missing values using 5-nearest-neighbor imputation")
    X <- handleMissingData(X, method = "impute.knn5")$data
    logger.completed()
  }
  logger.start(c("Performing dimension reduction"))
  idfBase <- NULL
  if (doTfIdf) {
    # binarize methylation calls
    X <- X < 0.4 # is unmethylated
    idfBase <- log(1 + ncol(X) / rowSums(X, na.rm = TRUE))
    logger.start("TF-IDF normalization")
    X <- tfIdf(X)
    logger.completed()
  }
  if (imputeMode == "pcimpute") {
    pcaCoord_it2 <- muRtools::getDimRedCoords.pca(t(X), components = 1:max(it2pcs), method = "prcomp_iter_impute_scbs")
  } else {
    pcaCoord_it2 <- muRtools::getDimRedCoords.pca(t(X), components = 1:max(it2pcs), method = "irlba_svd")
  }
  pcM <- pcaCoord_it2
  if (normPcs) {
    logger.info("Scaling SVDs")
    pcM <- ChrAccR::rowZscores(pcM, na.rm = TRUE)
  }
  pcM <- pcM[, it2pcs, drop = FALSE]

  paramL <- c(list(X = pcM), umapParams)
  umapCoord <- do.call(muRtools::getDimRedCoords.umap, paramL)
  umapRes <- attr(umapCoord, "umapRes")
  attr(umapCoord, "umapRes") <- NULL
  logger.completed()

  logger.start(c("Clustering"))
  sObj <- Seurat::CreateSeuratObject(dummyMat, project = "scMeth", min.cells = 0, min.features = 0, assay = "meth")
  sObj[["pca"]] <- Seurat::CreateDimReducObject(embeddings = pcM, key = "PC_", assay = "meth")
  sObj <- Seurat::FindNeighbors(sObj, reduction = "pca", assay = "meth", dims = 1:ncol(pcM), k.param = 30)
  clustRes <- Seurat::FindClusters(sObj, algorithm = 1, n.start = 100, n.iter = 10, resolution = it2clusterResolution, verbose = FALSE)
  clustAss <- factor(paste0("c", clustRes@active.ident), levels = paste0("c", levels(clustRes@active.ident)))
  names(clustAss) <- names(clustRes@active.ident)
  logger.info(c("Number of clusters found:", nlevels(clustAss)))
  logger.completed()
  logger.completed()
  res <- list(
    pcaCoord = pcaCoord_it2,
    pcs = it2pcs,
    idfBase = idfBase,
    umapCoord = umapCoord,
    umapRes = umapRes,
    clustAss = clustAss,
    regionGr = rGr,
    regionGr_unfiltered = regCoord,
    iterationData = list(
      iteration0 = list(
        mostVarRegs = it0mostVarRegs
      ),
      iteration1 = list(
        pcaCoord = pcaCoord_it1,
        clustAss = clustAss_it1,
        pcs = it1pcs,
        clusterResolution = it1clusterResolution,
        mostVarRegs = it1mostVarRegs
      )
    ),
    .params = list(
      doTfIdf = doTfIdf,
      normPcs = normPcs,
      umapParams = umapParams
    )
  )
  class(res) <- "iterativeLSIResultScMeth"
  return(res)
}
