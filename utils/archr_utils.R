cellTypeProportionPlot <- function(
    mat, scale = TRUE, center = FALSE, groupName = NULL,
    colorPalette = NULL, theme = theme_minimal(), order = NULL) {
  if (is.null(groupName)) {
    stop("Please provide a group column name!")
  }
  if (!class(mat) %in% c("matrix", "table")) {
    stop("Please provide a matrix!")
  }
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
    !requireNamespace("forcats", quietly = TRUE) ||
    !requireNamespace("reshape2", quietly = TRUE)) {
    stop("Please make sure ggplot2, forcats, and reshape2 packages are installed.")
  }
  if (scale) {
    mat <- scale(mat, center = center, scale = colSums(mat))
  }
  mat <- reshape2::melt(mat)
  colnames(mat) <- c("CellTypes", groupName, "FractionOfCells")
  mat$Group <- forcats::fct_reorder(mat[[groupName]], mat$FractionOfCells, mean, .desc = F)
  if (!is.null(order)) {
    mat$Group <- factor(mat$Group, levels = order, ordered = FALSE)
  }
  stacked <- ggplot(mat, aes(fill = CellTypes, y = FractionOfCells, x = Exposure)) +
    geom_bar(position = "fill", stat = "identity") +
    labs(x = groupName, y = "Cell Type Proportions") +
    theme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  if (!is.null(colorPalette)) {
    stacked <- stacked + scale_fill_manual(values = colorPalette)
  }
  return(stacked)
}


plotUniqueFragsvsTSS <- function(project) {
  df <- getCellColData(project, select = c("log10(nFrags)", "TSSEnrichment"))
  # plot QC
  p <- ggPoint(
    x = df[, 1],
    y = df[, 2],
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10 Unique Fragments",
    ylabel = "TSS Enrichment",
    xlim = c(log10(500), quantile(df[, 1], probs = 0.99)),
    ylim = c(0, quantile(df[, 2], probs = 0.99))
  ) + geom_hline(yintercept = 4, lty = "dashed") + geom_vline(xintercept = 3, lty = "dashed")
  return(p)
}

plotNiceArchRumap <- function(df, colorPalette, legendpos = "bottom") {
  if (!require(ggplot2)) {
    stop("ggplot2 is not installed")
  }
  if (class(df) != "data.frame") {
    stop("df must be a data frame")
  }
  if (!all(c("UMAP1", "UMAP2", "Group") %in% colnames(df))) {
    stop("df must have columns named UMAP1, UMAP2, and Group")
  }
  if (length(unique(df$Group)) > length(colorPalette)) {
    stop("colorPalette must have at least as many colors as there are unique values in df$Group")
  }

  umap <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = Group)) +
    geom_point() +
    scale_color_manual(values = colorPalette) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme_void() +
    theme(legend.position = legendpos)
  return(umap)
}

# These scripts are used to run the ArchR pipeline.
# They are mainly taken from ArchR's source code, with some modifications to make them work with exposure atlas data, or speeding up.

addPeakAnnotationsNew <- function(
    ArchRProj = NULL,
    regions = NULL,
    name = "Region",
    force = FALSE,
    logFile = createLogFile("addPeakAnnotations")) {
  ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  ArchR:::.validInput(input = regions, name = "regions", valid = c("grangeslist", "character"))
  ArchR:::.validInput(input = name, name = "name", valid = c("character"))
  ArchR:::.validInput(input = force, name = "force", valid = c("boolean"))
  ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))

  tstart <- Sys.time()
  ArchR:::.startLogging(logFile = logFile)
  ArchR:::.logThis(mget(names(formals()), sys.frame(sys.nframe())), "addPeakAnnotations Input-Parameters", logFile = logFile)

  if (name %in% names(ArchRProj@peakAnnotation)) {
    if (force) {
      ArchR:::.logMessage("peakAnnotation name already exists! Overriding.", verbose = TRUE, logFile = logFile)
    } else {
      ArchR:::.logStop("peakAnnotation name already exists! set force = TRUE to override!", logFile = logFile)
    }
  }

  if (inherits(regions, "GRanges")) {
    regionPositions <- GRangesList(region = regions)
  } else {
    if (is.null(names(regions))) {
      names(regions) <- paste0("Region_", seq_along(regions))
    }

    if (any(duplicated(names(regions)))) {
      stop("Found duplicated region names! Please make unique!")
    }

    regionPositions <- lapply(seq_along(regions), function(x) {
      if (inherits(regions[[x]], "GRanges")) {
        gr <- ArchR:::.validGRanges(regions[[x]])
      } else if (is.character(regions[[x]])) {
        gr <- tryCatch(
          {
            makeGRangesFromDataFrame(
              df = data.frame(data.table::fread(regions[[x]])),
              keep.extra.columns = FALSE,
              seqnames.field = "V1",
              start.field = "V2",
              end.field = "V3"
            )
          },
          error = function(y) {
            ArchR:::.logMessage(paste0("Could not successfully get region : ", regions[[x]]), verbose = TRUE, logFile = logFile)

            if (!file.exists(regions[[x]])) {
              ArchR:::.logStop(paste0("If region provided is a path it does not exist!"), logFile = logFile)
            }

            ArchR:::.logStop("Could not create GRanges from region", logFile = logFile)
          }
        )
      } else {
        .logStop("Unrecognized input in regions please input GRanges, GRangesList, or Paths to bed files!", logFile = logFile)
      }

      gr
    }) %>% GRangesList()

    names(regionPositions) <- names(regions)
  }

  #############################################################
  # Peak Overlap Matrix
  #############################################################
  peakSet <- getPeakSet(ArchRProj)
  if (is.null(peakSet)) {
    .logStop("peakSet is NULL. You need a peakset to run addMotifAnnotations! See addReproduciblePeakSet!", logFile = logFile)
  }
  allPositions <- unlist(regionPositions, use.names = FALSE)

  ArchR:::.logDiffTime("Creating Peak Overlap Matrix", t1 = tstart, verbose = TRUE, logFile = logFile)

  overlapRegions <- findOverlaps(peakSet, allPositions, ignore.strand = TRUE)
  if (length(overlapRegions) == 0) {
    stop("No Overlaps Found between regions and peak Matrix!")
  }
  ArchR:::.logThis(overlapRegions, "overlapRegions", logFile = logFile)

  regionMat <- Matrix::sparseMatrix(
    i = queryHits(overlapRegions),
    j = match(names(allPositions), names(regionPositions))[subjectHits(overlapRegions)],
    x = rep(TRUE, length(overlapRegions)),
    dims = c(length(peakSet), length(regionPositions))
  )
  colnames(regionMat) <- names(regionPositions)
  ArchR:::.logThis(regionMat, "regionMat", logFile = logFile)

  regionMat <- SummarizedExperiment::SummarizedExperiment(assays = SimpleList(matches = regionMat), rowRanges = peakSet)
  ArchR:::.logThis(regionMat, "regionSE", logFile = logFile)

  #############################################################
  # Filter Regions With No Matches
  #############################################################

  # Number of Overlaps
  nO <- Matrix::colSums(assay(regionMat))
  rF <- names(which(nO == 0))

  if (all(nO == 0)) {
    stop("No Overlaps Found! Please check your peakSet and genome!")
  }

  if (length(rF) > 0) {
    ArchR:::.logDiffTime(paste0("Filtering Region Annotations with 0 overlaps :\n\n ", paste(rF, collapse = ", "), "\n\n"), t1 = tstart, verbose = TRUE, logFile = logFile)
    # Filter
    regionPositions <- regionPositions[!(names(regionPositions) %in% rF)]
    regionMat <- regionMat[, names(regionPositions), drop = FALSE]
  } else {
    ArchR:::.logDiffTime(paste0("All Regions Overlap at least 1 peak!"), t1 = tstart, verbose = TRUE, logFile = logFile)
  }

  #############################################################
  # Summarize and Save
  #############################################################

  dir.create(file.path(getOutputDirectory(ArchRProj), "Annotations"), showWarnings = FALSE)
  savePositions <- file.path(getOutputDirectory(ArchRProj), "Annotations", paste0(name, "-Positions-In-Peaks.rds"))
  saveMatches <- file.path(getOutputDirectory(ArchRProj), "Annotations", paste0(name, "-Matches-In-Peaks.rds"))

  out <- SimpleList(
    regionMatches = regionMat,
    regionPositions = regionPositions,
    date = Sys.Date()
  )

  ArchRProj@peakAnnotation[[name]]$Name <- name
  ArchRProj@peakAnnotation[[name]]$Positions <- savePositions
  ArchRProj@peakAnnotation[[name]]$Matches <- saveMatches

  ArchR:::.safeSaveRDS(out, file.path(getOutputDirectory(ArchRProj), "Annotations", paste0(name, "-In-Peaks-Summary.rds")), compress = FALSE)
  ArchR:::.safeSaveRDS(out$regionPositions, savePositions, compress = FALSE)
  ArchR:::.safeSaveRDS(out$regionMatches, saveMatches, compress = FALSE)

  return(ArchRProj)
}


projectBulkATAC_vs2 <- function(
    ArchRProj = NULL,
    seATAC = NULL,
    reducedDims = "IterativeLSI",
    embedding = "UMAP",
    n = 250,
    verbose = TRUE,
    threads = getArchRThreads(),
    force = FALSE,
    logFile = createLogFile("projectBulkATAC")) {
  ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  ArchR:::.validInput(input = seATAC, name = "seATAC", valid = c("SummarizedExperiment"))
  ArchR:::.validInput(input = reducedDims, name = "reducedDims", valid = c("character"))
  ArchR:::.validInput(input = embedding, name = "embedding", valid = c("character", "null"))
  ArchR:::.validInput(input = n, name = "n", valid = c("integer"))
  ArchR:::.validInput(input = verbose, name = "verbose", valid = c("boolean"))
  ArchR:::.validInput(input = threads, name = "threads", valid = c("integer"))
  ArchR:::.validInput(input = force, name = "force", valid = c("boolean"))
  ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))

  tstart <- Sys.time()

  ArchR:::.startLogging(logFile = logFile)
  # ArchR:::.logThis(mget(names(formals()),sys.frame(sys.nframe())), "projectBulkATAC Input-Parameters", logFile = logFile)

  ##################################################
  # Reduced Dimensions
  ##################################################
  rD <- getReducedDims(ArchRProj, reducedDims = reducedDims, returnMatrix = FALSE)
  ArchR:::.logThis(names(rD), "reducedDimsNames", logFile = logFile)
  ArchR:::.logThis(rD[[1]], "reducedDimsMat", logFile = logFile)
  rDFeatures <- rD[[grep("Features", names(rD))]]
  if ("end" %in% colnames(rDFeatures)) {
    rDGR <- GRanges(seqnames = rDFeatures$seqnames, IRanges(start = rDFeatures$start, end = rDFeatures$end))
  } else {
    rDGR <- GRanges(seqnames = rDFeatures$seqnames, IRanges(start = rDFeatures$start, width = (rDFeatures$start) / (rDFeatures$idx - 1)))
  }
  ArchR:::.logThis(rDGR, "reducedDimsGRanges", logFile = logFile)
  subATAC <- subsetByOverlaps(seATAC, rDGR, ignore.strand = TRUE)
  subATAC <- subATAC[order(rowSums(as.matrix(ArchR:::.getAssay(subATAC, "counts"))), decreasing = TRUE), ]
  o <- DataFrame(findOverlaps(subATAC, rDGR, ignore.strand = TRUE))
  sumOverlap <- length(unique(o[, 2]))
  ArchR:::.logThis(o, "overlapATAC", logFile = logFile)

  if (sumOverlap == 0) {
    ArchR:::.logStop(paste0(
      "No overlaps between bulk ATAC data and reduce dimensions feature found.",
      "\nEither recreate counts matrix or most likely these data sets are incompatible!"
    ), logFile = logFile)
  }
  if ((sumOverlap / length(rDGR)) < 0.25) {
    if (force) {
      ArchR:::.logMessage("Less than 25% of the features are present in this bulk ATAC data set! Continuing since force = TRUE!", verbose = TRUE, logFile = logFile)
    } else {
      ArchR:::.logStop("Less than 25% of the features are present in this bulk ATAC data set! Set force = TRUE to continue!", logFile = logFile)
    }
  }
  ArchR:::.logMessage("Overlap Ratio of Reduced Dims Features = ", (sumOverlap / length(rDGR)), verbose = TRUE, logFile = logFile)

  o <- o[!duplicated(o$subjectHits), ]
  subATAC <- subATAC[o$queryHits, ]
  rownames(subATAC) <- paste0("f", o$subjectHits)
  ArchR:::.logThis(subATAC, "subsettedATAC", logFile = logFile)

  ##################################################
  # Create Bulk Matrix
  ##################################################
  bulkMat <- ArchR:::.safeSubset(
    mat = ArchR:::.getAssay(subATAC, "counts"),
    subsetRows = paste0("f", seq_along(rDGR))
  )
  ArchR:::.logThis(bulkMat, "bulkATACMat", logFile = logFile)

  ##################################################
  # Simulate and Project
  ##################################################
  depthN <- round(sum(rD$rowSm / rD$nCol))
  nRep <- 5
  n2 <- ceiling(n / nRep)
  ratios <- c(2, 1.5, 1, 0.5, 0.25) # range of ratios of number of fragments

  simRD <- ArchR:::.safelapply(seq_len(ncol(bulkMat)), function(x) {
    ArchR:::.logDiffTime(sprintf("Projecting Sample (%s of %s)", x, ncol(bulkMat)), t1 = tstart, verbose = verbose, logFile = logFile)
    counts <- bulkMat[, x]
    counts <- rep(seq_along(counts), counts)
    simMat <- lapply(seq_len(nRep), function(y) {
      ratio <- ratios[y]
      simMat <- matrix(sample(x = counts, size = ceiling(ratio * depthN) * n2, replace = TRUE), ncol = n2)
      simMat <- Matrix::summary(as(simMat, "dgCMatrix"))[, -1, drop = FALSE]
      simMat[, 1] <- simMat[, 1] + (y - 1) * n2
      simMat
    }) # %>%  Reduce("rbind", .)
    simMat <- rlist::list.rbind(simMat)
    simMat <- Matrix::sparseMatrix(i = simMat[, 2], j = simMat[, 1], x = rep(1, nrow(simMat)), dims = c(length(rDGR), n2 * nRep))
    simRD <- as.matrix(ArchR:::.projectLSI(simMat, LSI = rD, verbose = FALSE))
    rownames(simRD) <- paste0(colnames(bulkMat)[x], "#", seq_len(nrow(simRD)))
    simRD
  }, threads = threads) # %>% Reduce("rbind", .)

  simRD <- rlist::list.rbind(simRD)

  if (is.null(embedding)) {
    if (rD$scaleDims) {
      simRD <- ArchR:::.scaleDims(simRD)
    }
    out <- SimpleList(
      simulatedReducedDims = simRD
    )
    return(out)
  }
  ArchR:::.logThis(simRD, "simulatedReducedDims", logFile = logFile)

  ##################################################
  # Prep Reduced Dims
  ##################################################
  embedding <- getEmbedding(ArchRProj = ArchRProj, embedding = embedding, returnDF = FALSE)
  corCutOff <- embedding$params$corCutOff
  dimsToUse <- embedding$params$dimsToUse
  scaleDims <- embedding$params$scaleDims

  if (is.null(scaleDims)) {
    scaleDims <- rD$scaleDims
  }

  simRD <- ArchR:::.scaleDims(simRD)

  if (embedding$params$nc != ncol(simRD)) {
    if (is.null(dimsToUse)) {
      dimsToUse <- seq_len(ncol(rD[[1]]))
    }

    if (!is.null(corCutOff)) {
      if (scaleDims) {
        corToDepth <- rD$corToDepth$scaled
        dimsToUse <- dimsToUse[corToDepth < corCutOff]
      } else {
        corToDepth <- rD$corToDepth$none
        dimsToUse <- dimsToUse[corToDepth < corCutOff]
      }
    }

    if (embedding$params$nc != ncol(simRD)) {
      ArchR:::.logMessage("Error! Inconsistency found with matching LSI dimensions to those used in addUMAP or addTSNE",
        "\nReturning with simulated reduced dimension coordinates...",
        verbose = TRUE, logFile = logFile
      )
      out <- SimpleList(
        simulatedReducedDims = simRD
      )
      return(out)
    }

    simRD <- simRD[, dimsToUse, drop = FALSE]
  }

  ##################################################
  # Get Previous UMAP Model
  ##################################################
  umapModel <- ArchR:::.loadUWOT(embedding$params$uwotModel, embedding$params$nc)

  idx <- sort(sample(seq_len(nrow(rD[[1]])), min(nrow(rD[[1]]), 5000))) # Try to use 5000 or total cells to check validity
  rD2 <- getReducedDims(
    ArchRProj = ArchRProj,
    reducedDims = reducedDims,
    dimsToUse = embedding$params$dimsToUse,
    scaleDims = embedding$params$scaleDims,
    corCutOff = embedding$params$corCutOff
  )[idx, , drop = FALSE]

  ##################################################
  # Project UMAP
  ##################################################
  set.seed(1)
  threads2 <- max(floor(threads / 2), 1)
  simUMAP <- uwot::umap_transform(
    X = rbind(rD2, simRD),
    model = umapModel,
    verbose = TRUE,
    n_threads = threads2
  )
  rownames(simUMAP) <- c(rownames(rD2), rownames(simRD))
  ArchR:::.logThis(simUMAP, "simulatedUMAP", logFile = logFile)

  # Check if the projection matches using previous data
  c1 <- cor(simUMAP[rownames(rD2), 1], embedding[[1]][rownames(rD2), 1])
  c2 <- cor(simUMAP[rownames(rD2), 2], embedding[[1]][rownames(rD2), 2])
  if (min(c1, c2) < 0.8) {
    ArchR:::.logMessage("Warning projection correlation is less than 0.8 (R = ", round(min(c1, c2), 4), ").\nThese results may not be accurate because of the lack of heterogeneity in the single cell data.", verbose = TRUE, logFile = logFile)
  }

  dfUMAP <- embedding[[1]]
  colnames(dfUMAP) <- c("UMAP1", "UMAP2")
  colnames(simUMAP) <- c("UMAP1", "UMAP2")
  dfUMAP <- DataFrame(dfUMAP)
  dfUMAP$Type <- Rle("scATAC", lengths = nrow(dfUMAP))

  simUMAP <- DataFrame(simUMAP[rownames(simRD), , drop = FALSE])
  simUMAP$Type <- Rle(stringr::str_split(rownames(simUMAP), pattern = "#", simplify = TRUE)[, 1])

  out <- SimpleList(
    simulatedBulkUMAP = simUMAP,
    singleCellUMAP = dfUMAP,
    simulatedReducedDims = simRD
  )
  ArchR:::.endLogging(logFile = logFile)

  return(out)
}
