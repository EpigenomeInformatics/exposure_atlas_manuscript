jaccard_helper <- function(mat, i, j) {
  AiB <- mat[i, j]
  AuB <- sum(mat[i, ]) + sum(mat[, j]) - AiB
  AiB / AuB
}


Jaccard <- function(data) {
  data <- t(data)
  df <- data # to assign names to later
  jaccards <- matrix(data = NA, nrow = NROW(data), ncol = ncol(data))
  for (r in 1:NROW(data)) {
    for (c in 1:ncol(data)) {
      jaccards[r, c] <- jaccard_helper(data, r, c)
    }
  }

  colnames(jaccards) <- colnames(data)
  row.names(jaccards) <- row.names(data)
  pheatmap <- pheatmap::pheatmap(jaccards,
    show_rownames = T,
    cluster_rows = F,
    display_numbers = as.matrix(df)
  )
  return(pheatmap)
}

prettyOrderMat <- function(mat, scale = TRUE, cutOff = 1, lmat = NULL, clusterCols = TRUE) {
  # Reorder mat in a prettier way for plotting
  # Adapted from Jeff's ArchR .binarySort
  ###################################
  # mat = matrix (like) object to sort
  # scale = should mat be scaled before building logical mat
  # cutOff = cutoff for lmat
  # lmat = logical matrix for ordering rows (binary sorting)
  # clusterCols = should columns be clustered?
  mat <- as.matrix(mat)

  if (is.null(lmat)) {
    # Compute row Z-scores
    if (scale) {
      lmat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), `/`)
    } else {
      lmat <- mat
    }
    # Logical matrix of values passing cutoff
    lmat <- lmat >= cutOff
  }

  # Transpose:
  mat <- t(mat)
  lmat <- t(lmat)

  # Identify column ordering:
  if (clusterCols) {
    hc <- hclust(dist(mat))
    colIdx <- hc$order
    mat <- t(mat[colIdx, ])
    lmat <- t(lmat[colIdx, ])
  } else {
    mat <- t(mat)
    lmat <- t(lmat)
    hc <- NULL
  }

  # Identify row ordering:
  rowIdx <- do.call("order", c(as.data.frame(lmat)[seq_len(ncol(lmat))], list(decreasing = TRUE)))
  mat <- mat[rowIdx, ]

  return(list(mat = mat, hclust = hc))
}


# Heatmap wrapper:
BORHeatmap <- function(
    mat, # Data to plot (matrix or dataframe)
    limits = NULL, # Enforced limits for colormap (2 dimensional array)
    clusterCols = TRUE, # Should columns be clustered
    clusterRows = TRUE, # Should rows be clustered
    labelCols = FALSE, # Should columns be labeled
    labelRows = FALSE, # Should rows be labeled
    dataColors = NULL, # Colormap for plotting data
    dataColorMidPoint = NULL, # The data value to be the middle of the color map
    customRowLabel = NULL,
    customRowLabelIDs = NULL,
    customColLabel = NULL,
    customColLabelIDs = NULL,
    customLabelWidth = 0.15,
    useRaster = TRUE, # Should heatmap be rasterized
    rasterDevice = "CairoPNG",
    rasterQuality = 5, # Raster quality. Higher is {better?}
    fontSize = 6, # Font size for labels
    showColDendrogram = FALSE, # Should the column dendrogram be shown
    showRowDendrogram = FALSE, # Should the row dendrogram be shown
    borderColor = NA, # Color for lines between cells
    mapname = " ", # 'Name' to give heatmap
    legendTitle = " ", # Name of legend
    ...) {
  # Packages
  suppressPackageStartupMessages(require(ComplexHeatmap))
  suppressPackageStartupMessages(require(circlize))

  # Make sure mat is actually a matrix
  if (!is.matrix(mat)) {
    message("'mat' needs to be a matrix. Converting...")
    mat <- as.matrix(mat)
  }

  # Prepare color function
  if (!is.null(limits)) {
    ll <- limits[1]
    ul <- limits[2]
  } else {
    ll <- min(mat, na.rm = TRUE)
    ul <- max(mat, na.rm = TRUE)
  }
  # If no colormap provided, use solarExtra
  if (is.null(dataColors)) {
    dataColors <- c(
      "1" = "#3361A5", "2" = "#248AF3", "3" = "#14B3FF",
      "4" = "#88CEEF", "5" = "#C1D5DC", "6" = "#EAD397",
      "7" = "#FDB31A", "8" = "#E42A2A", "9" = "#A31D1D"
    )
  }
  dataColFun <- makeColFun(ll, ul, dataColors, midpoint = dataColorMidPoint)

  message("Preparing Heatmap...")
  hm <- Heatmap(
    # Main components:
    matrix = mat,
    name = mapname,
    col = dataColFun,

    # Legend options:
    heatmap_legend_param = list(
      color_bar = "continuous",
      legend_direction = "vertical",
      legend_width = unit(1, "cm"),
      title = legendTitle
    ),
    rect_gp = gpar(col = borderColor),

    # Column options:
    show_column_names = labelCols,
    cluster_columns = clusterCols,
    show_column_dend = showColDendrogram,
    clustering_method_columns = "ward.D2",
    # column_names_gp = gpar(fontsize = fontSize),

    # Row options:
    show_row_names = labelRows,
    cluster_rows = clusterRows,
    show_row_dend = showRowDendrogram,
    clustering_method_rows = "ward.D2",
    # row_names_gp = gpar(fontsize = fontSize),

    # Raster info:
    use_raster = useRaster,
    raster_device = rasterDevice,
    raster_quality = rasterQuality,
  )

  # Add row labels if provided:
  if (!is.null(customRowLabel)) {
    if (is.null(customRowLabelIDs)) {
      customRowLabelIDs <- rownames(mat)[customRowLabel]
    }
    hm <- hm + rowAnnotation(
      link = anno_mark(at = customRowLabel, labels = customRowLabelIDs, labels_gp = gpar(fontsize = fontSize)),
      width = unit(customLabelWidth, "cm") + max_text_width(customRowLabelIDs)
    )
  }

  return(hm)
}

# This is used primarily for making colormaps for ComplexHeatmap
makeColFun <- function(start, end, cmap, midpoint = NULL) {
  # Make a color ramp function from provided start and end breaks,
  # and optionally a midpoint
  cmapLen <- length(cmap)
  if (!is.null(midpoint)) {
    interpolate <- function(c1, c2, colorspace = "Lab") {
      rgb2hex(colorRamp(c(c1, c2), space = colorspace)(0.5))
    }
    if (length(cmap) %% 2 == 0) {
      # Interpolate middle colors if necessary to get midpoint
      preMidIdx <- floor(cmapLen / 2)
      midCol <- interpolate(cmap[preMidIdx], cmap[preMidIdx + 1])
      cmap <- c(cmap[1:preMidIdx], midCol, cmap[(preMidIdx + 1):cmapLen])
      cmapLen <- length(cmap)
    }
    midIdx <- ceiling(cmapLen / 2)
    breaks <- c(seq(start, midpoint, length.out = midIdx), seq(midpoint, end, length.out = midIdx)[2:midIdx])
  } else {
    breaks <- seq(start, end, length.out = cmapLen)
  }
  colorRamp2(breaks, cmap)
}

# Compute the z-score, helper for clusterATACzscores
computeZScore <- function(counts) {
  counts <- (counts - matrixStats::rowMeans2(counts)) / matrixStats::rowSds(counts)
  counts[base::is.nan(counts)] <- 0
  return(counts)
}

# helper function for atac
cutL0.5fc2Padj05 <- function(dm, padj = 0.05) {
  abs(dm[, "log2FoldChange"]) > 0.5 & dm[, "padj"] < padj
}

# Consensus K-means clustering
consensus_kmeans <- function(counts, k, km_repeats) {
  partition_list <- lapply(seq_len(km_repeats), function(i) {
    as.cl_hard_partition(kmeans(counts, k, iter.max = 50))
  })
  partition_list <- cl_ensemble(list = partition_list)
  partition_consensus <- cl_consensus(partition_list)
  as.vector(cl_class_ids(partition_consensus))
}

# Cutoff function for differential peaks
cutL2FCpadj <- function(data, lfc = 0.5, padj = 0.05, logpadj = FALSE) {
  if (logpadj) {
    if (any(c("log2FoldChange", "logPval")) %in% colnames(data)) {
      stop("Please provide a data.table with log2FoldChange and logPval columns!")
    }
    return(abs(data[, "log2FoldChange"]) > lfc & data[, "logPval"] < padj)
  } else {
    if (any(c("log2FoldChange", "padj")) %in% colnames(data)) {
      stop("Please provide a data.table with log2FoldChange and padj columns!")
    }
    return(abs(data[, "log2FoldChange"]) > lfc & data[, "padj"] < padj)
  }
}
