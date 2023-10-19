converToEPP <- function(granges, save = TRUE, filePath = NULL) {
  if (!class(granges) %in% c("data.frame", "GRanges", "GRangesList", "CompressedGRangesList")) {
    # check the type of the object is valid
    stop("Provided object doesn't contain a GRanges,GRangesList or data.frame!")
  }
  if (!isTRUE(save) || !isFALSE(save)) {
    if (!is.null(filePath) && is.character(filePath)) {
      save <- TRUE
    } else {
      save <- FALSE
      logger:log_info("Invalid save option is provided, setting save to FALSE.")
    }
  }
  if (is.data.frame(granges)) {
    granges <- GenomicRanges::makeGRangesFromDataFrame(granges, keep.extra.columns = TRUE)
  }
  if (class(granges) == "GRanges") {
    tiles <- EPP_helper(granges)
  }
  if (!save) {
    if (class(granges) %in% c("GRangesList", "CompressedGRangesList")) {
      tiles <- lapply(granges, EPP_helper(granges))
    }
    return(tiles)
  } else {
    if (!is.null(filePath) && is.character(filePath)) {
      if (!dir.exists(filePath)) {
        # create dir if not exist already
        dir.create(filePath)
      }
      if (class(granges) %in% c("GRangesList", "CompressedGRangesList")) {
        # loop over for each GRanges within list
        tiles <- lapply(names(granges), function(i) {
          if (!file.exists(paste0(filePath, "/", i, ".tsv"))) {
            logger::log_info("Processing ", i)
            gr <- EPP_helper(granges[[i]])
            ChrAccR:::cleanMem()
            logger::log_success()
            logger::log_info("Writing ", i, " into tsv file")
            data.table::fwrite(data.table::as.data.table(gr), paste0(filePath, "/", i, ".tsv"), TRUE, FALSE,
              sep = "\t", row.names = FALSE, col.names = FALSE, showProgress = FALSE
            )
            logger::log_success()
          }
        })
      } else {
        logger::log_success()
        logger::log_info("Writing EPP into tsv file")
        data.table::fwrite(data.table::as.data.table(tiles), paste0(filePath, "/sample_EPP.tsv"), TRUE, FALSE,
          sep = "\t", row.names = FALSE, col.names = FALSE, showProgress = isFALSE
        )
        logger::log_success()
      }
    }
  }
}

EPP_helper <- function(granges,ignoreStrand=TRUE) {
  tiles <- GenomicRanges::resize(granges, 3)
  #granges <- granges[!duplicated(granges)]
  hits <- IRanges::findOverlaps(tiles, granges,ignore.strand=ignoreStrand)
  agg <- aggregate(granges, hits, score = sum(score), coverage = sum(coverage))
  tiles <- data.table::as.data.table(tiles)
  tiles[, score2 := round((round(agg$score * agg$coverage, 0) / agg$coverage) * 1000, 0)]
  tiles[, score := paste0(round(agg$score * agg$coverage, 0), "/", agg$coverage)]
  tiles[, .(seqnames, start, end, score, score2, strand)]
  tiles <- tiles %>% dplyr::select(seqnames, start, end, score, score2, strand)
  return(tiles)
}

diagDivCellHeatmap <- function(ml, mr, col.l=NULL, col.r=NULL, name.l="Lower left", name.r="Upper right", ...){
	require(ComplexHeatmap)
	if (nrow(ml)!=nrow(mr)) stop("Numbers of rows of the two matrices must match")
	if (ncol(ml)!=ncol(mr)) stop("Numbers of columns of the two matrices must match")

	if (!is.function(col.l)) col.l <- getColorFun(ml, col.l)
	if (!is.function(col.r)) col.r <- getColorFun(mr, col.r)
	colPalLegend_l <- col.l
	if (is.character(ml) || is.factor(ml)){
		colPalLegend_l <- environment(col.l)[["mapVec"]]
	}
	colPalLegend_r <- col.r
	if (is.character(mr) || is.factor(mr)){
		colPalLegend_r <- environment(col.r)[["mapVec"]]
	}

	dummyColor <- circlize::colorRamp2(seq(0, 1, length.out=2), rep("grey", 2))
	res <- Heatmap(ml,
		col=colPalLegend_l, rect_gp=gpar(type="none"),
		cell_fun = function(j, i, x, y, w, h, fill) {
			#grid.rect(x=x, y=y, width=width, height=height, gp=gpar(col="grey", fill="grey"))
			grid.polygon(
             unit.c(x - 0.5*w, x - 0.5*w, x + 0.5*w), 
                unit.c(y - 0.5*h, y + 0.5*h, y - 0.5*h),
				gp=gpar(col=NA, fill=col.l(ml[i, j]))
			)
			grid.polygon(
                unit.c(x + 0.5*w, x + 0.5*w, x - 0.5*w), 
                unit.c(y + 0.5*h, y - 0.5*h, y + 0.5*h),
				gp=gpar(col=NA, fill=col.r(mr[i, j]))
			)
		},
		name=name.l,
    ...
	)
	# dummy heatmap for color legend
	dummyM <- matrix(rep(NA, length.out=nrow(mr)), nrow=nrow(ml), ncol=1)
	rownames(dummyM) <- rownames(ml)
	dummyHm <- Heatmap(dummyM, col=colPalLegend_r, width=unit(0, "mm"), name=name.r)

	return(res + dummyHm)
}
