processRawALLC <- function(
    file, blackList, excludeChr, context = "^CG", assembly = "hg38", mergeStrands = TRUE,
    region, singleCovOff, returnMethLev = TRUE) {
  # read the allc file
  meth <- data.table::fread(file) # muRtools::readTab(file, header = FALSE)

  if (ncol(meth) < 7) {
    meth <- data.table::data.table(matrix(nrow = length(region), ncol = 2))
    colnames(meth) <- c("methLev", "cov")
    meth_mean <- N_sites <- N_valid_sites <- N_valid_sites <- 0
    # regId <- rep(0, length(region))
    # TODO check if it is correct
  } else {
    colnames(meth) <- c("seqnames", "start", "strand", "context", "mc", "cov", "isMeth")
    N_sites <- NROW(meth) # Save the number of sites in the file
    if (!is.null(excludeChr)) {
      meth <- meth[!(meth[["seqnames"]] %in% excludeChr), ]
    }
    # meth <- meth[(meth[["context"]] %in% mcType), ]
    cntLen <- nchar(meth[, "context"])
    if (context == "^CG") {
      cntLen <- rep(2L, nrow(meth))
    }
    meth[, "end"] <- meth[, "start"]
    isNeg <- meth[, "strand"] == "-"
    meth <- as.data.frame(meth)
    meth[!isNeg, "end"] <- meth[!isNeg, "end"] + cntLen[!isNeg] - 1L
    meth[isNeg, "start"] <- meth[isNeg, "start"] - cntLen[isNeg] + 1L

    # meth <- muRtools::df2granges(meth,
    # ids = NULL,
    # chrom.col = "seqnames",
    # start.col = "start",
    # end.col = "end",
    # strand.col = "strand",
    # coord.format = "B1RI",
    # assembly = assembly,
    # doSort = TRUE,
    # adjNumChromNames = TRUE
    # )
    meth <- GenomicRanges::GRanges(
      seqnames = meth$seqnames,
      ranges = IRanges::IRanges(start = meth$start, end = meth$end),
      strand = meth$strand,
      mc = meth$mc,
      cov = meth$cov
    )
    if (mergeStrands) {
      tmp <- meth
      GenomicRanges::elementMetadata(tmp) <- NULL
      GenomicRanges::strand(tmp) <- "*"
      tmp <- unique(tmp)
      GenomicRanges::elementMetadata(tmp) <- data.frame(
        mc  = rep(as.integer(NA), length(tmp)),
        cov = rep(as.integer(NA), length(tmp))
      )
      oo <- GenomicRanges::findOverlaps(meth, tmp, ignore.strand = TRUE)
      mcS <- tapply(GenomicRanges::elementMetadata(meth)[, "mc"], S4Vectors::subjectHits(oo), sum)
      GenomicRanges::elementMetadata(tmp)[as.integer(names(mcS)), "mc"] <- mcS
      covS <- tapply(GenomicRanges::elementMetadata(meth)[, "cov"], S4Vectors::subjectHits(oo), sum)
      GenomicRanges::elementMetadata(tmp)[as.integer(names(covS)), "cov"] <- covS
      meth <- tmp
      rm(tmp)
    }

    N_valid_sites <- length(meth) # Save the number of sites after filtering
    olaps <- GenomicRanges::findOverlaps(region, meth, ignore.strand = TRUE)
    subject_hits <- S4Vectors::subjectHits(olaps)
    if (NROW(olaps) == 0) {
      meth <- data.table::data.table(matrix(nrow = length(region), ncol = 2))
      colnames(meth) <- c("methLev", "cov")
      meth_mean <- N_sites <- N_valid_sites <- N_valid_sites <- 0
    } else {
      meth <- data.table::as.data.table(data.frame(
        id = S4Vectors::queryHits(olaps),
        GenomicRanges::elementMetadata(meth)[subject_hits, ],
        stringsAsFactors = FALSE
      ))

      counts_sum <- meth[, lapply(.SD, sum, na.rm = TRUE), by = id]
      region <- as.data.frame(region) # data.table::as.data.table(region))
      meth <- data.table::data.table(
        regions = paste0(
          region[counts_sum$id, "seqnames"], "_",
          region[counts_sum$id, "start"], "to",
          region[counts_sum$id, "end"]
        ),
        counts_sum[, c(2:(ncol(counts_sum))), with = FALSE],
        stringsAsFactors = FALSE
      )[cov > 0 & cov < singleCovOff]

      if (returnMethLev) {
        region <- data.table::as.data.table(region)
        region <- region[, regions := paste0(seqnames, "_", start, "to", end)][, .(regions)]
        # regId <- region$regions

        meth[, methLev := mc / cov]
        # fill the missing regions
        meth <- merge(region, meth, by = "regions", all = TRUE, sort = FALSE)

        # meth[is.na(meth)] <- 0  # Set missing values to 0
        meth <- meth[, .(methLev, cov)]
        meth_mean <- mean(meth$methLev, na.rm = TRUE)
      } else {
        meth_mean <- mean(meth$mc, na.rm = TRUE) / mean(meth$cov, na.rm = TRUE)
      }
    }
  }

  sampleStats <- data.frame(
    fileId = file,
    N_sites_cov1 = sum(meth[, "cov"] == 1, na.rm = TRUE),
    N_sites_cov2plus = sum(meth[, "cov"] > 1, na.rm = TRUE),
    N_sites_cov5plus = sum(meth[, "cov"] > 4, na.rm = TRUE),
    covg_cum = mean(meth$cov, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
  sampleStats <- cbind(sampleStats, meth_mean = meth_mean, N_sites = N_sites, N_valid_sites = N_valid_sites)
  return(list(data = meth, stats = sampleStats)) # Return processed data and QC statistics
}

mcParseR <- function(pattern) {
  IUPAC_table <- c(
    "A" = "A",
    "T" = "T",
    "C" = "C",
    "G" = "G",
    "R" = "AG",
    "Y" = "CT",
    "S" = "GC",
    "W" = "AT",
    "K" = "GT",
    "M" = "AC",
    "B" = "CGT",
    "D" = "AGT",
    "H" = "ATC",
    "V" = "ACG",
    "N" = "ATCGN"
  )
  if (pattern == "CGN") {
    return(c("CGA", "CGC", "CGG", "CGT"))
  } else {
    positions <- vector()
    pattern <- toupper(pattern)
    pattern <- unlist(strsplit(pattern, split = ""))
    for (base in pattern) {
      pos <- as.character(IUPAC_table[base])
      if (nchar(pos) > 3) {
        positions <- c(positions, (unlist(strsplit(pos, split = ""))))
      } else {
        positions <- c(positions, pos)
      }
    }
    if (NROW(positions) > 3) {
      beg <- positions[1:2]
      pos_list <- vector()
      for (i in positions[-c(1:2)]) {
        pos_list <- c(pos_list, stringr::str_flatten(beg, i))
      }
      return(pos_list)
    } else {
      return(stringr::str_flatten(positions))
    }
  }
}
moonToSingleCellExperiment <- function(.Object, region, regionGR = NULL, isTransformed = FALSE) {
  if (is.null(region)) {
    stop("Please specify a region.")
  }
  available_mat <- getAvailableCountMatrices(.Object)
  if (!(paste0(region, ".methLev") %in% available_mat)) {
    stop("The specified region is not available in the moon Object.")
  }
  if (!region %in% getRegionTypes(.Object)) {
    stop("The specified region is not available in the moon Object.")
  }
  if (isTransformed) {
    assay <- list(
      mc = .Object@countTransform[[paste0(region, ".methLev")]],
      cov = .Object@countTransform[[paste0(region, ".cov")]]
    )
  } else {
    assay <- list(
      mc = .Object@counts[[paste0(region, ".methLev")]],
      cov = .Object@counts[[paste0(region, ".cov")]]
    )
  }
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = assay,
    colData = data.frame(samples = getSamples(.Object, region, isTransformed = isTransformed)),
    metadata = .Object@sampleAnnot,
    mainExpName = region
  )
  if (!is.null(regionGR) && class(regionGR) == "GRanges") {
    SummarizedExperiment::rowData(sce) <- regionGR
  }
  return(sce)
}

ImportALLC <- function(sampleAnnot, filePathCol, genomeAssembly, sampleIdCol = NULL,
                       indexed = FALSE, batchSize = 25, numThreads = 1,
                       regionSets = NULL, regionSetNames = NULL, context = "^CG",
                       promoterRegion = c(2000, 100),
                       blackList = NULL, singleCovOff = 99999,
                       excludeChr = c("chrX", "chrY", "chrM")) {
  assertthat::assert_that(!is.null(sampleAnnot) && is.data.frame(sampleAnnot),
    msg = "Please provide sampleAnnot as a data.frame."
  )
  assertthat::assert_that(filePathCol %in% colnames(sampleAnnot),
    msg = "Please provide valid filePathCol"
  )
  assertthat::assert_that(sampleIdCol %in% colnames(sampleAnnot),
    msg = "Please provide valid sampleIdCol"
  )
  assertthat::assert_that(!is.null(genomeAssembly) && is.character(genomeAssembly),
    msg = "Please provide a genome assembly using genomeAssembly parameter."
  )
  if (!is.null(excludeChr) && !is.character(excludeChr)) {
    message("Invalid excludeChr value detected. Ignoring parameter excludeChr.")
    excludeChr <- NULL
  }
  if (is.null(singleCovOff) || !is.numeric(singleCovOff)) {
    message("WARNING! Invalid singleCovOff input. Using the default parameters.")
    singleCovOff <- 99999
  }
  assertthat::assert_that(is.numeric(batchSize) || length(batchSize) > 1,
    msg = "Invalid batchSize input."
  )
  if (is.null(context) || !is.character(context)) {
    message("WARNING! Invalid context input. Using the default parameters.")
    context <- "^CG"
  }
  if (!is.numeric(numThreads) | length(numThreads) > 1) {
    message("WARNING! Invalid numThreads input. Using the default parameters.")
    numThreads <- 1
  }
  if (length(promoterRegion) != 2) {
    message("WARNING! Invalid promoterRegion input. Using the default parameters.")
    promoterRegion <- c(2000, 100)
  }
  if (batchSize > 100 && numThreads == 1) {
    message("WARNING: batchSize is larger than 100! And no parallelization registered.
            This may cause memory issues.")
  }
  # get sample path info from sampleannot
  inputFilenames <- as.character(sampleAnnot[[filePathCol]])
  sampleIds <- as.character(sampleAnnot[[sampleIdCol]])
  if (indexed) {
    extend <- ".tbi"
  } else {
    extend <- ".tsv"
  }
  # make sure we have the final paths
  inputFilenames <- sampleAnnot[grepl(extend, inputFilenames), ] %>%
    dplyr::select(dplyr::all_of(filePathCol))
  inputFilenames <- as.character(sampleAnnot[[filePathCol]])
  names(inputFilenames) <- sampleIds

  # check if all files exists
  logger::log_info("Checking if all input files exists.")
  if (!all(file.exists(inputFilenames))) {
    missingSamples <- sampleIds[!file.exists(inputFilenames)]
    stop("Missing input files for samples: ", paste(missingSamples, collapse = ", "))
  }
  logger::log_success("Located all of the samples.")
  if (is.null(regionSets)) {
    logger::log_info("Preparing default region sets: tiling5kbp and, promoter regions.")
    regionSets <- getDefaultRegions(genomeAssembly, promoterRegion, excludeChr)
    regionSetNames <- names(regionSets)
    logger::log_success()
  }

  if (!is.null(regionSets) && !is.null(regionSetNames) &&
    (is.character(regionSetNames))) {
    assertthat::assert_that(
      class(regionSets) %in% c(
        "data.frame", "GRanges",
        "list", "GRangesList"
      ),
      msg = "Please provide valid regionSets imput object."
    )

    if (class(regionSets) %in% c("list", "GRangesList")) {
      assertthat::assert_that(length(regionSetNames) == length(regionSets),
        msg = "Provided regionSetNames length doesn't match with the number of
             objects within the list."
      )
    }
    if (is.data.frame(regionSets) || class(regionSets) == "GRanges") {
      regionSets <- data.table::as.data.table(regionSets)
      regionSets <- regionSets[, regions := paste0(seqnames, ":", start, "-", end)]
      regionSets <- regionSets %>% distinct(regions, .keep_all = TRUE)
      regionSets <- dplyr::select(regionSets, !regions)
      regionSets <- list(GenomicRanges::makeGRangesFromDataFrame(regionSets, keep.extra.columns = TRUE))
    }
    names(regionSets) <- regionSetNames
  }

  # TODO add check for GRangesList of regions

  if (is.null(blackList)) { # if present, get the blacklist regions
    blackList <- artemis:::blackListeR(assembly = genomeAssembly, path = NULL)
  } else {
    blackList <- artemis:::blackListeR(
      assembly = genomeAssembly,
      path = blackList
    )
  }

  logger::log_info("Creating moon object")
  obj <- artemis:::moon(sampleAnnot, regionSetNames, genomeAssembly, regionSets)
  logger::log_success("Created a moon object.")

  for (rt in names(regionSets)) {
    logger::log_info("Including region set: ", rt)
    if (!dir.exists(paste0("data/", rt))) {
      dir.create(paste0("data/", rt))
    }
    logger::log_info("Initializing matrix for ", rt)

    # DelayedArray::setAutoRealizationBackend("HDF5Array") # do we need this?
    sinks <- lapply(c("methLev", "cov"), function(t) {
      # TODO convert this to data.table operations
      regIds <- as.data.frame(regionSets[[rt]]) %>%
        dplyr::mutate(regIds = paste0(seqnames, ":", start, "-", end)) %>%
        dplyr::select(regIds)

      sink <- HDF5Array::HDF5RealizationSink(
        dim = c(NROW(regionSets[[rt]]), length(inputFilenames)),
        dimnames = list(as.vector(regIds$regIds), sampleIds),
        type = ifelse(t == "methLev", "double", "integer"),
        filepath = paste0("data/", rt, "/", t, ".h5"),
        name = rt, level = 6
      )
      return(sink)
    })
    names(sinks) <- c("methLev", "cov")

    # Chunk the files into batches
    if (numThreads > 1) {
      chunkfiles <- artemis:::chunkeR(inputFilenames,
        batchSize = batchSize,
        numCores = numThreads
      )
    } else {
      chunkfiles <- artemis:::chunkeR(files = inputFilenames, batchSize = batchSize)
    }
    logger::log_info("Created ", length(chunkfiles), " chunks from inputs.")

    allStats <- data.frame() # initialize the sample & region statistics

    logger::log_info("Starting parallel backend.")
    # register parallel backend
    cl <- parallel::makeCluster(numThreads)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
    parallel::clusterExport(cl, varlist = c(
      "processRawALLC", "blackList", "excludeChr",
      "regionSets", "rt", "singleCovOff", "context", "genomeAssembly"
    ), envir = environment())

    # set the grid
    grid <- DelayedArray::ArbitraryArrayGrid(list(
      NROW(regionSets[[rt]]),
      cumsum(lengths(chunkfiles))
    ))
    for (chunk in seq_along(chunkfiles)) {
      logger::log_info("Prosessing chunk ", chunk, " of ", length(chunkfiles))
      reads <- parallel::parLapply(cl, chunkfiles[[chunk]],
        fun = processRawALLC, blackList = blackList,
        excludeChr = excludeChr, # mcType = mcType,
        region = regionSets[[rt]], context = context,
        assembly = genomeAssembly,
        singleCovOff = singleCovOff
      )

      # combine the results into a single data.table
      mreads <- lapply(c("methLev", "cov"), function(x) {
        lapply(names(reads), function(f) {
          reads[[f]]$data %>%
            dplyr::select(dplyr::all_of(x))
        }) %>%
          rlist::list.cbind() %>%
          data.table::as.data.table()
      }) # TODO find a better way to do this
      names(mreads) <- c("methLev", "cov")

      # combine the statistics
      stats <- data.table::as.data.table(rlist::list.rbind(lapply(names(reads), function(f) {
        reads[[f]]$stats
      })))
      allStats <- rbind(allStats, stats)

      # combine the region statistics
      # regs <- data.table::as.data.table(rlist::list.cbind(lapply(names(reads), function(f) {
      #  reads[[f]]$regStats
      # })))
      # regionStats <- cbind(regionStats, regs)

      # write the results to the sink
      for (t in c("methLev", "cov")) {
        DelayedArray::write_block(
          block = as.matrix(mreads[[t]]),
          viewport = grid[[as.integer(chunk)]], sink = sinks[[t]]
        )
      }
      logger::log_success("Prosessed chunk ", chunk, " of ", length(chunkfiles))
      # clean memory
      artemis:::cleanMem()
      rm(reads, stats, mreads)
    }
    logger::log_info("Saving sample statistics")
    write.table(allStats, file = paste0("data/", rt, "/sampleStats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    # colnames(regionStats) <- basename(colnames(regionStats))
    # write.table(regionStats, file = paste0("data/", rt, "/regionStats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  }
  # close the sinks
  for (t in c("methLev", "cov")) {
    DelayedArray::close(sinks[[t]]) # close the sink
    obj@counts[[rt]][[t]] <- as(sinks[[t]], "DelayedArray") # convert to DelayedArray
  }
  obj@counts <- unlist(obj@counts)
  logger::log_success("Created matrices for the region set: ", rt)
  artemis:::cleanMem()
  if (length(names(regionSets)) == 1) {
    data.table::setnames(allStats, old = "fileId", new = filePathCol)
    obj@sampleAnnot <- merge(sampleAnnot, allStats, by = filePathCol)
  }
  return(obj)
}
