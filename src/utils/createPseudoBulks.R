# A function to create pseudobulks per one or two group
# Irem B. Gunduz
#
createPseudoBulks <- function(sampleAnnot, filePathCol, sampleIdCol = NULL, numThreads = 1, mcType = "CGN",
                              groupName, groupName2, fileType = "allc", singleCovOff = 99999, singleBP = FALSE,
                              indexed = FALSE, excludeChr = c("chrX", "chrY", "chrM", "chrL"), outputDir = NULL,
                              regionGR = NULL) {
  assertthat::assert_that(isFALSE(singleBP) || isTRUE(singleBP),
    msg = "Please provide TRUE/FALSE value for singleBP"
  )
  if (!singleBP && class(regionGR) != "GRanges") {
    stop("Please provide regionSet as GRanges object")
  }
  assertthat::assert_that(!is.null(sampleAnnot) && is.data.frame(sampleAnnot),
    msg = "Please provide sampleAnnot as a data.frame."
  )
  assertthat::assert_that(filePathCol %in% colnames(sampleAnnot),
    msg = "Please provide valid filePathCol input"
  )
  assertthat::assert_that(sampleIdCol %in% colnames(sampleAnnot),
    msg = "Please provide valid sampleIdCol input"
  )
  assertthat::assert_that(groupName %in% colnames(sampleAnnot),
    msg = "Please provide valid groupName input"
  )
  assertthat::assert_that(!is.null(fileType),
    msg = "Please provide valid fileType input"
  )
  if (!is.null(excludeChr) && !is.character(excludeChr)) {
    message("Invalid excludeChr value detected. Ignoring parameter excludeChr.")
    excludeChr <- NULL
  }
  if (is.null(singleCovOff) || !is.numeric(singleCovOff)) {
    message("WARNING! Invalid singleCovOff input. Using the default parameters.")
    singleCovOff <- 99999
  }
  if (!is.numeric(numThreads) | length(numThreads) > 1) {
    message("WARNING! Invalid numThreads input. Using the default parameters.")
    numThreads <- 1
  }
  if (!is.null(outputDir) && !is.character(outputDir)) {
    message("Invalid outputDir submitted. Switching to default settings.")
    outputDir <- NULL
  }
  if (!is.null(outputDir) && is.character(outputDir) && !dir.exists(outputDir)) {
    message("Couldn't locate outputDir. Switching to default settings.")
    outputDir <- NULL
  }
  if ((is.null(mcType) || length(mcType) != 1) && fileType == "allc") {
    stop("Only one mcType must be provided for allc files!")
  }
  if (is.null(outputDir) || !is.character(outputDir)) {
    if (!dir.exists("data")) {
      dir.create("data")
    }
    if (!dir.exists("data/pseudoBulks")) {
      dir.create("data/pseudoBulks")
    }
    message("Created a directory at data/pseudoBulks")
  }
  # start parallization
  cl <- parallel::makeCluster(numThreads)
  doParallel::registerDoParallel(cl)
  on.exit(parallel::stopCluster(cl))

  for (group in unique(sampleAnnot[[groupName]])) {
    sannot <- sampleAnnot[sampleAnnot[[groupName]] %in% group]
    if (!is.null(groupName2)) {
      for (group2 in unique(sampleAnnot[[groupName2]])) {
        if (!file.exists(paste0(outputDir, "/", group, "_", group2, ".bedGraph"))) {
          logger::log_info("Processing ", group2, " for ", group)
          files <- sannot[sannot[[groupName2]] %in% group2][[filePathCol]]
          if (length(files) > 100) {
            logger::log_info("Found ", length(files), " files. Reading the information may take a while...")
          } else {
            logger::log_info("Found ", length(files), " files.")
          }
          if (length(files) > 0) {
            if (fileType == "allc" && !indexed) {
              # aggregate allcs
              aggregateALLCSamples(
                files, excludeChr, singleCovOff, group, group2,
                mcParseR(mcType), outputDir, singleBP, regionGR, cl
              )
            }
          }
        }
      }
    } else {
      sannot <- sampleAnnot[sampleAnnot[[groupName]] %in% group]
      files <- sannot[[filePathCol]]
      if (length(files) > 0 && !file.exists(paste0(outputDir, group, ".bedGraph"))) {
        if (fileType == "allc" && !indexed) {
          # aggregate allcs
          aggregateALLCSamples(files, excludeChr, singleCovOff, group,
            group2 = NULL,
            mcParseR(mcType), outputDir, singleBP, regionGR, cl
          )
        }
      }
    }
  }
  logger::log_success("Successfully created all the samples.")
}
# helper for pseudobulks, reads files, converts to data.table to GRanges and exports as a bed file
# Irem B. Gunduz
#
aggregateALLCSamples <- function(files, excludeChr, singleCovOff, group, group2,
                                 mcType, outputDir, singleBP, regionGR = NULL, cl) {
  if (singleBP) {
    parallel::clusterExport(cl, varlist = c(
      "excludeChr", "mcType", "singleCovOff", "%>%", "%in%"
    ), envir = environment())
  } else {
    parallel::clusterExport(cl, varlist = c(
      "excludeChr", "mcType", "singleCovOff", "%>%", "%in%", "regionGR"
    ), envir = environment())
  }
  if (length(files) > 50) {
    dir.create(paste0(outputDir, "/pseudobulkTemp")) # create a temp dir to save files
    cfiles <- artemis:::chunkeR(files, size = 50) # create chuncks from the files
    names(cfiles) <- sapply(1:length(cfiles), function(i) paste0("chunk_", i))
    logger::log_info("Created ", length(cfiles), " chunks from input.")
    for (files in names(cfiles)) {
      logger::log_info("Processing ", files)
      if (singleBP) {
        meth <- parallel::parLapply(cl, cfiles[[files]], fun = function(file) {
          tab <- data.table::fread(file) %>%
            dplyr::filter(V4 %in% mcType, !V1 %in% excludeChr, V6 < singleCovOff)
          tab <- tab[grep("chr", tab$V1), ] # this will return only main chr
          tab <- tab %>%
            dplyr::rename(seqnames = V1, start = V2, mc = V5, cov = V6) %>%
            dplyr::select(seqnames, start, mc, cov)
          return(tab)
        })
      } else {
        meth <- parallel::parLapply(cl, files, fun = function(file) {
          tab <- artemis:::processRawALLC(file,
            blackList = NULL, excludeChr = excludeChr,
            region = regionGR, mcType = mcType, resizeWidth = NULL,
            singleCovOff = singleCovOff, returnMethLev = FALSE
          )$data
          return(tab)
        })
      }
      logger::log_info("Summarising info for ", files)
      meth <- rlist::list.rbind(meth)
      meth <- meth[, lapply(.SD, sum, na.rm = T), by = .(seqnames, start), .SDcols = c("mc", "cov")] %>%
        dplyr::filter(cov < singleCovOff)
      arrow::write_parquet(x = meth, sink = paste0(outputDir, "/pseudobulkTemp/", files, ".gz.parquet"), compression = "gzip", compression_level = 5)
      rm(meth)
    }
    rm(cfiles, files)
    logger::log_info("Retrieved the information. Now summarizing...")
    meth <- arrow::open_dataset(paste0(outputDir, "/pseudobulkTemp/"), unify_schemas = TRUE) %>%
      dplyr::group_by(seqnames, start) %>%
      dplyr::summarise(mc = sum(mc), cov = sum(cov)) %>%
      dplyr::filter(cov < singleCovOff) %>%
      dplyr::mutate(cov = cov - mc) %>%
      dplyr::collect()
  } else {
    if (singleBP) {
      meth <- parallel::parLapply(cl, files, fun = function(file) {
        tab <- data.table::fread(file) %>%
          dplyr::filter(V4 %in% mcType, !V1 %in% excludeChr, V6 < singleCovOff)
        tab <- tab[grep("chr", tab$V1), ] # this will return only main chr
        tab <- tab %>%
          dplyr::rename(seqnames = V1, start = V2, mc = V5, cov = V6) %>%
          dplyr::select(seqnames, start, mc, cov)
        return(tab)
      })
    } else {
      meth <- parallel::parLapply(cl, files, fun = function(file) {
        tab <- artemis:::processRawALLC(file,
          blackList = NULL, excludeChr = excludeChr,
          region = regionGR, mcType = mcType, resizeWidth = NULL,
          singleCovOff = singleCovOff, returnMethLev = FALSE
        )
        return(tab)
      })
    }
    meth <- rlist::list.rbind(meth)
    meth <- meth[, lapply(.SD, sum, na.rm = T), by = .(seqnames, start), .SDcols = c("mc", "cov")] %>%
      dplyr::filter(cov < singleCovOff)
  }
  makeGRanges(meth, group, group2, outputDir)
}

# A helper for pseudobulks,converts and exports GRanges objects
#
makeGRanges <- function(meth, group, group2, outputDir) {
  # Add empty strand column
  meth$strand <- rep("*", NROW(meth$start))
  meth <- meth %>% dplyr::select(seqnames, start, strand, mc, cov)

  if (is.null(group2)) {
    logger::log_success("Created a GRanges object for ", group, ". Now writing into a bedGraph file.")
    data.table::fwrite(meth, paste0(outputDir, "/", group, ".bedGraph"), TRUE, FALSE,
      sep = "\t", row.names = FALSE, col.names = FALSE, showProgress = TRUE
    )
    logger::log_success("Wrote ", group, " pseudobulks into a bedGraph file.")
  } else {
    logger::log_success("Created a GRanges object for ", group2, " of ", group, ". Now writing into a bedGraph file.")
    data.table::fwrite(meth, paste0(outputDir, "/", group, "_", group2, ".bedGraph"), TRUE, FALSE,
      sep = "\t", row.names = FALSE, col.names = FALSE, showProgress = TRUE
    )
    logger::log_success("Wrote ", group2, " for ", group, " pseudobulks into a bedGraph file.")
  }
  logger::log_info("Cleaning up memory and/or temp files...")
  artemis:::cleanMem()
  if (dir.exists(paste0(outputDir, "/pseudobulkTemp/"))) {
    unlink(list.files(paste0(outputDir, "/pseudobulkTemp/"), full.names = TRUE), recursive = TRUE) # Remove the temp files.
  }
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
