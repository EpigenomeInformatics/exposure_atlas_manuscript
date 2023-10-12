addPeakAnnotationsNew <- function(
  ArchRProj = NULL,
  regions = NULL,
  name = "Region",
  force = FALSE,
  logFile = createLogFile("addPeakAnnotations")
  ){

  ArchR:::.validInput(input = ArchRProj, name = "ArchRProj", valid = c("ArchRProj"))
  ArchR:::.validInput(input = regions, name = "regions", valid = c("grangeslist", "character"))
  ArchR:::.validInput(input = name, name = "name", valid = c("character"))
  ArchR:::.validInput(input = force, name = "force", valid = c("boolean"))
  ArchR:::.validInput(input = logFile, name = "logFile", valid = c("character"))

  tstart <- Sys.time()
  ArchR:::.startLogging(logFile = logFile)
  ArchR:::.logThis(mget(names(formals()),sys.frame(sys.nframe())), "addPeakAnnotations Input-Parameters", logFile = logFile)

  if(name %in% names(ArchRProj@peakAnnotation)){
    if(force){
      ArchR:::.logMessage("peakAnnotation name already exists! Overriding.", verbose = TRUE, logFile = logFile)
    }else{
      ArchR:::.logStop("peakAnnotation name already exists! set force = TRUE to override!", logFile = logFile)
    }
  }

  if(inherits(regions, "GRanges")){

    regionPositions <- GRangesList(region = regions)

  }else{

    if(is.null(names(regions))){
      names(regions) <- paste0("Region_", seq_along(regions))
    }

    if(any(duplicated(names(regions)))){
      stop("Found duplicated region names! Please make unique!")
    }

    regionPositions <- lapply(seq_along(regions), function(x){
      
      if(inherits(regions[[x]], "GRanges")){

          gr <- ArchR:::.validGRanges(regions[[x]])

      }else if(is.character(regions[[x]])){

        gr <- tryCatch({
          makeGRangesFromDataFrame(
            df = data.frame(data.table::fread(regions[[x]])), 
            keep.extra.columns = FALSE,
            seqnames.field = "V1",
            start.field = "V2",
            end.field = "V3"
          )
        }, error = function(y){

          ArchR:::.logMessage(paste0("Could not successfully get region : ", regions[[x]]), verbose = TRUE, logFile = logFile)

          if(!file.exists(regions[[x]])){
            ArchR:::.logStop(paste0("If region provided is a path it does not exist!"), logFile = logFile)
          }
          
          ArchR:::.logStop("Could not create GRanges from region", logFile = logFile)

        })

      }else{
        
        .logStop("Unrecognized input in regions please input GRanges, GRangesList, or Paths to bed files!", logFile = logFile)
      
      }

      gr

    }) %>% GRangesList

    names(regionPositions) <- names(regions)

  }

  #############################################################
  # Peak Overlap Matrix
  #############################################################
  peakSet <- getPeakSet(ArchRProj)
  if(is.null(peakSet)){
    .logStop("peakSet is NULL. You need a peakset to run addMotifAnnotations! See addReproduciblePeakSet!", logFile = logFile)
  }
  allPositions <- unlist(regionPositions, use.names=FALSE)

  ArchR:::.logDiffTime("Creating Peak Overlap Matrix", t1 = tstart, verbose = TRUE, logFile = logFile)

  overlapRegions <- findOverlaps(peakSet, allPositions, ignore.strand=TRUE)
  if(length(overlapRegions) == 0){
    stop("No Overlaps Found between regions and peak Matrix!")
  }
  ArchR:::.logThis(overlapRegions, "overlapRegions", logFile = logFile)

  regionMat <- Matrix::sparseMatrix(
    i = queryHits(overlapRegions),
    j = match(names(allPositions),names(regionPositions))[subjectHits(overlapRegions)],
    x = rep(TRUE, length(overlapRegions)),
    dims = c(length(peakSet), length(regionPositions))
  )
  colnames(regionMat) <- names(regionPositions)
  ArchR:::.logThis(regionMat, "regionMat", logFile = logFile)

  regionMat <- SummarizedExperiment::SummarizedExperiment(assays=SimpleList(matches = regionMat), rowRanges = peakSet)
  ArchR:::.logThis(regionMat, "regionSE", logFile = logFile)

  #############################################################
  # Filter Regions With No Matches
  #############################################################

  #Number of Overlaps
  nO <- Matrix::colSums(assay(regionMat))
  rF <- names(which(nO == 0))

  if(all(nO == 0)){
    stop("No Overlaps Found! Please check your peakSet and genome!")
  }

  if(length(rF) > 0){
   ArchR:::.logDiffTime(paste0("Filtering Region Annotations with 0 overlaps :\n\n ", paste(rF, collapse=", "), "\n\n"), t1 = tstart, verbose = TRUE, logFile = logFile)
    #Filter
    regionPositions <- regionPositions[!(names(regionPositions) %in% rF)]
    regionMat <- regionMat[,names(regionPositions),drop=FALSE]
  }else{
    ArchR:::.logDiffTime(paste0("All Regions Overlap at least 1 peak!"), t1 = tstart, verbose = TRUE, logFile = logFile)
  }  

  #############################################################
  # Summarize and Save
  #############################################################

  dir.create(file.path(getOutputDirectory(ArchRProj), "Annotations"), showWarnings=FALSE)
  savePositions <- file.path(getOutputDirectory(ArchRProj), "Annotations", paste0(name,"-Positions-In-Peaks.rds"))
  saveMatches <- file.path(getOutputDirectory(ArchRProj), "Annotations", paste0(name,"-Matches-In-Peaks.rds"))

  out <- SimpleList(
      regionMatches = regionMat,
      regionPositions = regionPositions,
      date = Sys.Date()
    )

  ArchRProj@peakAnnotation[[name]]$Name <- name
  ArchRProj@peakAnnotation[[name]]$Positions <- savePositions
  ArchRProj@peakAnnotation[[name]]$Matches <- saveMatches

  ArchR:::.safeSaveRDS(out, file.path(getOutputDirectory(ArchRProj),  "Annotations", paste0(name,"-In-Peaks-Summary.rds")), compress = FALSE)
  ArchR:::.safeSaveRDS(out$regionPositions, savePositions, compress = FALSE)
  ArchR:::.safeSaveRDS(out$regionMatches, saveMatches, compress = FALSE)

  return(ArchRProj)

}