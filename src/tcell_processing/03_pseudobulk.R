#!/usr/bin/env Rscript

#####################################################################
# 06_pseudobulk.R
# created on 2023-08-24 by Irem Gunduz
# Call peaks then create pseudobulk files for each cell type
#####################################################################

set.seed(12) # set seed
outputDir <- "/icbb/projects/igunduz/Tcell_subset/"

# load libraries
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(muLogR)
  library(BSgenome.Hsapiens.UCSC.hg38)
})

addArchRThreads(threads = 30) # set the cores
addArchRGenome("hg38") # set the reference genome
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)

# add group coverages
project <- ArchR::addGroupCoverages(
  ArchRProj = project,
  groupBy = "tcellClust_final"
)

# find MACS2
pathToMacs2 <- findMacs2()

# call peaks for each group
project <- ArchR::addReproduciblePeakSet(
  ArchRProj = project,
  groupBy = "tcellClust_final",
  pathToMacs2 = pathToMacs2
)

# add peak matrix and peak set
project <- addPeakMatrix(project)

saveRDS(getPeakSet(project), paste0(outputDir, "/global_peaklist.rds"))
saveArchRProject(project, outputDirectory = outputDir, load = FALSE)

# extract info for cells,cell types and sample
ArchRprojectAnnot <- data.frame(
  project$predictedGroup_tcell, project$Sample,
  project$cellNames, project$sample_exposure_group
) %>%
  dplyr::rename(
    CellType = project.predictedGroup_tcell,
    Samples = project.Sample,
    cellNames = project.cellNames,
    sample_exposure_group = project.sample_exposure_group
  ) %>%
  na.omit()

cells <- sort(unique(ArchRprojectAnnot$CellType))

fragments <- lapply(X = cells, FUN = function(cell) {
  lapply(unique(ArchRprojectAnnot$Samples), function(sample) {
    cell_names <- ArchRprojectAnnot %>%
      dplyr::filter(CellType == paste0(cell)) %>%
      dplyr::filter(Samples == paste0(sample)) %>%
      dplyr::select(cellNames) %>%
      as.vector()

    # find the indices of the cellNames
    idx <- BiocGenerics::which(ArchRprojectAnnot$cellNames %in% cell_names$cellNames)

    # extract the cells names
    cell_names <- project$cellNames[idx]
    if (NROW(cell_names) > 50) {
      logger.info(c("Extracting fragments for the cell type: ", paste0(cell), " and, sample ", paste(sample)))
      # extract the fragmens for each sample
      getFragmentsFromArrow(
        ArrowFile = paste0(
          paste0(outputDir, "/ArrowFiles/"),
          paste0(sample, ".arrow")
        ),
        cellNames = cell_names, verbose = FALSE
      )
    } else {
      logger.info("Filtered a sample that has less than 50 cells.")
    }
  }) %>% `names<-`(base::unique(ArchRprojectAnnot$Samples))
}) %>%
  `names<-`(cells)
logger.info("Extracting the samples is finished. Saving them as BED files.")

# remove NULL elements from the fragment list
fragments <- purrr::map(fragments, ~ purrr::compact(.)) %>% purrr::keep(~ length(.) != 0)

# create a directory for BED files
outputDir <- "/icbb/projects/igunduz/DARPA_analysis/BedFiles_Tcells/"
if (!dir.exists(outputDir)) {
  dir.create(outputDir)
  logger.info("Created new directory: ", outputDir)
}

# save the fragments as BED file within a specific folder seperated by cells
lapply(cells,
  FUN = function(cell) {
    if (!dir.exists(paste0(outputDir, cell))) {
      # create a folder for the cell
      dir.create(paste0(outputDir, cell))
      message("Created new directory: ", paste0(outputDir, cell, "/"))
      # save the fragments files of samples
      lapply(names(fragments[[cell]]), function(sample) {
        message("Creating a BED file for the cell type: ", paste0(cell), " and, sample ", paste(sample))
        # export fragments as bed files
        rtracklayer::export.bed(
          object = fragments[[cell]][[sample]],
          con = paste0(paste0(outputDir, cell, "/"), sample, ".bed")
        )
      })
    }
  }
)
#########