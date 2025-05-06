#################################################################
# 02_aggregate_meth.R
# created on 2023-08-29 written by Irem Gunduz
# Aggregate single-cell methylation over peak regions
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(artemis)
  library(DelayedArray)
  library(HDF5Array)
  library(readr)
  library(SummarizedExperiment)
  library(muLogR)
  library(GenomicRanges)
})
source("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/utils/allc_utils.R")
wdir <- "/icbb/projects/igunduz/DARPA_analysis/artemis_031023/"
if (!dir.exists(wdir)) {
  dir.create(wdir)
}
sdir <- "/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData"
sampleAnnot <- data.table::fread("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/sample_annot.tsv")

logger.start("Preparing sample annotation")
regionSets <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/regionsGR.rds")
regionSets <- data.table::as.data.table(regionSets)
regionSets <- regionSets[!(regionSets$seqnames %in% c("chrX", "chrY", "chrM")), ]
regionSets <- GenomicRanges::makeGRangesFromDataFrame(regionSets, keep.extra.columns = TRUE)
saveRDS(regionSets, "/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/regionsGR.rds")

regionSetNames <- "archr_peaks"
logger.completed()
logger.info(c("Number of regions:", length(regionSets)))

subset <- ImportALLC(sampleAnnot,
  filePathCol = "allC_FilePath_fixed",
  genomeAssembly = "hg38",
  sampleIdCol = "cellId",
  indexed = FALSE,
  batchSize = 300,
  numThreads = 30,
  regionSets, regionSetNames,
  # mcType = "CGN",
  blackList = NULL, singleCovOff = 99999
)

saveRDS(subset, "/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/artemis_obj.rds")

#################################################################

logger.start("Create summarized experiment object")
data <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/artemis_obj.rds")

# create a summarised experiment object
methSe <- moonToSummarizedExperiment(data, regionSetNames, isTransformed = F, regionGR = regionSets)
saveRDS(methSe, "/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/methSe.rds")

logger.completed()
#################################################################
