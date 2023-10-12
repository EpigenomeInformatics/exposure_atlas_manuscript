suppressPackageStartupMessages({
	library(dplyr)
	library(artemis)
	library(muLogR)
	#library(data.table)
	library(SummarizedExperiment)
    library(DelayedMatrixStats)
	library(HDF5Array)
})
set.seed(43)
meth_cellfilter_nsites_range <- c(5e5, 4e6)
filter_region_ncells <- 0.215

logger.start("Quality check")
logger.start("Loading data")
#create a summarised experiment object
methSe <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/methSe.rds")
logger.completed()


logger.start("Checking QC stats and removing low quality samples")
#check QC stats
qc_stats <- as.data.frame(methSe@metadata)
filteredCellIdx <- rep(TRUE, nrow(qc_stats))
nSites <- qc_stats[, "N_valid_sites"]
survive <- nSites >= meth_cellfilter_nsites_range[1] & nSites <= meth_cellfilter_nsites_range[2]
logger.info(paste0("", sum(survive), " of ", sum(filteredCellIdx), " (", round(100*sum(survive)/sum(filteredCellIdx), 2), "%) cells retained after filtering for nSites"))
filteredCellIdx <- filteredCellIdx & survive
qc_stats <- qc_stats[filteredCellIdx,]
methSe <- methSe[,filteredCellIdx] 
nCells <- nrow(qc_stats)
logger.completed()


logger.info("Retrieving the matrix")
methSeL <- assay(methSe, "mc")

logger.start("Binaring the matrix")
methSeL <- !is.na(methSeL) & methSeL < 0.4
logger.completed()

logger.start("Filtering low quality regions")
keep <- DelayedMatrixStats::rowSums2(methSeL) 
regIdmc <- keep > 0
logger.info(c("Removing ", sum(!regIdmc), " of ", nrow(methSeL), " (", round(100 * sum(!regIdmc) / nrow(methSeL), 2), "%) regions because they are unmethylated in all cells"))
logger.info(paste0("Retained ", sum(regIdmc), " of ", nrow(methSeL), " (", round(100 * sum(regIdmc) / nrow(methSeL), 2), "%) regions"))
methSe <- methSe[regIdmc, ]
logger.completed()


logger.start("Saving filtered regions for later usage...")
rr <- data.table::as.data.table(readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/regionsGR.rds"))
rr <- as.data.frame(rr)
rr <- rr[regIdmc, ]


logger.start("Filtering based on coverage")
regCov <- assay(methSe, "cov")
minNcells <- ceiling(filter_region_ncells * nCells)
logger.start("Computing cell coverages")
regCov <- regCov > 0 #!is.na(regCov) & regCov > 0
keep <- DelayedMatrixStats::rowSums2(regCov,na.rm=TRUE)
logger.completed()
regIdx <- keep >= minNcells
nRem <- sum(!regIdx)
nKeep <- sum(regIdx)
logger.info(c("Removing", nRem, " of ", length(regIdx),paste0("(", round(100*nRem/length(regIdx), 2), "%)"),
	"regions", "because they are covered in fewer than",minNcells, "cells",paste0("(",filter_region_ncells*100, "%)")))
logger.info(c("Retained", nKeep, "regions"))
methSe <- methSe[regIdx, ]
logger.completed()
rr <- rr[regIdx, ]
saveRDS(rr, "/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/regionsGR_filtered.rds")
logger.completed()

logger.info("Saving filtered annotation dataset...")
saveRDS(qc_stats, "/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/cellAnnot_meth.rds")
logger.info("Saving filtered summarised experiment object")
saveRDS(methSe, "/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/methSe_filtered.rds")
logger.completed()



