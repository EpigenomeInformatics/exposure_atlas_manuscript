#!/usr/bin/env Rscript

#####################################################################
# 12_cell_type_chraccr.R
# created on 2023-10-25 by Irem Gunduz
# Run vanilla ChrAccR analysis for cell-type level comparisons
#####################################################################

suppressPackageStartupMessages({
  library(ArchR)
  library(ChrAccR)
  library(dplyr)
  library(muLogR)
  library(muRtools)
})
set.seed(12) # set seed
outputDir <- "/icbb/projects/igunduz/archr_project_011023"
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)

# read the sample annotation
sampleannot <- read.delim("/icbb/projects/igunduz/sampleannot.tsv") %>%
  dplyr::filter(sample_exposure_type %in% c("HIV", "OP"))
sampleannot$fragmentFiles <- gsub(x = sampleannot$fragmentFiles, pattern = ".bed", replacement = ".tsv.gz")

# set directory for the output
bedDir <- "/icbb/projects/igunduz/DARPA_analysis/BedFiles_final/Cell_Types"
rundir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_2023-10-27_cell_type"
# paste0("/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_", Sys.Date(), "_cell_type")
# if (!dir.exists(rundir)) dir.create(rundir)


cellAnnot_atac <- read.delim("/icbb/projects/igunduz/artemis_run/final_moon_240823/cellColData.tsv")
rownames(cellAnnot_atac) <- cellAnnot_atac$cellId_archr

colname_samplId_atac <- "sample_sampleId_cminid"
colname_samplId_meth <- "CommonMinID"

cellAnnot_meth <- read.delim("/icbb/projects/igunduz/DARPA/allc_full_sampleannot_290823.tsv")
rownames(cellAnnot_meth) <- cellAnnot_meth[, "Cell_UID"]
cellAnnot_meth <- cellAnnot_meth %>%
  dplyr::filter(sample_exposure_type %in% c("COVID", "FLU", "HIV", "OP"))


logger.start("Finding common samples between METH-ATAC")
sampleIds_atac <- sort(unique(cellAnnot_atac[, colname_samplId_atac]))
sampleIds_meth <- sort(unique(cellAnnot_meth[, colname_samplId_meth]))
summarizeSetOverlap(sampleIds_atac, sampleIds_meth, set1name = "ATAC", set2name = "BS", doVenn = FALSE)
sampleIds_int <- intersect(sampleIds_atac, sampleIds_meth)

nCells_perSample_atac <- sort(table(cellAnnot_atac[cellAnnot_atac[, colname_samplId_atac] %in% sampleIds_int, colname_samplId_atac]))
nCells_perSample_meth <- sort(table(cellAnnot_meth[cellAnnot_meth[, colname_samplId_meth] %in% sampleIds_int, colname_samplId_meth]))
# exclude samples with too few cells
idx <- nCells_perSample_atac[sampleIds_int] >= 200 & nCells_perSample_meth[sampleIds_int] >= 50
logger.info(paste0("excluding ", sum(!idx), " of ", length(idx), " (", round(100 * sum(!idx) / length(idx), 2), "%) samples because they do not contain enough cells"))
sampleIds_int <- sampleIds_int[idx]

idx <- cellAnnot_atac[, colname_samplId_atac] %in% sampleIds_int
cellIds_atac <- rownames(cellAnnot_atac)[idx]
logger.info(paste0("Retained ", sum(idx), " of ", nrow(cellAnnot_atac), " (", round(100 * sum(idx) / nrow(cellAnnot_atac), 2), "%) ATAC cells for shared samples"))
idx <- cellAnnot_meth[, colname_samplId_meth] %in% sampleIds_int
cellIds_meth <- rownames(cellAnnot_meth)[idx]
logger.info(paste0("Retained ", sum(idx), " of ", nrow(cellAnnot_meth), " (", round(100 * sum(idx) / nrow(cellAnnot_meth), 2), "%) METH cells for shared samples"))

cellAnnot_atac <- cellAnnot_atac[cellIds_atac, ]
cellAnnot_meth <- cellAnnot_meth[cellIds_meth, ]

logger.completed()

shared_atac_samples <- unique(cellAnnot_atac$sample_sampleId)
sampleannot <- dplyr::filter(sampleannot, sampleIdCol %in% shared_atac_samples)
cells <- c("B_mem", "B_naive", "DC", "Mono_CD14", "Mono_CD16", "NK_CD16", "T_mem_CD8", "T_mem_CD4", "T_mix", "T_naive", "T_mait")

# get the fragment files for the cell types
all_pseudobulks <- list.files(bedDir, full.names = T)
# all_pseudobulks <- rlist::list.rbind(unlist(lapply(cells,function(cell){list.files(paste0(bedDir,cell),full.names=T)})))
# all_pseudobulks <- as.vector(all_pseudobulks)

# reassign the fragment file paths
new_sampleannot <- lapply(cells, function(cell) {
  sampleannot$cell_type <- rep(cell, nrow(sampleannot))
  sampleannot$full_sample_path <- paste0(bedDir, "/", cell, "_", sampleannot$fragmentFiles)
  sampleannot$full_sample_path <- gsub(x = sampleannot$full_sample_path, pattern = ".tsv.gz", replacement = ".tsv.gz.bed")
  return(sampleannot)
})
new_sampleannot <- rlist::list.rbind(new_sampleannot)
new_sampleannot <- dplyr::filter(new_sampleannot, full_sample_path %in% all_pseudobulks)
new_sampleannot$chraccr_path <- gsub(x = new_sampleannot$full_sample_path, pattern = paste0(bedDir, "/"), replacement = "")
new_sampleannot <- dplyr::select(new_sampleannot, full_sample_path, chraccr_path, cell_type, sample_exposure_type)

diffCompNames <- character(0)

# Loop through the cells and create comparisons
for (i in 1:length(cells)) {
  for (j in (i + 1):length(cells)) {
    if (!is.na(cells[i]) && !is.na(cells[j]) && cells[i] != cells[j]) {
      comparison <- paste(cells[i], "vs", cells[j], "[cell_type]")
      diffCompNames <- c(diffCompNames, comparison)
    }
  }
}


# lapply(diff, function(diffCompNames) {
# get the cell names
# cells <- unlist(strsplit(diff, " vs "))
# cells <- c(gsub(" \\[.*\\]", "", cells[1]),gsub(" \\[.*\\]", "", cells[2]) )

# reassign the fragment file paths
# sampleannot <- dplyr::filter(new_sampleannot,cell_type %in% cells)

# get the peak set
peaks <- getPeakSet(project)
regionSetList <- list(
  archr_peaks = sort(peaks) # ,
  # tiling200bp = muRtools::getTilingRegions("hg38", width = 200L, onlyMainChrs = TRUE)
)
# set configuration elements
setConfigElement("annotationColumns", c("full_sample_path", "sample_exposure_type", "cell_type"))
setConfigElement("differentialColumns", c("cell_type"))
# setConfigElement("filteringCovgCount", 1L)
setConfigElement("filteringSexChroms", TRUE)
# setConfigElement("filteringCovgReqSamples", 0.005)
setConfigElement("differentialCutoffL2FC", 0.5)
setConfigElement("normalizationMethod", "quantile")
setConfigElement("differentialCompNames", diffCompNames)
setConfigElement("lolaDbPaths", "/icbb/projects/igunduz/annotation/lolaDB/hg38/")

# if the rundir exist continue with existing analysis
if (!file.exists(rundir)) {
  # run ChrAccR on the aggregated fragment files
  ChrAccR::run_atac(
    anaDir = rundir, genome = "hg38",
    input = "full_sample_path", sampleAnnot = dplyr::distinct(new_sampleannot),
    sampleIdCol = "full_sample_path", regionSets = regionSetList
  )
} else {
  # run ChrAccR on the aggregated fragment files
  ChrAccR::run_atac(
    anaDir = rundir, genome = "hg38",
    sampleAnnot = dplyr::distinct(new_sampleannot),
    sampleIdCol = "full_sample_path", regionSets = regionSetList
  )
}
# })


#####################################################################
