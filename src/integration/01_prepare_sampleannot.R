#!/usr/bin/env Rscript

#####################################################################
# 01_prepare_sampleannot.R
# created on 2023-08-29 written by Irem Gunduz
# Preparing sample annotation for integration
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(ArchR)
  library(dplyr)
  library(muLogR)
  library(muRtools)
})

outputDir <- "/icbb/projects/igunduz/ATAC_processed_final/"
cellAnnot_atac <- read.delim("/icbb/projects/igunduz/artemis_run/final_moon_240823/cellColData.tsv")
rownames(cellAnnot_atac) <- cellAnnot_atac$cellId_archr

colname_samplId_atac <- "sample_sampleId_cminid"
colname_samplId_meth <- "CommonMinID"

cellAnnot_meth <- read.delim("/icbb/projects/igunduz/DARPA/allc_full_sampleannot_290823.tsv")
rownames(cellAnnot_meth) <- cellAnnot_meth[, "Cell_UID"]
cellAnnot_meth <- cellAnnot_meth %>%
  dplyr::filter(sample_exposure_type %in% c("COVID", "FLU", "HIV", "OP"))


logger.start("Finding common samples")
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

# Save data frame as TSV
write.table(cellAnnot_meth, "/icbb/projects/igunduz/DARPA_analysis/artemis_031023/sample_annot.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

logger.completed()

rawdir <- "/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/"
if (!dir.exists(rawdir)) {
  dir.create(rawdir)
}
logger.start("Saving Archr peak regions...")
ap_ha <- ArchR::loadArchRProject("/icbb/projects/igunduz/archr_project_011023/", force = T, showLogo = F)
regGr <- ArchR::getPeakSet(ap_ha)
saveRDS(regGr, paste0(rawdir, "regionsGR.rds"))
logger.completed()
#####################################################################
