#!/usr/bin/env Rscript

#####################################################################
# 07_mTFR_run.R
# created on 2023-08-24 by Irem Gunduz
# Run methylTFR clustered motifs analysis using pseudobulks
#####################################################################

set.seed(42)
logger::log_info("Loading libraries...")
suppressPackageStartupMessages({
  library(GenomicRanges)
  library(dplyr)
  library(methylTFRAnnotationHg38)
  library(methylTFR)
})

shareddir <- "/icbb/projects/igunduz/mTFR_bias_fix_v3/all_pseudobulks_181124/"
if (!dir.exists(shareddir)) {
  dir.create(shareddir)
}

# extract unique cells
cells <- data.table::fread("/icbb/projects/igunduz/DARPA_analysis/methyltfr_041023/allc_shared_samples_ClustResults_121123/sample_annotation_unique.tsv") # %>%

# get the unique cell types
cells <- names(table(cells$cell_type))
cells <- sort(cells)

# Remove Other-cells
cells <- cells[!grepl("Other", cells)]
motifSetList <- c("altius", "jaspar2020_distal")[2]

if (!file.exists("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")) {
  distal <- fread("/icbb/projects/share/annotations/lolaDB/hg38/EnsemblRegBuildBP/regions/regionSet_1.bed", header = FALSE)
  distal$V6 <- stringr::str_replace(distal$V6, ".", "*")
  distal <- GRanges(
    seqnames = distal$V1,
    ranges = IRanges(start = distal$V2, end = distal$V3),
    strand = distal$V6
  )
  saveRDS(distal, "/icbb/projects/igunduz/annotation/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
} else {
  distal <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
}

gc_dist <- getGenomeGC()
for (motifset in motifSetList) {
  out.dir <- paste0(shareddir, paste0(motifset, "/"))
  if (!dir.exists(out.dir)) {
    dir.create(out.dir)
  }
  distal <- if (motifset != "jaspar2020_distal") {
    NULL
  } else {
    distal
  }
  tfset <- if (motifset == "jaspar2020_distal") {
    "jaspar2020"
  } else {
    motifset
  }
  # tfset <- if(motifset == "altius_distal"){"altius"} else {motifset}
  logger::log_info("Motif set: ", motifset)
  logger::log_info("Out dir: ", out.dir)
  logger::log_info("Loading the TF binding sites, GC freqs and GC dist")
  gcfreqs <- getGCfreq(motifSet = tfset)
  tf_bindsites <- getTFbindsites(motifSet = tfset)
  logger::log_info("Start running methylTFR...")

  logger::log_info("Number of motifs in gcfreqs: ", length(gcfreqs))
  epp_dir <- "/icbb/projects/igunduz/DARPA_analysis/bedToEPP/EPP_shared_samples_201123/"

  logger::log_info("Start running methylTFR...")
  for (cell in cells) {
    if (!file.exists(paste0(out.dir, cell, "_deviations.RDS"))) {
      afile <- read.table(paste0(epp_dir, cell, "/sample_methylation_summary.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      afile <- dplyr::filter(afile)
      afile$bedFile <- paste0(epp_dir, cell, "/", afile$bedFile)
      write.table(afile, file = paste0(epp_dir, cell, "/sannot_subset.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
      logger::log_info("Running methylTFR for in ", cell)
      # deviation score matrix
      deviations <- run_methyltfr(
        annfile = paste0(epp_dir, cell, "/sannot_subset.tsv"),
        full_path = TRUE,
        threads = 31,
        chunkSize = 10,
        tf_bindsites = tf_bindsites,
        gcfreqs = gcfreqs,
        gc_dist = gc_dist,
        filetype = "EPP",
        ignoreStrand = TRUE,
        enhancer = distal,
        sampleColName = "bedFile",
        cov_threshold = 1
      )
      # save deviations as RDS
      saveRDS(deviations, paste0(out.dir, cell, "_deviations.RDS"))
      rm(deviations)
      ChrAccR:::cleanMem()
      logger::log_info("Finished running methylTFR for ", cell)
    }
  }
}

#####################################################################
