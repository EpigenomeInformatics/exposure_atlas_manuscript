#!/usr/bin/env Rscript

#####################################################################
# 04_run_methylTFR_clust.R
# created on 2023-08-24 by Irem Gunduz
# Run methylTFR clustered motifs analysis using pseudobulks
#####################################################################

set.seed(42)
logger::log_info("Loading libraries...")
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(methylTFRAnnotationHg38))
suppressPackageStartupMessages(library(methylTFR))

# extract unique cells
cells <- data.table::fread("/icbb/projects/igunduz/DARPA_analysis/methyltfr_041023/allc_shared_samples_ClustResults_121123/sample_annotation_unique.tsv") # %>%

# get the unique cell types
cells <- names(table(cells$cell_type))
cells <- sort(cells)
motifSetList <- c("jaspar2020", "altius")


shareddir <- "/icbb/projects/igunduz/DARPA_analysis/mtfr_results_241123/shared/updated/distal/"
if (!dir.exists(shareddir)) {
  dir.create(shareddir)
}


for (motifset in motifSetList) {
  out.dir <- paste0(shareddir, paste0(motifset, "/"))
  if (!dir.exists(out.dir)) {
    dir.create(out.dir)
  }

  logger::log_info("Cell types: ", paste0(cells, collapse = ", "))
  logger::log_info("Out dir: ", out.dir)
  logger::log_info("Loading the TF binding sites, GC freqs and GC dist")

  gcfreqs <- getGCfreq(motifSet = motifset)
  gc_dist <- getGenomeGC()
  tf_bindsites <- getTFbindsites(motifSet = motifset)

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
        threads = 30,
        chunkSize = 20,
        tf_bindsites = tf_bindsites,
        gcfreqs = gcfreqs,
        gc_dist = gc_dist,
        filetype = "EPP",
        ignoreStrand = TRUE,
        enhancer = NULL
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
