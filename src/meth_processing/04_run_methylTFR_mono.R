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
motifSetList <- c("altius", "jaspar2020")

# extract unique cells
cells <- data.table::fread("/icbb/projects/igunduz/DARPA/allc_sample_annot_final.csv") %>%
  dplyr::select(!V1) %>%
  dplyr::filter(!cell_type %in% c("Other-cell", "Tc-Eff", "Th-Eff"))

# get the unique cell types
cells <- names(table(cells$cell_type))
cells <- sort(cells)[2]

out.dir <- "/icbb/projects/igunduz/DARPA_analysis/mtfr_results_241123/shared/updated/covid_mono/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
out.dir <- "/icbb/projects/igunduz/DARPA_analysis/mtfr_results_241123/shared/updated/covid_mono/distal/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}

comparasions <- c("C19_ctrl", "C19_sev")
logger::log_info("Cell types: ", paste0(cells, collapse = ", "))
logger::log_info("Condition: ", paste0(comparasions, collapse = ", "))
logger::log_info("Out dir: ", out.dir)
logger::log_info("Loading the TF binding sites, GC freqs and GC dist")


if (!file.exists("/icbb/projects/igunduz/annotation/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")) {
  distal <- fread("/icbb/projects/share/annotations/lolaDB/hg38/EnsemblRegBuildBP/regions/regionSet_1.bed", header = FALSE)
  distal$V6 <- str_replace(distal$V6, ".", "*")
  distal <- GRanges(
    seqnames = distal$V1,
    ranges = IRanges(start = distal$V2, end = distal$V3),
    strand = distal$V6
  )
  saveRDS(distal, "/icbb/projects/igunduz/annotation/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
} else {
  distal <- readRDS("/icbb/projects/igunduz/annotation/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
}

for (motifset in motifSetList) {
  gcfreqs <- getGCfreq(motifSet = motifset)
  gc_dist <- getGenomeGC()
  tf_bindsites <- getTFbindsites(motifSet = motifset)

  logger::log_info("Number of NAs in gc_dist: ", sum(is.na(gc_dist)))
  epp_dir <- "/icbb/projects/igunduz/DARPA_analysis/bedToEPP/EPP_271023/"

  logger::log_info("Start running methylTFR...")
  for (cell in cells) {
    for (comp in comparasions) {
      if (!file.exists(paste0(out.dir, cell, "_", comp, "_deviations.RDS"))) {
        afile <- read.table(paste0(epp_dir, cell, "/sample_annotation_unique.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        afile <- dplyr::filter(afile, condition == comp)
        write.table(afile, file = paste0(epp_dir, cell, "/sannot_subset.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
        logger::log_info("Running methylTFR for ", comp, " in ", cell)
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
          enhancer = distal
        )
        # save deviations as RDS
        saveRDS(deviations, paste0(out.dir, cell, "_", comp, "_", motifset, "_deviations.RDS"))
        rm(deviations)
        ChrAccR:::cleanMem()
        logger::log_info("Finished running methylTFR for ", cell, " in ", comp, "for motifset ", motifset)
      }
    }
  }
}


#####################################################################
