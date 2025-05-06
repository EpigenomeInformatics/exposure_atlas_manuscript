#!/usr/bin/env Rscript

#####################################################################
# 06_C19_mTFR.R
# created on 12-11-24 by Irem Gunduz
# Run methylTFR for covid pseudobulks
#####################################################################

set.seed(42)
suppressPackageStartupMessages({
  library(data.table)
  library(methylTFRAnnotationHg38)
  library(dplyr)
  library(muLogR)
  library(methylTFR)
  library(GenomicRanges)
})

# extract unique cells
cells <- data.table::fread("/icbb/projects/igunduz/DARPA/allc_sample_annot_final.csv") %>%
  dplyr::select(!V1) %>%
  dplyr::filter(!cell_type %in% c("Other-cell", "Tc-Eff", "Th-Eff"))

# get the unique cell types
cells <- names(table(cells$cell_type))
cells <- sort(cells)[2]

shareddir <- "/icbb/projects/igunduz/methylTFR_081124/covid_mtfr_061224_c1/"
if (!dir.exists(shareddir)) {
  dir.create(shareddir)
}
comparasions <- c("C19_ctrl", "C19_sev")
logger::log_info("Cell types: ", paste0(cells, collapse = ", "))
logger::log_info("Condition: ", paste0(comparasions, collapse = ", "))
logger::log_info("Out dir: ", shareddir)
logger::log_info("Loading the TF binding sites, GC freqs and GC dist")
motifSetList <- c("jaspar2020_distal", "altius")[1]


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

logger::log_info("Start running methylTFR...")
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
  gcfreqs <- readRDS("/icbb/projects/igunduz/annotation/jaspar2020_distal_motif_gcfreq_vs2.rds")
  tf_bindsites <- getTFbindsites(motifSet = tfset)
  logger::log_info("Number of motifs in gcfreqs: ", length(gcfreqs))
  logger::log_info("Start running methylTFR...")
  gc_dist <- readRDS("/icbb/projects/igunduz/annotation/genome_wide_GC.RDS")
  for (comp in comparasions) {
    for (cell in cells) {
      if (!file.exists(paste0(out.dir, cell, "_", comp, "_", motifset, "_deviations.RDS"))) {
        afile <- read.table(paste0("/icbb/projects/igunduz/DARPA_analysis/bedToEPP/EPP_271023/", cell, "/sample_annotation_unique.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
        afile <- dplyr::filter(afile, condition == comp)
        write.table(afile, file = paste0("/icbb/projects/igunduz/DARPA_analysis/bedToEPP/EPP_271023/", cell, "/sannot_subset.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
        logger::log_info("Running methylTFR for ", comp, " in ", cell)

        # deviation score matrix
        deviations <- run_methyltfr(
          annfile = paste0("/icbb/projects/igunduz/DARPA_analysis/bedToEPP/EPP_271023/", cell, "/sannot_subset.tsv"),
          full_path = TRUE,
          threads = 30,
          sampleColName = "bedFile",
          chunkSize = 20,
          tf_bindsites = tf_bindsites,
          gcfreqs = gcfreqs,
          gc_dist = gc_dist,
          filetype = "EPP",
          ignoreStrand = TRUE,
          enhancer = distal,
          cov_threshold = 1
        )
        # save deviations as RDS
        saveRDS(deviations, paste0(out.dir, cell, "_", comp, "_", motifset, "_deviations.RDS"))
        rm(deviations)
        ChrAccR:::cleanMem()
        logger::log_info("Finished running methylTFR for ", cell, " in ", comp)
      }
    }
  }
}

#####################################################################
