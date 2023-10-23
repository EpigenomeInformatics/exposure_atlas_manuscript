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
cells <- data.table::fread("/icbb/projects/igunduz/DARPA/allc_sample_annot_final.csv") %>%
  dplyr::select(!V1) %>%
  dplyr::filter(!cell_type %in% c("Other-cell", "Tc-Eff", "Th-Eff"))

# get the unique cell types
cells <- names(table(cells$cell_type))
cells <- sort(cells)

# data.dir <- "/icbb/projects/igunduz/DARPA/data/pseudoBulks/perSample"
out.dir <- "/icbb/projects/igunduz/DARPA_analysis/methyltfr_041023/ClustResults/"
if (!dir.exists(out.dir)) {
  dir.create(out.dir)
}
# comparasions <- c("C19_mild","C19_ctrl","C19_sev","HIV_acu","HIV_ctrl","HIV_chr","OP_low","OP_med","OP_high","Influenza_ctrl","Influenza_d30")
comparasions <- c("HIV_acu", "HIV_ctrl", "HIV_chr", "OP_low", "OP_med", "OP_high")

logger::log_info("Cell types: ", paste0(cells, collapse = ", "))
logger::log_info("Condition: ", paste0(comparasions, collapse = ", "))
logger::log_info("Out dir: ", out.dir)
logger::log_info("Loading the TF binding sites, GC freqs and GC dist")

gcfreqs <- getGCfreq(motifSet = "altius")
gc_dist <- getGenomeGC()
tf_bindsites <- getTFbindsites(motifSet = "altius")

logger::log_info("Number of motifs in gcfreqs: ", length(gcfreqs))
if (length(gcfreqs) != length(tf_bindsites)) {
  logger::log_info("Number of motifs in gcfreqs and tf_bindsites are not equal")
  tf_bindsites <- tf_bindsites[names(gcfreqs)]
}
logger::log_info("Number of NAs in gc_dist: ", sum(is.na(gc_dist)))


logger::log_info("Start running methylTFR...")
for (comp in comparasions) {
  for (cell in cells) {
    if (!file.exists(paste0(out.dir, cell, "_", comp, "_deviations.RDS"))) {
      afile <- read.table(paste0("/icbb/projects/igunduz/DARPA_analysis/bedToEPP/EPP/", cell, "/sample_annotation_unique.tsv"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
      afile <- dplyr::filter(afile, condition == comp)
      write.table(afile, file = paste0("/icbb/projects/igunduz/DARPA_analysis/bedToEPP/EPP/", cell, "/sannot_subset.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
      logger::log_info("Running methylTFR for ", comp, " in ", cell)
      # deviation score matrix
      deviations <- run_methyltfr(
        annfile = paste0("/icbb/projects/igunduz/DARPA_analysis/bedToEPP/EPP/", cell, "/sannot_subset.tsv"),
        full_path = TRUE,
        threads = 30,
        chunkSize = 10,
        tf_bindsites = tf_bindsites,
        gcfreqs = gcfreqs,
        gc_dist = gc_dist,
        filetype = "EPP"
      )
      # save deviations as RDS
      saveRDS(deviations, paste0(out.dir, cell, "_", comp, "_deviations.RDS"))
      rm(deviations)
      ChrAccR:::cleanMem()
      logger::log_info("Finished running methylTFR for ", cell, " in ", comp)
    }
  }
}


#####################################################################
