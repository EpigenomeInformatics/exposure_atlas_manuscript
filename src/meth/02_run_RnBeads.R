#!/usr/bin/env Rscript

#####################################################################
# 02_run_RnBeads.R
# created on 2023-08-24 by Irem Gunduz
# Run RnBeads vanilla analysis using pseudobulks
#####################################################################

set.seed(42)
suppressPackageStartupMessages({
  library(RnBeads)
  library(dplyr)
  library(grid)
})
# extract unique cells
cells <- data.table::fread("/icbb/projects/igunduz/DARPA/allc_sample_annot_final.csv") %>%
  dplyr::select(!V1) %>%
  dplyr::filter(!cell_type %in% c("Other-cell", "Tc-Eff", "Th-Eff"))

# get the unique cell types
cells <- names(table(cells$cell_type))
cells <- sort(cells)

# Directory where your data is located
data.dir <- "/icbb/projects/igunduz/DARPA/data/pseudoBulks/perSample"

# Directory where the output should be written to
analysis.dir <- "/icbb/projects/igunduz/DARPA_analysis/rnbeads_run_061023/"

if (!dir.exists(analysis.dir)) {
  dir.create(analysis.dir)
}

comparasions <- c(
  "HIV_acu_vs_Ctrl", "Influenza_ctrl_vs_d30", "HIV_chr_vs_Ctrl",
  "OP_high_vs_low", "OP_high_vs_med", "OP_low_vs_med"
)

# add archr peaks dataframe to the annotations
# global_peaklist <- readRDS("/icbb/projects/igunduz/DARPA/peak_tables/global_peaklist.rds") %>%
global_peaklist <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/regionsGR.rds") %>%
  data.table::as.data.table() %>%
  dplyr::select(seqnames, start, end) %>%
  as.data.frame()

colnames(global_peaklist) <- c("Chromosome", "Start", "End")
rnb.set.annotation(type = "archr_peaks", regions = global_peaklist, assembly = "hg38")

for (cell in cells) {
  # create cell-type specific analysis directory
  if (!dir.exists(paste0(analysis.dir, cell))) {
    dir.create(paste0(analysis.dir, cell))
  }
  for (comp in comparasions) {
    logger.info(paste0("Running RnBeads for ", cell, " ", comp))
    if (!dir.exists(paste0(analysis.dir, cell, "/", comp))) {
      dir.create(paste0(analysis.dir, cell, "/", comp))
      analysis.dir <- paste0(analysis.dir, cell, "/", comp)
      sample.annotation <- paste0(data.dir, "/", cell, "_", comp, "_sampleannot.tsv")

      # Directory where the report files should be written to
      report.dir <- file.path(analysis.dir, "reports")
      rnb.initialize.reports(report.dir)

      # Multiprocess
      num.cores <- 30
      parallel.setup(num.cores)

      # set the project specific options
      rnb.options(
        analysis.name = paste0(cell),
        assembly = "hg38",
        import.table.separator = "\t",
        region.aggregation = "sum",
        import.default.data.type = "data.dir",
        import.bed.style = "bismarkCytosine",
        filtering.sex.chromosomes.removal = TRUE,
        disk.dump.big.matrices = TRUE,
        strand.specific = FALSE,
        enforce.memory.management = TRUE,
        differential.comparison.columns = comp,
        differential.enrichment.lola = FALSE, # this is not working
        identifiers.column = "bedFile",
        region.types = c("archr_peaks")
      )

      ## Data import
      data.source <- c(data.dir, sample.annotation, 3)
      result <- rnb.run.import(data.source = data.source, data.type = "bs.bed.dir", dir.reports = report.dir)
      rnb.set <- result$rnb.set

      ## Quality Control
      rnb.run.qc(rnb.set, report.dir)

      ## Preprocessing
      rnb.set <- rnb.run.preprocessing(rnb.set, dir.reports = report.dir)$rnb.set

      # save the object
      save.rnb.set(rnb.set, paste0(report.dir, "/data_import_data/rnb.set_preprocessed"), archive = FALSE)

      ## Differential methylation
      rnb.run.differential(rnb.set, report.dir)
    }
  }
}
