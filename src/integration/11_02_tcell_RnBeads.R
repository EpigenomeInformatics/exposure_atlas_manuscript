#!/usr/bin/env Rscript

#####################################################################
# 08_tcell_RnBeads.R
# created on 26-05-25 by Irem Gunduz
# Run RnBeads vanilla analysis using pseudobulks
#####################################################################

set.seed(42)
suppressPackageStartupMessages({
  library(RnBeads)
  library(dplyr)
  library(grid)
})
# extract unique cells
cells <- data.table::fread("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/sample_annots/allc_sample_annot_final.csv") %>%
  dplyr::select(!V1) %>%
  dplyr::filter(!cell_type %in% c("Other-cell", "Tc-Eff", "Th-Eff"))

# get the unique cell types
cells <- names(table(cells$cell_type))
cells <- sort(cells)[4:7]

# Directory where your data is located
data.dir <- "/icbb/projects/igunduz/DARPA/data/pseudoBulks/perSample"

# Directory where the output should be written to
analysis.dir <- "/icbb/projects/igunduz/DARPA_analysis/rnbeads_Tcell_run_260525/"

if (!dir.exists(analysis.dir)) {
  dir.create(analysis.dir)
}

comparasions <- c(
  "Tc-Naive_vs_Tc-Mem", "Th-Naive_vs_Th-Mem"
)

for (comp in comparasions[2]) {
  logger.info(paste0("Running RnBeads for ", comp))
  dir.create(paste0(analysis.dir, "/", comp))
  analysis.dir <- paste0(analysis.dir, "/", comp)
  sample.annotation <- paste0(data.dir, "/", comp, ".tsv")

  # Directory where the report files should be written to
  report.dir <- file.path(analysis.dir, "reports")
  rnb.initialize.reports(report.dir)

  # Multiprocess
  num.cores <- 30
  parallel.setup(num.cores)

  # set the project specific options
  rnb.options(
    analysis.name = paste0(comp),
    assembly = "hg38",
    import.table.separator = "\t",
    region.aggregation = "sum",
    import.default.data.type = "data.dir",
    import.bed.style = "bismarkCytosine",
    filtering.sex.chromosomes.removal = TRUE,
    disk.dump.big.matrices = TRUE,
    strand.specific = FALSE,
    enforce.memory.management = TRUE,
    differential.comparison.columns = "cell_type",
    differential.enrichment.lola = FALSE, # this is not working
    identifiers.column = "bedFile",
    region.types = c("promoters", "genes")
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
