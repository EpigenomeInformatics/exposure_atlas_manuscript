#!/usr/bin/env Rscript

#####################################################################
# 02_run_RnBeads.R
# created on 08-04-25 by Irem Gunduz
# Run RnBeads vanilla analysis using pseudobulks for Mono_1, Mono_2 covid
#####################################################################
# Load the functions
set.seed(42)
suppressPackageStartupMessages({
  library(RnBeads)
  library(dplyr)
  library(grid)
  library(LOLA)
})

# extract unique cells
cells <- data.table::fread("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/sample_annots/allc_sample_annot_final.csv") %>%
  dplyr::select(!V1) %>%
  dplyr::filter(!cell_type %in% c("Other-cell", "Tc-Eff", "Th-Eff"))

# get the unique cell types
cells <- names(table(cells$cell_type))
cells <- sort(cells)[2]

# Directory where your data is located
data.dir <- "/icbb/projects/igunduz/DARPA/data/pseudoBulks/perSample"

# Directory where the output should be written to
# analysis.dir <- "/icbb/projects/igunduz/DARPA/Generated/RnBeadsRuns/"
analysis.dir <- "/icbb/projects/igunduz/DARPA_analysis/rnbeads_C19_monosub_080425_V2/"
if (!dir.exists(analysis.dir)) {
  dir.create(analysis.dir)
}

comparasions <- c(
  "C19_sev_vs_Ctrl", "C19_mild_vs_Ctrl"
)
mono_1 <- c(
  "CoV_S_S15_D1",
  "CoV_S_S11_D3",
  "CoV_S_S7_D1",
  "CoV_S_S11_D1"
)
# Create two groups: severe_mono1 and severe_mono2

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

      # Re-organize annotation to create sub monocytes
      sample.annotation <- read.delim(sample.annotation)

      # Update C19_sev_vs_Ctrl column based on mono_1 list and condition
      sample.annotation$C19_sev_vs_Ctrl <- ifelse(sample.annotation$condition == "CommercialControl_healthy",
        "C19_ctrl",
        ifelse(sample.annotation$Common_Minimal_Informative_ID %in% mono_1,
          "C19_sev_mono1",
          "C19_sev_mono2"
        )
      )

      # Save the new sample annotation
      write.table(sample.annotation,
        file = paste0(data.dir, "/", cell, "_", comp, "_sampleannot2.tsv"),
        sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE
      )

      sample.annotation <- paste0(data.dir, "/", cell, "_", comp, "_sampleannot2.tsv")

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
        differential.comparison.columns.all.pairwise = comp,
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

#####################################################################
# LOLA Analysis
#######################################################################

logger.start("Running LOLA for Covid-19 in Monocytes")
analysis.dir <- "/icbb/projects/igunduz/DARPA_analysis/rnbeads_C19_monosub_080425_V2/Monocyte/C19_sev_vs_Ctrl/"
cell <- cells[2] # Monocyte

# Load the LOLA database
lolaDb_path <- "/icbb/projects/share/annotation/lolaDB/hg38/"

logger.info("Loading RnBeads objects")
diffMeth <- load.rnb.diffmeth(paste0(analysis.dir, "/reports/differential_methylation_data/differential_rnbDiffMeth/"))
rnb_set <- load.rnb.set(paste0(analysis.dir, "/reports/data_import_data/rnb.set_preprocessed/"))

# Run LOLA
res <- performLolaEnrichment.diffMeth(rnb_set, diffMeth, lolaDb_path)
logger.info("Saving results")
saveRDS(res, paste0(analysis.dir, "/reports/differential_methylation_data/differential_rnbDiffMeth/lola_results.rds"))
logger.completed()

#####################################################################
