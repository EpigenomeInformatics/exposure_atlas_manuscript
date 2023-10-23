#!/usr/bin/env Rscript

#####################################################################
# 02_run_RnBeads.R
# created on 2023-08-24 by Irem Gunduz
# Run RnBeads vanilla analysis using pseudobulks
#####################################################################

set.seed(42)
suppressPackageStartupMessages(library(RnBeads))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(grid))

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
# analysis.dir <- "/icbb/projects/igunduz/DARPA/Generated/RnBeadsRuns/"
analysis.dir <- "/icbb/projects/igunduz/DARPA_analysis/rnbeads_run_061023/"

if (!dir.exists(analysis.dir)) {
  dir.create(analysis.dir)
}

comparasions <- c(
  "C19_mild_vs_Ctrl", "C19_sev_vs_Ctrl",
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

      ## Data export
      # rnb.options(export.to.csv = FALSE, export.to.bed = TRUE)
      # rnb.run.tnt(rnb.set, report.dir)
      # rnb.export.to.trackhub(rnb.set, out.dir, reg.type = "sites", data.type = "bigBed")
      # rnb.execute.tnt(rnb.set,"/icbb/projects/igunduz/scmeth/bedtoEPP_300923/EPP/",TRUE,region.types='sites')
      ## Exploratory analysis
      # rnb.run.exploratory(rnb.set, report.dir)

      ## Differential methylation
      rnb.run.differential(rnb.set, report.dir)
    }
  }
}

logger.start("Running LOLA for Covid-19 in Monocytes")
analysis.dir <- "/icbb/projects/igunduz/DARPA_analysis/rnbeads_run_061023/Monocyte/C19_sev_vs_Ctrl/"
cell <- cells[2] # Monocyte

logger.info("Loading RnBeads objects")
diffMeth <- load.rnb.diffmeth(paste0(analysis.dir, "/reports/differential_methylation_data/differential_rnbDiffMeth/"))
rnb_set <- load.rnb.set(paste0(analysis.dir, "/reports/data_import_data/rnb.set_preprocessed/"))

# logger.info("Loading LOLA database")
lolaDb_path <- "/icbb/projects/igunduz/annotation/lolaDB/hg38/"

# Run LOLA
res <- performLolaEnrichment.diffMeth(rnb_set, diffMeth, lolaDb_path)
logger.info("Saving results")
saveRDS(res, paste0(analysis.dir, "/reports/differential_methylation_data/differential_rnbDiffMeth/lola_results.rds"))
source("/icbb/projects/igunduz/exposure_atlas_manuscript/src/utils/plot_utils.R")

# lola meth
outputDir <- "/icbb/projects/igunduz/DARPA_analysis/rnbeads_run_061023/Monocyte/C19_sev_vs_Ctrl/reports/differential_methylation_data/differential_rnbDiffMeth/lola_results.rds"
lolaDb <- loadLolaDbs(lolaDb_path)
p <- lolaVolcanoPlotC19("Mono_CD14", lolaDb, outputDir, pValCut = 1.5, region = "archrPeaks", database = "TF_motif_clusters", signifCol = "qValue", cnames = c("loss", "gain"), n = 3, colorpanel = c())$plot
ggsave(plot = p, "/icbb/projects/igunduz/Figures/TF_motif_clusters_METH.pdf", width = 20, height = 20)

p <- lolaVolcanoPlotC19("Mono_CD14", lolaDb, outputDir, pValCut = 1.5, region = "archrPeaks", database = "TF_motifs_vert", signifCol = "qValue", cnames = c("loss", "gain"), n = 3, colorpanel = c())$plot
ggsave(plot = p, "/icbb/projects/igunduz/Figures/TF_motif_vert_METH.pdf", width = 20, height = 20)

p <- lolaVolcanoPlotC19("Mono_CD14", lolaDb, outputDir, pValCut = 1.5, region = "archrPeaks", database = "TF_motifs", signifCol = "qValue", cnames = c("loss", "gain"), n = 3, colorpanel = c())$plot
ggsave(plot = pp, "/icbb/projects/igunduz/Figures/TF_motifs_METH.pdf", width = 20, height = 20)
