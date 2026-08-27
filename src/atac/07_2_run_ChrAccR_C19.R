#!/usr/bin/env Rscript

#####################################################################
# 08_run_ChrAccR_covid.R
# created on 2023-08-24 by Irem Gunduz
# Run vanilla ChrAccR analysis considering batch effects for Covid samples
#####################################################################

suppressPackageStartupMessages({
  library(ArchR)
  library(ChrAccR)
  library(dplyr)
})

set.seed(12) # set seed
outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023"
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)

# read the sample annotation
sampleannot <- read.delim("/icbb/projects/igunduz/sampleannot.tsv") %>%
  dplyr::filter(sample_exposure_group %in% c("C19_mild", "C19_mod", "C19_sev", "C19_ctrl"))
sampleannot$fragmentFiles <- gsub(x = sampleannot$fragmentFiles, pattern = ".bed", replacement = ".tsv.gz")

# read the batch info
batch <- data.table::fread("/scratch/icbb/igunduz/irem_github/exposure_atlas_manuscript/sample_annots/ATAC_metadata_covid.csv") %>%
  dplyr::mutate(fragmentFiles = paste0(arrow_name, "_fragments.tsv.gz")) %>%
  dplyr::mutate(processing_date = as.factor(processing_date))

# add batch info to sampleannot
sampleannot <- merge(sampleannot, batch, by = "fragmentFiles")

# set directory for the output
outputDir <- "/icbb/projects/igunduz/DARPA_analysis/BedFiles_final/"
#rundir <- paste0("/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_", Sys.Date())
rundir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/"
if (!dir.exists(rundir)) dir.create(rundir)

# set the cell types and the comparisons
cells <- c("Mono_CD16", "Mono_CD14", "NK_CD16", "T_mem_CD8", "T_mem_CD4", "T_naive")


diffCompNames <- c(
  "C19_mild vs C19_ctrl [sample_exposure_group]",
  "C19_mod vs C19_ctrl [sample_exposure_group]",
  "C19_sev vs C19_ctrl [sample_exposure_group]"
)

lapply(cells, function(cell) {
  # reassign the fragment file paths
  sampleannot$fragmentFiles2 <- paste0(
    paste0(outputDir, cell, "/"),
    gsub(x = sampleannot$fragmentFiles, pattern = ".tsv.gz", replacement = ".tsv.gz.bed")
  )
  # check the bed files
  beds <- Sys.glob(file.path(paste0(outputDir, cell), "*.bed"))
  beds <- beds[beds %in% sampleannot$fragmentFiles2]
  sampleannot <- sampleannot %>%
    dplyr::filter(fragmentFiles2 %in% beds) %>%
    dplyr::mutate(beds = beds)
  sampleannot$beds <- beds


  # Filter the original table
  filtered_table <- table(sampleannot$sample_exposure_group)
  filtered_table <- filtered_table[filtered_table > 2]
  filtered_names <- names(filtered_table)

  # filter the sample annotation
  sampleannot <- sampleannot[sampleannot$sample_exposure_group %in% filtered_names, ]

  diffCompNames <- diffCompNames[sapply(strsplit(diffCompNames, " vs "), function(x) {
    before <- gsub(" \\[.*\\]", "", x[1]) # Remove text within square brackets
    after <- gsub(" \\[.*\\]", "", x[2]) # Remove text within square brackets
    before_in_filtered <- before %in% filtered_names
    after_in_filtered <- after %in% filtered_names
    before_in_filtered && after_in_filtered
  })]

  # get peaks from the ArchR project and subset based on cell type
  peaks <- getPeakSet(project)
  regionSetList <- list(
    archr_peaks = sort(peaks)
  )

  # set configuration elements
  setConfigElement("differentialAdjColumns", "processing_date")
  setConfigElement("differentialColumns", c("sample_exposure_group"))
  setConfigElement("annotationColumns", c("sampleIdCol", "sample_exposure_type", "sample_exposure_group"))
  setConfigElement("filteringSexChroms", TRUE)
  setConfigElement("differentialCutoffL2FC", 0.5)
  setConfigElement("normalizationMethod", "quantile")
  setConfigElement("differentialCompNames", diffCompNames)
  setConfigElement("lolaDbPaths", "/icbb/projects/igunduz/annotation/lolaDB/hg38/")

  message("Running vanilla analysis for ", paste0(cell))
  # if the rundir exist continue with existing analysis
  if (!file.exists(paste0(rundir, cell))) {
    # run ChrAccR on the aggregated fragment files
    ChrAccR::run_atac(
      anaDir = paste0(rundir, cell), genome = "hg38",
      input = "beds", sampleAnnot = sampleannot,
      sampleIdCol = "sampleIdCol", regionSets = regionSetList
    )
  } else {
    # run ChrAccR on the aggregated fragment files
    ChrAccR::run_atac(
      anaDir = paste0(rundir, cell), genome = "hg38",
      sampleAnnot = sampleannot,
      sampleIdCol = "sampleIdCol", regionSets = regionSetList
    )
  }
})


#####################################################################
# R3.4: direct moderate-vs-severe differential accessibility
# Loads the EXISTING processed ChrAccR DsATAC object (no re-run of import/
# filtering/normalization) and runs ChrAccR's own differential module for the
# C19_sev vs C19_mod comparison, so it uses the same normalization, model and
# cutoffs as the control comparisons above.
#####################################################################
cell_ms <- "Mono_CD14"
# existing ChrAccR run directory (see src/atac/08_1_chraccr_plots.R)
existing_run <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/"
plot_dir <- "/icbb/projects/igunduz/finalize_echo_050824/C19/"

# load the processed object
ds <- ChrAccR::loadDsAcc(paste0(existing_run, cell_ms, "/data/dsATAC_filtered/"))

# configure the differential module for the direct moderate-vs-severe test,
# mirroring the settings used for the control comparisons above
setConfigElement("differentialColumns", c("sample_exposure_group"))
setConfigElement("differentialCompNames", c("C19_sev vs C19_mod [sample_exposure_group]"))
setConfigElement("differentialAdjColumns", "processing_date")
setConfigElement("differentialCutoffL2FC", 0.5)
setConfigElement("filteringSexChroms", TRUE)

# run into a FRESH analysis dir: run_atac_differential skips if a differential
# report already exists in anaDir, so we do not reuse the control-comparison dir
ms_dir <- paste0(existing_run, cell_ms, "_mod_vs_sev/")
if (!dir.exists(ms_dir)) dir.create(ms_dir, recursive = TRUE)
run_atac_differential(ds, ms_dir)

# read the resulting DAR table. We search for the diffTab file recursively so we can find it even if the subdir structure changes in future versions of ChrAccR.
source("/scratch/icbb/igunduz/irem_github/exposure_atlas_manuscript/utils/helpers.R") # cutL0.5fc2Padj05

diff_files <- list.files(ms_dir, pattern = "diffTab.*\\.tsv$", recursive = TRUE, full.names = TRUE)
message("diffTab files found under ", ms_dir, ": ", paste(basename(diff_files), collapse = ", "))
stopifnot(length(diff_files) >= 1)

dm <- read.delim(diff_files[1])
isDiff <- cutL0.5fc2Padj05(dm[, c("log2FoldChange", "padj")])
isDiff[is.na(isDiff)] <- FALSE
mod_sev_tab <- data.frame(
  Chromosome = dm$chrom,
  Start = dm$chromStart + 1, # BED 0-based -> 1-based, matches get_dar() ids
  End = dm$chromEnd,
  log2FC = dm$log2FoldChange,
  padj = dm$padj,
  isDiff = isDiff,
  stringsAsFactors = FALSE
)
mod_sev_tab$id <- paste0(mod_sev_tab$Chromosome, ":", mod_sev_tab$Start, "-", mod_sev_tab$End)

message("moderate-vs-severe DARs (|log2FC|>0.5, padj<0.05): ", sum(mod_sev_tab$isDiff))
write.table(mod_sev_tab,
  file = paste0(plot_dir, "diffTab_mod_vs_sev_archrPeaks.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
saveRDS(mod_sev_tab, paste0(plot_dir, "mod_vs_sev_DAR.rds"))

#####################################################################
