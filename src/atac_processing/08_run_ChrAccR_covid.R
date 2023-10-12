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
outputDir <- "/icbb/projects/igunduz/archr_project_011023"
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)

# read the sample annotation
sampleannot <- read.delim("/icbb/projects/igunduz/sampleannot.tsv") %>%
  dplyr::filter(sample_exposure_group %in% c("C19_mild", "C19_mod", "C19_sev", "C19_ctrl"))
sampleannot$fragmentFiles <- gsub(x = sampleannot$fragmentFiles, pattern = ".bed", replacement = ".tsv.gz")

# read the batch info
batch <- data.table::fread("/icbb/projects/igunduz/DARPA/ATAC_metadata_covid.csv") %>%
  dplyr::mutate(fragmentFiles = paste0(arrow_name, "_fragments.tsv.gz")) %>%
  dplyr::mutate(processing_date = as.factor(processing_date))

# add batch info to sampleannot
sampleannot <- merge(sampleannot, batch, by = "fragmentFiles")

# set directory for the output
outputDir <- "/icbb/projects/igunduz/DARPA_analysis/BedFiles_final/"
# rundir <- paste0("/icbb/projects/igunduz/DARPA/Generated/ChrAccRuns_covid_030723/")
rundir <- paste0("/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_", Sys.Date(), "/")
if(!dir.exists(rundir)) dir.create(rundir)

# set the cell types and the comparisons
#cells <- c("B_mem", "B_naive", "Mono_CD14","Mono_CD16", "NK_CD16", "T_mem_CD8", "T_mem_CD4", "T_naive" "T_mix","T_naive")
cells <- c("Mono_CD14", "NK_CD16", "T_mem_CD8", "T_mem_CD4","T_naive")


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
  # peaks <- peaks[grep(paste0("^",cell), peaks$GroupReplicate), ]
  regionSetList <- list(
    archr_peaks = sort(peaks)#,
    #tiling200bp = muRtools::getTilingRegions("hg38", width = 200L, onlyMainChrs = TRUE)
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
