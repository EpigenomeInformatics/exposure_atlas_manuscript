#!/usr/bin/env Rscript

#####################################################################
# 07_run_ChrAccR.R
# created on 2023-08-24 by Irem Gunduz
# Run vanilla ChrAccR analysis
#####################################################################

suppressPackageStartupMessages({
  library(ArchR)
  library(ChrAccR)
  library(dplyr)
})
set.seed(12) # set seed
outputDir <- "/icbb/projects/igunduz/archr_project_011023"
project <- ArchR::loadArchRProject(outputDir, force = T)

# read the sample annotation
sampleannot <- read.delim("/icbb/projects/igunduz/sampleannot.tsv") %>%
  dplyr::filter(!sample_exposure_group %in% c("C19_mild", "C19_mod", "C19_sev", "C19_ctrl"))
sampleannot$fragmentFiles <- gsub(x = sampleannot$fragmentFiles, pattern = ".bed", replacement = ".tsv.gz")

# set directory for the output
outputDir <- "/icbb/projects/igunduz/DARPA_analysis/BedFiles_Tcells/"
rundir <- paste0("/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_TCELL_subset", Sys.Date(), "/")
if (!dir.exists(rundir)) dir.create(rundir)
# rundir <- "/icbb/projects/igunduz/DARPA/Generated/ChrAccRuns_280623/"

# set the cell types and the comparisons
# cells <- c("B_mem", "B_naive", "Mono_CD14","Mono_CD16", "NK_CD16", "T_mem_CD8", "T_mem_CD4", "T_naive", "T_mix","T_naive","T_mait")
cells <- c("Terminal_TE","Intermediate_TEx","Early_TEx","T_mem_CD8", "T_mem_CD4", "Activated_CD4T", "CD4T_mem", "CD4T_naive", "CD8T_mem", "CD8T_naive",  "Effector_CD8T",  "Other_T",  "Tfh1", "Tfh2", "Th1", "Th17", "Treg1", "Treg2", "Treg3", "Treg4")
diffCompNames <- c(
  "HIV_ctrl vs HIV_chr [sample_exposure_group]",
  "HIV_ctrl vs HIV_acu [sample_exposure_group]",
  "Influenza_d3 vs Influenza_ctrl [sample_exposure_group]",
  "Influenza_d30 vs Influenza_ctrl [sample_exposure_group]",
  "Influenza_d6 vs Influenza_ctrl [sample_exposure_group]",
  "OP_high vs OP_low [sample_exposure_group]",
  "OP_high vs OP_med [sample_exposure_group]",
  "OP_med vs OP_low [sample_exposure_group]"
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
  filtered_table <- filtered_table[filtered_table >= 3]
  filtered_names <- names(filtered_table)

  # filter the sample annotation
  sampleannot <- sampleannot[sampleannot$sample_exposure_group %in% filtered_names, ]
  # filtered_names <- filtered_names[filtered_names!= "C19_ctrl"]

  # Find the matching elements in diffCompNames
  diffCompNames <- diffCompNames[sapply(strsplit(diffCompNames, " vs "), function(x) {
    before <- gsub(" \\[.*\\]", "", x[1]) # Remove text within square brackets
    after <- gsub(" \\[.*\\]", "", x[2]) # Remove text within square brackets
    before_in_filtered <- before %in% filtered_names
    after_in_filtered <- after %in% filtered_names
    before_in_filtered && after_in_filtered
  })]

  # peaks <- readRDS("/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/rawData/regionsGR.rds")
  # peaks <- peaks[grep(paste0("^",cell), peaks$GroupReplicate), ]
  peaks <- getPeakSet(project)
  regionSetList <- list(
    archr_peaks = sort(peaks) # ,
    # tiling200bp = muRtools::getTilingRegions("hg38", width = 200L, onlyMainChrs = TRUE)
  )
  # set configuration elements
  setConfigElement("annotationColumns", c("sampleIdCol", "sample_exposure_type", "sample_exposure_group"))
  setConfigElement("differentialColumns", c("sample_exposure_group"))
  # setConfigElement("filteringCovgCount", 1L)
  setConfigElement("filteringSexChroms", TRUE)
  # setConfigElement("filteringCovgReqSamples", 0.005)
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
