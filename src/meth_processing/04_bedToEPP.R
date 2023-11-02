suppressPackageStartupMessages({
library(dplyr)
library(RnBeads)
library(data.table)
library(muLogR)
})

# extract unique cells
cells <- data.table::fread("/icbb/projects/igunduz/DARPA/allc_sample_annot_final.csv") %>%
  dplyr::select(!V1) %>%
  #dplyr::filter(!cell_type == "B-cell")
  dplyr::filter(!cell_type %in% c("Other-cell","Tc-Eff","Th-Eff"))

# get the unique cell types
cells <- names(table(cells$cell_type))
cells <- sort(cells)[2]

# load the functions
source("/icbb/projects/igunduz/exposure_atlas_manuscript/src/utils/methyltfr_utils.R")

comparasions <- c(
 "C19_sev_vs_Ctrl", "C19_mild_vs_Ctrl", 
  "HIV_acu_vs_Ctrl", "Influenza_ctrl_vs_d30","HIV_chr_vs_Ctrl", 
  "OP_high_vs_low","OP_high_vs_med","OP_low_vs_med"
)[1]

logger.info("Setting up directories...")
rnbeads.dir <- "/icbb/projects/igunduz/DARPA_analysis/rnbeads_run_061023/"
out.dir <- "/icbb/projects/igunduz/DARPA_analysis/bedToEPP_new/"
if (!dir.exists(out.dir)) {dir.create(out.dir)}
#out.dir <- "/icbb/projects/igunduz/DARPA_analysis/bedToEPP/EPP/"
#if (!dir.exists(out.dir)) {dir.create(out.dir)}
data.dir <- "/icbb/projects/igunduz/DARPA/data/pseudoBulks/perSample"

for (cell in cells) {
  for (comp in comparasions) {
    if (!dir.exists(paste0(out.dir , cell))) {
      dir.create(paste0(out.dir , cell))
    }
    if (!dir.exists(paste0(out.dir , cell, "/", comp))) {
      dir.create(paste0(out.dir , cell, "/", comp))
    }
    # Directory where the output should be written to
    analysis.dir <- paste0(rnbeads.dir, cell, "/", comp)
    sample.annotation <- paste0(data.dir, "/", cell, "_", comp, "_sampleannot.tsv")
    logger.info(c("Checking if all EPP files exists for ",cell," - ",comp))
    sannot <- data.table::fread(sample.annotation) %>%
      dplyr::mutate(files = paste0(bedFile, ".tsv")) %>%
      dplyr::select(files)

    #check if all EPP files exists  
    if (!all(file.exists(paste0(out.dir ,cell,"/",comp,"/",sannot$files)))) {
      logger.start(c("Creating EPP files for ",cell," - ",comp))
      # Directory where the report files should be written to
      report.dir <- file.path(analysis.dir, "reports")

      rnb.set <- RnBeads::load.rnb.set(paste0(report.dir, "/data_import_data/rnb.set_preprocessed"))
      result <- RnBeads::rnb.RnBSet.to.GRangesList(rnb.set, "sites")
      rm(rnb.set)
      # convert sites to EPP object
      epps <- converToEPP(result, filePath = paste0(out.dir , cell, "/", comp))

      # read sample annotation and fix bedfile names
      sampleannot <- data.table::fread(paste0(data.dir, "/", cell, "_", comp, "_sampleannot.tsv")) %>%
        dplyr::mutate(bedFile = paste0(cell_type, "_", Common_Minimal_Informative_ID, ".bedGraph.tsv"))%>%
        dplyr::arrange(desc(row_number()))

      filePath <- paste0(out.dir , cell, "/", comp)

      # write the sample annotation to folder
      write.table(sampleannot, file = paste0(filePath, "/sample_methylation_summary.tsv"), row.names = FALSE, sep = "\t")
      logger.completed()
      rm(sampleannot,result)
      ChrAccR:::cleanMem()
    }else{
      rm(sannot)
    }
  }
}
