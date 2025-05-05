#!/usr/bin/env Rscript

#####################################################################
# 01_meth_pseudobulks.R
# created on 2023-08-24 by Irem Gunduz
# Create pseudobulks for the methylation data
#####################################################################
suppressPackageStartupMessages({
  library(dplyr)
})
set.seed(12)
source("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/utils/createPseudoBulks.R")

# read the sample annotation
sampleAnnot <- data.table::fread("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/sample_annots/allc_sample_annot_final.csv") %>%
  dplyr::select(!V1) %>%
  dplyr::filter(!cell_type == "Other-cell")

perSample <- "DARPA/data/pseudoBulks/perSample"
if (!dir.exists(perSample)) {
  dir.create(perSample)
}
outputDir <- "/icbb/projects/igunduz/DARPA/data/pseudoBulks/perSample"

# Create single base-pair pseudobulks
createPseudoBulks(sampleAnnot,
  filePathCol = "allC_FilePathfull",
  sampleIdCol = "Common_Minimal_Informative_ID",
  numThreads = 6, mcType = "CGN", groupName = "cell_type", singleBP = TRUE,
  groupName2 = "Common_Minimal_Informative_ID", fileType = "allc", singleCovOff = 99999,
  indexed = FALSE, excludeChr = c("chrX", "chrY", "chrM", "chrL"), outputDir = outputDir
)


# get the unique cell types
cells <- names(table(sampleAnnot$cell_type))

for (cell in cells) {
  # prepare cell type specific sample annotation for RnBeads analysis
  sampleann <- sampleAnnot %>%
    dplyr::filter(cell_type == cell) %>%
    dplyr::select(Common_Minimal_Informative_ID, condition, cell_type) %>%
    distinct(.keep_all = TRUE) %>%
    dplyr::mutate(bedFile = paste0(cell_type, "_", Common_Minimal_Informative_ID, ".bedGraph"))

  # arrange condition groups for differential analysis
  sampleann$C19_mild_vs_Ctrl <- ifelse(sampleann$condition == "CommercialControl_healthy", "C19_ctrl",
    ifelse(sampleann$condition == "COVID_mild", "C19_mild", NA)
  )
  sampleann$C19_sev_vs_Ctrl <- ifelse(sampleann$condition == "CommercialControl_healthy", "C19_ctrl",
    ifelse(sampleann$condition == "COVID_severe", "C19_sev", NA)
  )
  sampleann$HIV_acu_vs_Ctrl <- ifelse(sampleann$condition == "HIV_pre", "HIV_ctrl",
    ifelse(sampleann$condition == "HIV_acute", "HIV_acu", NA)
  )
  sampleann$HIV_chr_vs_Ctrl <- ifelse(sampleann$condition == "HIV_pre", "HIV_ctrl",
    ifelse(sampleann$condition == "HIV_chronic", "HIV_chr", NA)
  )
  sampleann$OP_high_vs_low <- ifelse(sampleann$condition == "OP_high", "OP_high",
    ifelse(sampleann$condition == "OP_low", "OP_low", NA)
  )
  sampleann$OP_high_vs_med <- ifelse(sampleann$condition == "OP_high", "OP_high",
    ifelse(sampleann$condition == "OP_medium", "OP_med", NA)
  )
  sampleann$OP_low_vs_med <- ifelse(sampleann$condition == "OP_low", "OP_low",
    ifelse(sampleann$condition == "OP_medium", "OP_med", NA)
  )
  sampleann$Influenza_ctrl_vs_d30 <- ifelse(sampleann$condition == "FLU_healthy", "Influenza_ctrl",
    ifelse(sampleann$condition == "FLU_day30", "Influenza_d30", NA)
  )

  # check if all files exists
  logger::log_info("Checking if all input files exists.")
  if (!all(file.exists(paste0(outputDir, "/", sampleann$bedFile)))) {
    missingSamples <- sampleann$Common_Minimal_Informative_ID[!file.exists(paste0(outputDir, "/", sampleann$bedFile))]
    stop("Missing input files for samples: ", paste(missingSamples, collapse = ", "))
  }
  logger::log_success("Located all of the samples.")
  # write sampleannot as tsv
  write.table(sampleann, file = paste0(outputDir, "/", cell, "_sampleannot.tsv"), quote = FALSE, sep = "\t", col.names = T, row.names = FALSE)
  logger::log_success("Created sample annotation for ", cell)
}

#####################################################################
comparasions <- c(
  "C19_mild_vs_Ctrl", "C19_sev_vs_Ctrl", "HIV_acu_vs_Ctrl", "Influenza_ctrl_vs_d30",
  "HIV_chr_vs_Ctrl", "OP_high_vs_low", "OP_high_vs_med", "OP_low_vs_med"
)

for (comp in comparasions) {
  # write sampleannot as tsv
  sampleann %>%
    dplyr::select(Common_Minimal_Informative_ID, condition, bedFile, cell_type, dplyr::all_of(comp)) %>%
    na.omit(comp) %>%
    arrange(desc(row_number())) %>%
    write.table(file = paste0(outputDir, "/", cell, "_", comp, "_sampleannot.tsv"), quote = FALSE, sep = "\t", col.names = T, row.names = FALSE)
}
#####################################################################
