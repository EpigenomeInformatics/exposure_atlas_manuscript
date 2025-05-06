#!/usr/bin/env Rscript

#####################################################################
# new_pseudobulks.R
# created on 08-04-25 by Irem Gunduz
# Create pseudobulks for the methylation data based on cell-type
# for methylTFR usage
#####################################################################

suppressPackageStartupMessages({
  library(muLogR)
  library(dplyr)
  library(data.table)
})
source("/icbb/projects/igunduz/exposure_atlas_manuscript/src/utils/createPseudoBulks.R")
outputDir <- "/icbb/projects/igunduz/new_covid19_pseudobulks_090425"
if (!dir.exists(outputDir)) {
  dir.create(outputDir)
}

logger.start("Organizing metadata...")
annot_meth <- fread("/icbb/projects/igunduz/DARPA/allc_sample_annot_final.csv")
annot_meth <- dplyr::filter(annot_meth, condition %in% c("COVID_severe", "CommercialControl_healthy"))
annot_meth <- dplyr::filter(annot_meth, cell_type == "Monocyte")
cellAnnot_meth <- annot_meth
cellAnnot_meth <- dplyr::select(cellAnnot_meth, c("allC_FilePathfull", "Common_Minimal_Informative_ID", "cell_type", "condition"))
cellAnnot_meth$allC_FilePath_fixed <- gsub("/icbb/projects/igunduz/DARPA/", "/icbb/projects/share/datasets/ECHO/", cellAnnot_meth$allC_FilePathfull)
logger.completed()

outliers <- c(
  "CoV_S_S8_D1",
  "CoV_S_S15_D7",
  "Ctrl_1_F_White_45yo",
  "Ctrl_2_F_AfAm_43yo"
)

# Remove outliers from the cellAnnot_meth data frame
cellAnnot_meth <- cellAnnot_meth %>%
  filter(!Common_Minimal_Informative_ID %in% outliers)

mono_1 <- c(
  "CoV_S_S15_D1",
  "CoV_S_S11_D3",
  "CoV_S_S7_D1",
  "CoV_S_S11_D1"
)

# Create two groups: severe_mono1 and severe_mono2
severe_mono1 <- cellAnnot_meth %>%
  filter(Common_Minimal_Informative_ID %in% mono_1) %>%
  mutate(condition = "severe_mono1")

severe_mono2 <- cellAnnot_meth %>%
  filter(!Common_Minimal_Informative_ID %in% mono_1) %>%
  mutate(condition = "severe_mono2")

control <- cellAnnot_meth %>%
  filter(condition == "CommercialControl_healthy") %>%
  mutate(condition = "control")

# Combine the two groups
combined_annot <- bind_rows(severe_mono1, severe_mono2)

# Combine the control group
combined_annot <- bind_rows(combined_annot, control)

# Create cell-type pseudobulks for each group
createPseudoBulks(
  sampleAnnot = data.table::as.data.table(combined_annot),
  filePathCol = "allC_FilePath_fixed",
  sampleIdCol = "Common_Minimal_Informative_ID", groupName2 = NULL,
  numThreads = 30, mcType = "CGN", groupName = "condition", singleBP = TRUE,
  fileType = "allc", singleCovOff = 99999,
  indexed = FALSE, excludeChr = c("chrX", "chrY", "chrM", "chrL"),
  outputDir = outputDir
)

# Create sample annotation
new_annot <- data.frame(
  bedFile = paste0(outputDir, "/", unique(combined_annot$cell_type), ".bedGraph"),
  Exposure = unique(combined_annot$condition)
)
write.table(new_annot, file = paste0(outputDir, "/new_annot.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#####################################################################
