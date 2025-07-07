#!/usr/bin/env Rscript

#####################################################################
# 11_01_tcell_pseudobulks.R
# created on 26-05-25 by Irem Gunduz
# Create pseudobulks for T cells
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
})
set.seed(12)

cell_annot_meth <- "/icbb/projects/igunduz/DARPA_analysis/artemis_031023/sample_annot.tsv"
cell_annot_meth <- fread(cell_annot_meth) %>%
    filter(cell_type %in% c("Tc-Naive", "Tc-Mem", "Th-Naive", "Th-Mem")) 

# Remove _Rep1 and _Rep2 from Common_Minimal_Informative_ID
cell_annot_meth$Common_Minimal_Informative_ID <- gsub("_Rep[12]$", "", cell_annot_meth$Common_Minimal_Informative_ID)
cell_annot_meth$PB_id <- paste0(cell_annot_meth$cell_type, "_", cell_annot_meth$Common_Minimal_Informative_ID,".bedGraph")
cell_annot_meth$PB_id <- paste0("/icbb/projects/igunduz/DARPA/data/pseudoBulks/perSample/", cell_annot_meth$PB_id)

# Check if all pseudobulks are existing
pb_dir <- "/icbb/projects/igunduz/DARPA/data/pseudoBulks/perSample"
pb_files <- list.files(pb_dir, pattern = "\\.bedGraph$", full.names = TRUE)
missing_pbs <- setdiff(cell_annot_meth$PB_id,pb_files)
paste0("Missing pseudobulks: ", paste(missing_pbs, collapse = ", "))


# Subset Tc cells and keep all relevant columns, ensuring unique rows for `bedFile`
tcell_annot <- cell_annot_meth %>%
  filter(cell_type %in% c("Tc-Naive", "Tc-Mem")) %>%
  mutate(bedFile = PB_id) %>%  # Create a new column `bedFile` from `PB_id`
  dplyr::select(bedFile, cell_type, Common_Minimal_Informative_ID) %>%  # Keep only the desired columns
  distinct(bedFile, .keep_all = TRUE)  # Ensure rows are unique for `bedFile`

# Subset Th cells and ensure unique rows for `bedFile`
th_cell_annot <- cell_annot_meth %>%
  filter(cell_type %in% c("Th-Naive", "Th-Mem")) %>%
  mutate(bedFile = PB_id) %>%  # Create a new column `bedFile` from `PB_id`
  dplyr::select(bedFile, cell_type, Common_Minimal_Informative_ID) %>%  # Keep only the desired columns
  distinct(bedFile, .keep_all = TRUE)  # Ensure rows are unique for `bedFile`

# Remove "/icbb/projects/igunduz/DARPA/data/pseudoBulks/perSample/" from bedFile
tcell_annot$bedFile <- gsub("/icbb/projects/igunduz/DARPA/data/pseudoBulks/perSample/", "", tcell_annot$bedFile)
th_cell_annot$bedFile <- gsub("/icbb/projects/igunduz/DARPA/data/pseudoBulks/perSample/", "", th_cell_annot$bedFile)

# Write the results as tsv files
write.table(tcell_annot, 
            file = "/icbb/projects/igunduz/DARPA/data/pseudoBulks/perSample/Tc-Naive_vs_Tc-Mem.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

write.table(th_cell_annot, 
            file = "/icbb/projects/igunduz/DARPA/data/pseudoBulks/perSample/Th-Naive_vs_Th-Mem.tsv", 
            sep = "\t", row.names = FALSE, quote = FALSE)

####################################################################################################