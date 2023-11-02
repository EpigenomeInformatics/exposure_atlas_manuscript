suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# Extract unique cells
cells <- data.table::fread("/icbb/projects/igunduz/DARPA/allc_sample_annot_final.csv") %>%
  dplyr::select(!V1) %>%
  dplyr::filter(!cell_type %in% c("Other-cell","Tc-Eff","Th-Eff"))

# Get the unique cell types
cells <- names(table(cells$cell_type))
cells <- sort(cells)
epp_dir <- "/icbb/projects/igunduz/DARPA_analysis/bedToEPP/EPP_271023/"

# All comparisons
comparisons <- c(
  "C19_sev_vs_Ctrl", 
  "C19_mild_vs_Ctrl", 
  "HIV_chr_vs_Ctrl",
  "HIV_acu_vs_Ctrl",
  "Influenza_ctrl_vs_d30",
  "OP_low_vs_med", 
  "OP_high_vs_med",
  "OP_high_vs_low"
)

for (cell in cells) {
  # Initialize an empty data frame to store sample information
  all_samples <- data.frame(sample_name = character(0), cell_type = character(0))
  for (comp in comparisons) {
    sample_dir <- paste0(epp_dir, cell, "/", comp)
    sample_ann <- "sample_methylation_summary.tsv"
    annfile <- file.path(sample_dir, sample_ann)
    if (endsWith(annfile, ".csv")) {
      samples <- read.table(annfile, header = TRUE, sep = ",", stringsAsFactors = FALSE)
    } else {
      samples <- read.table(annfile, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    }
    samples$bedFile <- sub(x=samples$bedFile, ".tsv", ".bed")

    files_list <- file.path(sample_dir, samples[, "bedFile"])
    samples$bedFile <-  files_list
    # Add the cell type to the samples data frame
    samples$cell_type <- cell
    samples <- samples %>% dplyr::select(Common_Minimal_Informative_ID,cell_type, bedFile, as.name(comp))
    colnames(samples) <- c("sample_name","cell_type", "bedFile","condition")
    
    # Merge the samples into the all_samples data frame
    all_samples <- dplyr::bind_rows(all_samples, samples)
  }
  # Remove duplicate sample names while keeping the first occurrence
  unique_samples <- all_samples %>% distinct(sample_name, .keep_all = TRUE)
  #unique_samples <- unique_samples %>% dplyr::filter(!condition %in% c("Influenza_d30","Influenza_ctrl","HIV_acu","HIV_ctrl","OP_low","OP_med"))
  #unique_samples <- unique_samples %>% dplyr::filter(condition %in% c("HIV_acu","HIV_ctrl","HIV_chr"))
  # Save the unique samples to a sample annotation file
  write.table(unique_samples, file = paste0(epp_dir,cell,"/sample_annotation_unique.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
