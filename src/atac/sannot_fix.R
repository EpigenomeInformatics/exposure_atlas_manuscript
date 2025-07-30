## Load Libraries
suppressPackageStartupMessages({library(ArchR)
library(dplyr)
library(readr)
library(openxlsx)
})
set.seed(12) # set seed

# Load ArchR project
outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/"
project <- ArchR::loadArchRProject(outputDir, force = T)
archr_samples <- unique(project@cellColData$Sample)
archr_samples <- archr_samples %>%
  sub("_fragments\\.tsv\\.gz$", "", .) %>%
  sub("\\.tsv\\.gz$", "", .) %>%
  sub("_fragments$", "", .)

# Extract cellColData
sample_annot <- read_tsv("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/sample_annots/sampleAnnot_5x_v202201.tsv")

# Remove rows where sampleId starts with "BA"
sample_annot <- sample_annot %>%
  filter(!grepl("^BA", sampleId))

# Save as Excel with one sheet
write.xlsx(
  list("Sample_Metadata" = sample_annot),
  file = "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/sample_annots/figures/TableS1.xlsx",
  rowNames = FALSE
)

# Load libraries
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(openxlsx)
  library(stringr)
})

# File paths
fig_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/sample_annots/figures"
new_s6_path <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/Mono_CD14/reports/differential_data/diffTab_3_archrPeaks.tsv"
tableS1_path <- file.path(fig_dir, "TableS1.xlsx")

# List and load S2–S5, S6 (which becomes S7)
supp_files <- list.files(fig_dir, full.names = TRUE, pattern = "Supplementary_table_.*\\.(csv|tsv)$")

# Sheet names: S2–S5 and S6 becomes S7
sheet_map <- c(
  "Supplementary_table_2" = "Table S2",
  "Supplementary_table_3" = "Table S3",
  "Supplementary_table_4" = "Table S4",
  "Supplementary_table_5" = "Table S5",
  "Supplementary_table_6" = "Table S7" # old S6 becomes S7
)

# Descriptions for index sheet
desc_map <- c(
  "Table S1" = "Sample metadata of the scATAC dataset",
  "Table S2" = "Results of pairwise Wilcoxon test for chromVAR deviation scores across different cell types",
  "Table S3" = "List of differentially accessible genes in different clusters within CD8+ T cells",
  "Table S4" = "List of differentially accessible peak regions in different clusters within CD8+ T cells",
  "Table S5" = "Differential gene activity and gene expression table of protein-coding genes for COVID-19 severe vs control in CD14+ monocytes",
  "Table S6" = "List of differentially accessible peak regions for COVID-19 severe vs control in CD14+ monocytes",
  "Table S7" = "Results of pairwise Wilcoxon test between one vs other cell type manner for methylTFR and chromVAR z-scores"
)

# Read and clean S2–S5, S7
table_list <- list()
for (file_path in supp_files) {
  file_name <- basename(file_path)
  base <- tools::file_path_sans_ext(file_name)
  sheet_name <- sheet_map[[base]]

  df <- if (grepl("\\.tsv$", file_path)) read_tsv(file_path, show_col_types = FALSE) else read_csv(file_path, show_col_types = FALSE)

  colnames(df) <- colnames(df) %>%
    str_replace_all("[\"']", "") %>%
    str_replace_all("\\.", "_") %>%
    str_replace_all("_+", "_") %>%
    str_trim()

  table_list[[sheet_name]] <- df
}

# Read Table S1
table_list[["Table S1"]] <- read.xlsx(tableS1_path)

# Read NEW Table S6 (ArchR peaks)
df_s6 <- read_tsv(new_s6_path, show_col_types = FALSE)
colnames(df_s6) <- colnames(df_s6) %>%
  str_replace_all("[\"']", "") %>%
  str_replace_all("\\.", "_") %>%
  str_replace_all("_+", "_") %>%
  str_trim()
table_list[["Table S6"]] <- df_s6

# Reorder all sheets according to the desired order
ordered_names <- paste0("Table S", 1:7)
table_list <- table_list[ordered_names]

# Create index sheet
index_df <- data.frame(
  Table = ordered_names,
  Description = unname(desc_map[ordered_names]),
  stringsAsFactors = FALSE
)
table_list <- c(list("Index" = index_df), table_list)

# Write to Excel
write.xlsx(
  table_list,
  file = file.path(fig_dir, "All_Supplementary_Tables.xlsx"),
  rowNames = FALSE
)

message("✅ Supplementary workbook created at: All_Supplementary_Tables.xlsx")
