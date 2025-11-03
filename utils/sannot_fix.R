## Load Libraries
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(readr)
  library(openxlsx)
  library(stringr)
})
set.seed(12)

# Load ArchR project
outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/"
project <- ArchR::loadArchRProject(outputDir, force = TRUE)
archr_samples <- unique(project@cellColData$Sample)
archr_samples <- archr_samples %>%
  sub("_fragments\\.tsv\\.gz$", "", .) %>%
  sub("\\.tsv\\.gz$", "", .) %>%
  sub("_fragments$", "", .)

# Load and filter sample metadata
sample_annot <- read_tsv("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/sample_annots/sampleAnnot_5x_v202201.tsv")
sample_annot <- sample_annot %>%
  filter(!grepl("^BA", sampleId))
sample_annot$exposure_group[sample_annot$exposure_group == "Influenza_d30"] <- "Influenza_d28"
sample_annot$exposure_grouping[sample_annot$exposure_grouping == "day30"] <- "day28"

# Extract QC metrics per cell and summarize by sample
qc_metrics <- as.data.frame(project@cellColData)
qc_summary <- qc_metrics %>%
  group_by(Sample) %>%
  summarise(
    n_cells = n(),
    mean_TSS = mean(TSSEnrichment, na.rm = TRUE),
    median_TSS = median(TSSEnrichment, na.rm = TRUE),
    mean_fragments = mean(log10(nFrags), na.rm = TRUE),
    median_fragments = median(log10(nFrags), na.rm = TRUE)
  )


# Load and clean sample annotations
sample_annot <- read_tsv("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/sample_annots/sampleAnnot_5x_v202201.tsv")
sample_annot <- sample_annot %>%
  filter(!grepl("^BA", sampleId))
sample_annot$exposure_group[sample_annot$exposure_group == "Influenza_d30"] <- "Influenza_d28"
sample_annot$exposure_grouping[sample_annot$exposure_grouping == "day30"] <- "day28"

# Merge sample annotations with QC metrics
sample_annot_qc <- left_join(sample_annot, qc_summary, by = c("arrow_name" = "Sample"))

# Save as Table S1
fig_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/sample_annots/figures"
tableS1_path <- file.path(fig_dir, "TableS1.xlsx")
write.xlsx(list("Sample_Metadata" = sample_annot_qc), file = tableS1_path, rowNames = FALSE)

# Supplementary tables (S2–S5, S6 becomes S7)
supp_files <- list.files(fig_dir, full.names = TRUE, pattern = "Supplementary_table_.*\\.(csv|tsv)$")
sheet_map <- c(
  "Supplementary_table_2" = "Table S2",
  "Supplementary_table_3" = "Table S3",
  "Supplementary_table_4" = "Table S4",
  "Supplementary_table_5" = "Table S5",
  "Supplementary_table_6" = "Table S7" # old S6 is now S7
)

# Description map
desc_map <- c(
  "Table S1" = "Sample metadata of the scATAC dataset",
  "Table S2" = "Results of pairwise Wilcoxon test for chromVAR deviation scores across different cell types",
  "Table S3" = "List of differentially accessible genes in different clusters within CD8+ T cells",
  "Table S4" = "List of differentially accessible peak regions in different clusters within CD8+ T cells",
  "Table S5" = "Differential gene activity and gene expression table of protein-coding genes for COVID-19 severe vs control in CD14+ monocytes",
  "Table S6" = "List of differentially accessible peak regions for COVID-19 severe vs control in CD14+ monocytes",
  "Table S7" = "Results of pairwise Wilcoxon test between one vs other cell type manner for methylTFR and chromVAR z-scores",
  "Table S8" = "List of differentially methylated peak regions in CD14+ monocytes"
)

# Read S2–S5 and S7 (renamed S6)
table_list <- list()
for (file_path in supp_files) {
  file_name <- basename(file_path)
  base <- tools::file_path_sans_ext(file_name)
  sheet_name <- sheet_map[[base]]
  df <- if (grepl("\\.tsv$", file_path)) {
    read_tsv(file_path, show_col_types = FALSE)
  } else {
    read_csv(file_path, show_col_types = FALSE)
  }
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
new_s6_path <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/Mono_CD14/reports/differential_data/diffTab_3_archrPeaks.tsv"
df_s6 <- read_tsv(new_s6_path, show_col_types = FALSE)
colnames(df_s6) <- colnames(df_s6) %>%
  str_replace_all("[\"']", "") %>%
  str_replace_all("\\.", "_") %>%
  str_replace_all("_+", "_") %>%
  str_trim()
df_s6_sig <- df_s6 %>%
  filter(!is.na(padj) & padj < 0.05)
table_list[["Table S6"]] <- df_s6_sig

# Read Table S8 (DMRs)
dmrs_path <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/sample_annots/Mono_CD14_dmr.tsv"
df_s8 <- read_tsv(dmrs_path, show_col_types = FALSE)
colnames(df_s8) <- colnames(df_s8) %>%
  str_replace_all("[\"']", "") %>%
  str_replace_all("\\.", "_") %>%
  str_replace_all("_+", "_") %>%
  str_trim()

df_s8_sig <- df_s8 %>%
  filter(!is.na(isDiff) & isDiff == TRUE)  
table_list[["Table S8"]] <- df_s8_sig

# Final sheet order
ordered_names <- paste0("Table S", 1:8)
table_list <- table_list[ordered_names]

# Create index sheet
index_df <- data.frame(
  Table = ordered_names,
  Description = unname(desc_map[ordered_names]),
  stringsAsFactors = FALSE
)
table_list <- c(list("Index" = index_df), table_list)

df_s5 <- table_list[["Table S5"]]

# Keep rows where isDiff1 or isDiff2 is TRUE
df_s5_filtered <- df_s5 %>%
  filter(isDiff_1 == TRUE | isDiff_2 == TRUE)

# Overwrite in the list
table_list[["Table S5"]] <- df_s5_filtered

# Write to Excel
write.xlsx(
  table_list,
  file = file.path(fig_dir, "All_Supplementary_Tables.xlsx"),
  rowNames = FALSE
)

message("✅ Supplementary workbook with Tables S1–S8 created at: All_Supplementary_Tables.xlsx")
