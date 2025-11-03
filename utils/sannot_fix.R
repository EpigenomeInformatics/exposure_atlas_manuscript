# -------------------------------------------------------------------------
# Create Supplementary Workbook with Tables S1–S9
# Includes bulk RNA-seq differential expression (Table S7)
# -------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(readr)
  library(openxlsx)
  library(stringr)
})
set.seed(12)

# -------------------------------------------------------------------------
# Paths
# -------------------------------------------------------------------------
outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/"
fig_dir   <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/sample_annots/figures"
sannot_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/sample_annots/"

# -------------------------------------------------------------------------
# Load ArchR project
# -------------------------------------------------------------------------
project <- ArchR::loadArchRProject(outputDir, force = TRUE)

# -------------------------------------------------------------------------
# Sample annotation and QC
# -------------------------------------------------------------------------
sample_annot <- read_tsv(file.path(sannot_dir, "sampleAnnot_5x_v202201.tsv"), show_col_types = FALSE) %>%
  filter(!grepl("^BA", sampleId))
sample_annot$exposure_group[sample_annot$exposure_group == "Influenza_d30"] <- "Influenza_d28"
sample_annot$exposure_grouping[sample_annot$exposure_grouping == "day30"] <- "day28"

qc_metrics <- as.data.frame(project@cellColData)
qc_summary <- qc_metrics %>%
  group_by(Sample) %>%
  summarise(
    n_cells = n(),
    mean_TSS = mean(TSSEnrichment, na.rm = TRUE),
    median_TSS = median(TSSEnrichment, na.rm = TRUE),
    mean_fragments = mean(log10(nFrags), na.rm = TRUE),
    median_fragments = median(log10(nFrags), na.rm = TRUE),
    .groups = "drop"
  )

sample_annot_qc <- left_join(sample_annot, qc_summary, by = c("arrow_name" = "Sample"))
tableS1_path <- file.path(fig_dir, "TableS1.xlsx")
write.xlsx(list("Sample_Metadata" = sample_annot_qc), file = tableS1_path, rowNames = FALSE)

# -------------------------------------------------------------------------
# Load other supplementary tables (S2–S5)
# -------------------------------------------------------------------------
supp_files <- list.files(fig_dir, full.names = TRUE, pattern = "Supplementary_table_.*\\.(csv|tsv)$")
sheet_map <- c(
  "Supplementary_table_2" = "Table S2",
  "Supplementary_table_3" = "Table S3",
  "Supplementary_table_4" = "Table S4",
  "Supplementary_table_5" = "Table S5"
)

table_list <- list()
for (file_path in supp_files) {
  base <- tools::file_path_sans_ext(basename(file_path))
  sheet_name <- sheet_map[[base]]
  if (is.null(sheet_name)) next
  df <- if (grepl("\\.tsv$", file_path)) read_tsv(file_path, show_col_types = FALSE) else read_csv(file_path, show_col_types = FALSE)
  colnames(df) <- colnames(df) %>%
    str_replace_all("[\"']", "") %>%
    str_replace_all("\\.", "_") %>%
    str_replace_all("_+", "_") %>%
    str_trim()
  table_list[[sheet_name]] <- df
}

table_list[["Table S1"]] <- read.xlsx(tableS1_path)

# -------------------------------------------------------------------------
# Table S6: differential ArchR peaks
# -------------------------------------------------------------------------
new_s6_path <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/Mono_CD14/reports/differential_data/diffTab_3_archrPeaks.tsv"
df_s6 <- read_tsv(new_s6_path, show_col_types = FALSE)
colnames(df_s6) <- colnames(df_s6) %>%
  str_replace_all("[\"']", "") %>%
  str_replace_all("\\.", "_") %>%
  str_replace_all("_+", "_") %>%
  str_trim()
table_list[["Table S6"]] <- df_s6 %>% filter(!is.na(padj) & padj < 0.05)

# -------------------------------------------------------------------------
# Table S8: DMRs
# -------------------------------------------------------------------------
dmrs_path <- file.path(sannot_dir, "Mono_CD14_dmr.tsv")
df_s8 <- read_tsv(dmrs_path, show_col_types = FALSE)
colnames(df_s8) <- colnames(df_s8) %>%
  str_replace_all("[\"']", "") %>%
  str_replace_all("\\.", "_") %>%
  str_replace_all("_+", "_") %>%
  str_trim()
table_list[["Table S9"]] <- df_s8 %>% filter(!is.na(isDiff) & isDiff == TRUE)

# -------------------------------------------------------------------------
# Filter Table S5 (gene-level diff)
# -------------------------------------------------------------------------
if ("Table S5" %in% names(table_list)) {
  df_s5 <- table_list[["Table S5"]]
  if (all(c("isDiff_1", "isDiff_2") %in% names(df_s5))) {
    table_list[["Table S5"]] <- df_s5 %>% filter(isDiff_1 == TRUE | isDiff_2 == TRUE)
  }
}

# -------------------------------------------------------------------------
# Table S7: bulk monocyte DE (filtered)
# -------------------------------------------------------------------------
bulk_path <- file.path(sannot_dir, "gene_expression_protein_coding_diffs.rds")
bulk.mono.de <- readRDS(bulk_path)

if (all(c("p_val_adj", "avg_log2FC") %in% names(bulk.mono.de))) {
  sig_genes <- bulk.mono.de %>% filter(p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
} else if (all(c("FDR_rna", "Log2FC_rna") %in% names(bulk.mono.de))) {
  sig_genes <- bulk.mono.de %>% filter(FDR_rna < 0.05 & abs(Log2FC_rna) > 0.5)
} else if (all(c("padj", "log2FoldChange") %in% names(bulk.mono.de))) {
  sig_genes <- bulk.mono.de %>% filter(padj < 0.05 & abs(log2FoldChange) > 0.5)
} else {
  stop("bulk.mono.de lacks recognizable p-value and logFC columns.")
}
table_list[["Table S7"]] <- sig_genes

# -------------------------------------------------------------------------
# Table S8: Wilcoxon methylTFR/chromVAR (old S7)
# -------------------------------------------------------------------------
if (!is.null(table_list[["Table S7"]])) {
  # If there was an older Table S7 from supplementary files, move to S8
  table_list[["Table S8"]] <- table_list[["Table S7"]]
}

# -------------------------------------------------------------------------
# Order and descriptions
# -------------------------------------------------------------------------
ordered_names <- paste0("Table S", 1:9)

desc_map <- c(
  "Table S1" = "Sample metadata of the scATAC dataset",
  "Table S2" = "Results of pairwise Wilcoxon test for chromVAR deviation scores across different cell types",
  "Table S3" = "List of differentially accessible genes in different clusters within CD8+ T cells",
  "Table S4" = "List of differentially accessible peak regions in different clusters within CD8+ T cells",
  "Table S5" = "Differential gene activity and gene expression table of protein-coding genes for COVID-19 severe vs control in CD14+ monocytes",
  "Table S6" = "List of differentially accessible peak regions for COVID-19 severe vs control in CD14+ monocytes",
  "Table S7" = "Bulk RNA-seq differential expression results for CD14+ monocytes (filtered: adj p < 0.05 & |log2FC| > 0.5)",
  "Table S8" = "Results of pairwise Wilcoxon test between one vs other cell type manner for methylTFR and chromVAR z-scores",
  "Table S9" = "List of differentially methylated peak regions in CD14+ monocytes"
)

index_df <- data.frame(
  Table = ordered_names,
  Description = unname(desc_map[ordered_names]),
  stringsAsFactors = FALSE
)

# -------------------------------------------------------------------------
# Final write
# -------------------------------------------------------------------------
final_table_list <- c(list("Index" = index_df), table_list[ordered_names])
out_path <- file.path(fig_dir, "All_Supplementary_Tables.xlsx")
write.xlsx(final_table_list, file = out_path, rowNames = FALSE)

message("✅ Supplementary workbook with Tables S1–S9 created at: ", out_path)