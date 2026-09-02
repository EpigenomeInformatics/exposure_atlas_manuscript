#!/usr/bin/env Rscript

#####################################################################
# 12_assemble_supplementary_tables.R
# Single source of truth for the supplementary-table workbook + numbering.
#
# Run LAST, after the scripts that produce the individual tables:
#   - 01_v2_quality_control.R  -> All_Supplementary_Tables_updated.xlsx
#         (Table S1 with the recovered metadata / sex-prediction columns, and
#          the per-PC covariate association sheet)
#   - 04_2_cellprops.R         -> figures/cell_prop_wilcox_bonferroni.csv
#         (cell-type composition statistics)
#   - sample_annots/allc_sample_annot_final.csv
#         (snmC-seq per-cell annotation, pre-filter; subset below via
#          cellAnnot_meth.rds)
#   - src/atac/07_3_confounder_adjusted_DARs.R -> figures/confounder_adjusted_DAR_summary.csv
#     src/atac/07_6_covariate_balance.R          -> figures/covariate_balance_summary.csv
#         (covariate balance and the adjusted-vs-unadjusted DAR comparison)
#   - src/meth/01b_meth_pseudobulk_qc.R -> figures/meth_pseudobulk_qc.csv
#         and figures/meth_pseudobulk_qc_by_group.csv
#
# Numbering is declared once in `final_layout`; every rename derives from it.
# Reads the pre-shift workbook and writes a fresh file, so re-running is safe.
#
# FINAL NUMBERING (ATAC, then snmC-seq, then the combined-modality tables)
# Sequential throughout: no lettered sub-tables, every sheet carries its own
# number and an explicit label on the Index sheet.
#
#   S1    Sample metadata of the scATAC-seq dataset                 (was S1)
#   S2    Per-PC covariate association, before/after Harmony       (was S1B)
#   S3    Cluster-to-cell-type annotation mapping                   (was S2)
#   S4    Cell-type composition statistics                              [NEW]
#   S5    Variance partition of TF motif activity, exposure model  (was S3A)
#   S6    Variance partition of TF motif activity, condition model (was S3B)
#   S7    Pairwise differential TF motif activity across cell types  (was S4)
#   S8    Differentially accessible genes in CD8+ T-cell clusters    (was S5)
#   S9    Differentially accessible cCREs in CD8+ T-cell clusters    (was S6)
#   S10   Differential gene activity & expression, COVID-19 sev vs ctrl (was S7)
#   S11   Differentially accessible cCREs, COVID-19 severe vs control (was S8)
#   S12   Covariate balance between compared groups                     [NEW]
#   S13   Confounder-adjusted DARs vs unadjusted, sensitivity           [NEW]
#   S14   Bulk RNA-seq differential expression, CD14+ monocytes      (was S9)
#   S15   Variance partition, COVID-19 severity model               (was S3C)
#   S16   snmC-seq per-cell annotation                                  [NEW]
#   S17   snmC-seq per-sample summary (donor x timepoint x cell type)   [NEW]
#   S18   snmC-seq pseudobulk QC, per pseudobulk                        [NEW]
#   S19   snmC-seq pseudobulk QC, per exposure group x cell type        [NEW]
#   S20   Pairwise Wilcoxon (methylTFR vs chromVAR) across cell types (was S10)
#   S21   Differentially methylated cCREs, COVID-19 CD14+ monocytes  (was S11)
#
# The COVID-19 severity variance partition (old S3C) is placed after the
# COVID-19 differential tables rather than beside the other two models.

#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)  # pivot_wider() for the snmC per-sample cell counts
  library(tibble) # tribble() for the layout table
  library(openxlsx)
})

repo_dir <- "/scratch/icbb/igunduz/irem_github/exposure_atlas_manuscript"
fig_dir  <- file.path(repo_dir, "figures")
annot    <- file.path(repo_dir, "sample_annots")

supp_updated <- file.path(annot, "All_Supplementary_Tables_updated.xlsx")
supp_master  <- file.path(annot, "All_Supplementary_Tables.xlsx")
supp_in      <- if (file.exists(supp_updated)) supp_updated else supp_master
supp_out     <- file.path(annot, "All_Supplementary_Tables_final.xlsx")
comp_csv     <- file.path(fig_dir, "cell_prop_wilcox_bonferroni.csv")
allc_csv     <- file.path(annot, "allc_sample_annot_final.csv")
pbqc_csv     <- file.path(fig_dir, "meth_pseudobulk_qc.csv")
pbqc_grp_csv <- file.path(fig_dir, "meth_pseudobulk_qc_by_group.csv")

stopifnot(file.exists(supp_in), file.exists(comp_csv), file.exists(allc_csv))
if (!all(file.exists(pbqc_csv, pbqc_grp_csv))) {
  stop("Missing the methylation pseudobulk QC tables. Run ",
    "src/meth/01b_meth_pseudobulk_qc.R first.")
}
message("Reading workbook: ", supp_in)

## ---- 1. Old -> new numbering -----------------------------------------------
# `old` is the incoming sheet name, NA for new tables
final_layout <- tibble::tribble(
  ~new,          ~old,          ~description,
  "Table S1",    "Table S1",    "Sample metadata of the scATAC-seq dataset",
  "Table S2",    "Table S1B",   paste(
    "Association of the IterativeLSI and Harmony dimensions with known covariates",
    "at the sample level, before and after batch correction"),
  "Table S3",    NA,            paste(
    "Association of the IterativeLSI and Harmony dimensions with known covariates",
    "at the sample x cell-type pseudobulk level, the only level at which cell type",
    "varies, before and after batch correction"),
  "Table S4",    "Table S1C",   paste(
    "Correlation of each IterativeLSI and Harmony dimension with per-cell TSS",
    "enrichment, fragment number and FRIP, across all cells and after subtracting",
    "each sample's mean from its cells"),
  "Table S5",    "Table S1D",   paste(
    "Confound structure among the design and technical variables: cross-tabulation",
    "of cohort against supplier and per-sample filtering thresholds, and of the",
    "COVID-19 severity arms against processing batch"),
  "Table S6",    "Table S1E",   paste(
    "Performance of the molecular sex classifier: training-set size, leave-one-out",
    "accuracy with its confidence interval, agreement with the held-out chrY",
    "threshold, and the marker thresholds used"),
  "Table S7",    "Table S2",    paste(
    "Cluster-to-cell-type annotation mapping: Jaccard similarity between each",
    "graph-based scATAC cluster and predicted scRNA-seq labels, with final",
    "assigned cell type and cluster size"),
  "Table S8",    NA,            paste(
    "Cell-type composition statistics: per-sample cell-type proportions compared",
    "between each exposure group and its matched within-cohort control",
    "(two-sided Wilcoxon rank-sum test, Bonferroni-adjusted)"),
  "Table S9",    "Table S3A",   paste(
    "Variance partition of TF motif activity into cell-type, exposure and donor",
    "components (exposure model: ~ Cell_Type + Exposure + Donor, all cohorts)"),
  "Table S10",   "Table S3B",   paste(
    "Variance partition of TF motif activity into cell-type, stage/severity and",
    "donor components (condition model: ~ Cell_Type + Condition + Donor,",
    "all cohorts)"),
  "Table S11",   NA,            paste(
    "Variance in TF motif activity attributable to the exposure condition within",
    "each cell type, computed separately for each cohort",
    "(~ Condition + Donor, donor fitted where the design supports it)"),
  "Table S12",   "Table S4",    paste(
    "Pairwise differential TF motif activity (Wilcoxon test of chromVAR deviation",
    "scores) across different cell types"),
  "Table S13",   "Table S5",    "Differentially accessible genes in different clusters within CD8+ T cells",
  "Table S14",   "Table S6",    "Differentially accessible cCREs in different clusters within CD8+ T cells",
  "Table S15",   "Table S7",    paste(
    "Differential gene activity and gene expression of protein-coding genes",
    "for COVID-19 severe vs control in CD14+ monocytes"),
  "Table S16",   "Table S8",    "Differentially accessible cCREs for COVID-19 severe vs control in CD14+ monocytes",
  "Table S17",   NA,            paste(
    "Balance of the per-sample quality-control metrics between the groups compared",
    "in each within-cohort differential analysis: median and maximum standardised",
    "mean difference, the median expected under the null at these group sizes, and",
    "the number of comparisons in which the metric is imbalanced by permutation"),
  "Table S18",   NA,            paste(
    "Differentially accessible regions re-called with the per-sample quality-control",
    "metrics included in the differential model, compared with the unadjusted calls,",
    "per cell type, comparison and adjustment set: region counts, overlap, direction",
    "concordance and the correlation of fold changes"),
  "Table S19",   "Table S9",    paste(
    "Bulk RNA-seq differential expression results for CD14+ monocytes",
    "(filtered: adj p < 0.05 and |log2FC| > 0.5)"),
  "Table S20",   "Table S3C",   paste(
    "Variance partition of TF motif activity into cell-type, severity and donor",
    "components within the COVID-19 cohort only",
    "(severity model: ~ Cell_Type + Condition + Donor)"),
  "Table S21",   NA,            paste(
    "Per-cell annotation of the snmC-seq dataset: cell identifier, assigned cell",
    "type, donor sex and age, condition, mapping and coverage statistics, and",
    "global mCG / mCH / CCC methylation rates"),
  "Table S22",   NA,            paste(
    "Per-sample summary of the snmC-seq dataset: cohort, subject, timepoint,",
    "exposure level, donor sex and age, sequencing libraries, and the number of",
    "profiled nuclei in each FACS-sorted cell population, with per-sample means of",
    "the mapping, coverage and methylation-rate metrics"),
  "Table S23",   NA,            paste(
    "Quality-control statistics for each snmC-seq pseudobulk (one per cell type",
    "x sample): donor, exposure group, cells pooled, per-cell sequencing depth",
    "and mapping metrics, CpGs called, CpG coverage and its distribution, global",
    "mCG level, and the CCC rate as the bisulfite non-conversion estimate"),
  "Table S24",   NA,            paste(
    "snmC-seq pseudobulk quality control summarised per exposure group and cell",
    "type: number of pseudobulks and donors, cells pooled, CpGs called and global",
    "mCG level"),
  "Table S25",   "Table S10",   paste(
    "Pairwise Wilcoxon test between one vs other cell type manner for methylTFR",
    "and chromVAR z-scores"),
  "Table S26",   "Table S11",   "Differentially methylated cCREs in CD14+ monocytes"
)

## ---- 2. New table content ---------------------------------------------------
# (a) cell-type composition statistics
comp_supp <- read.csv(comp_csv, stringsAsFactors = FALSE) %>%
  dplyr::transmute(
    `Cell type`                   = cell_type,
    `Group 1`                     = group1,
    `Group 2`                     = group2,
    `n (group 1)`                 = n1,
    `n (group 2)`                 = n2,
    `Median proportion (group 1)` = signif(median1, 3),
    `Median proportion (group 2)` = signif(median2, 3),
    `p (Wilcoxon rank-sum)`       = signif(p_value, 3),
    `p (Bonferroni-adjusted)`     = signif(p_adj, 3)
  )

# (b) snmC-seq per-cell annotation (Table S16). Drops the row-number column and the absolute
# allC_FilePathfull; the relative allC_FilePath is kept.
allc_raw <- read.csv(allc_csv, stringsAsFactors = FALSE, check.names = FALSE)
drop_cols <- c("", "X", "allC_FilePathfull")
allc_supp <- allc_raw[, !colnames(allc_raw) %in% drop_cols, drop = FALSE]
message("snmC-seq annotation as read: ", nrow(allc_supp), " cells x ",
        ncol(allc_supp), " columns (dropped: ",
        paste(intersect(drop_cols, colnames(allc_raw)), collapse = ", "), ")")

# Restrict to the analysed cell set. The CSV is pre-filter (115 samples,
# 81,907 cells); the Results report 22,488. Filter chain, in src/integration/:
#   01_prepare_sampleannot.R  COVID/FLU/HIV/OP                -> 105 samples
#                             intersect with ATAC sample IDs  ->  39 samples
#                             >= 200 ATAC and >= 50 snmC      ->  38, 27,225 cells
#   Other-cell dropped (7 FACS populations)                   ->  38, 23,981 cells
#   03_quality_check.R        N_valid_sites in [5e5, 4e6]     ->      22,488 cells
# The last step needs N_valid_sites, which the CSV does not carry, so use the
# saved post-QC annotation where available.
meth_qc_rds <- "/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/cellAnnot_meth.rds"
n_cells_expected <- 22488

if (file.exists(meth_qc_rds)) {
  qc_cells <- readRDS(meth_qc_rds)
  cell_ids <- if ("Cell_UID" %in% colnames(qc_cells)) qc_cells[["Cell_UID"]] else rownames(qc_cells)
  cell_ids <- unique(stats::na.omit(as.character(cell_ids)))
  n_before <- nrow(allc_supp)
  allc_supp <- allc_supp[allc_supp$Cell_UID %in% cell_ids, , drop = FALSE]
  message("Restricted to the post-QC cell set from cellAnnot_meth.rds: ",
          nrow(allc_supp), " of ", n_before, " cells")
  if (nrow(allc_supp) != n_cells_expected) {
    warning("Expected ", n_cells_expected, " cells, got ", nrow(allc_supp),
            "; check cellAnnot_meth.rds against the Results.")
  }
} else {
  warning("cellAnnot_meth.rds not found at ", meth_qc_rds,
          "; falling back to the repo filters, which omit the N_valid_sites step. ",
          "Cell count will exceed the ", n_cells_expected, " in the Results.")
  atac_cc <- file.path(annot, "cellColData.tsv.gz")
  stopifnot(file.exists(atac_cc))
  cc <- read.delim(gzfile(atac_cc), stringsAsFactors = FALSE)
  keep <- allc_supp$Common_Minimal_Informative_ID %in%
    intersect(unique(cc$sample_sampleId_cminid),
              unique(allc_supp$Common_Minimal_Informative_ID))
  allc_supp <- allc_supp[keep, , drop = FALSE]
  n_atac <- table(cc$sample_sampleId_cminid)
  n_meth <- table(allc_supp$Common_Minimal_Informative_ID)
  ok <- names(n_meth)[n_meth >= 50 &
                        n_atac[names(n_meth)] >= 200 &
                        !is.na(n_atac[names(n_meth)])]
  allc_supp <- allc_supp[allc_supp$Common_Minimal_Informative_ID %in% ok, , drop = FALSE]
  allc_supp <- allc_supp[allc_supp$cell_type != "Other-cell", , drop = FALSE]
  message("Fallback subset: ", nrow(allc_supp), " cells, ",
          length(unique(allc_supp$Common_Minimal_Informative_ID)), " samples")
}

message("snmC-seq annotation for S16: ", nrow(allc_supp), " cells, ",
        length(unique(allc_supp$Common_Minimal_Informative_ID)), " samples")
if (nrow(allc_supp) > 1e5) {
  message("NB this sheet is large; Excel's row limit is 1,048,576 so it fits, ",
          "but the workbook will be slow to open.")
}

# (b2) snmC-seq per-sample summary: cells per donor, timepoint and cell type.
# Aggregated from allc_supp so it cannot disagree with the per-cell sheet.
# Subject / timepoint / exposure level are parsed from the sample identifier:
#   CoV_S_S11_D1          severity (S / nS), subject, day
#   Ctrl_10_M_White_39yo  index, sex, ethnicity, age
#   Flu_S10_D28           subject, day
#   HIV_S1_Pre            subject, stage (Pre / Acu / Cro)
#   OP_S12_High_D28       subject, exposure level, day (or NonFarm)
ct_levels <- c("B-cell", "Monocyte", "NK-cell",
               "Th-Naive", "Th-Mem", "Th-Eff",
               "Tc-Naive", "Tc-Mem", "Tc-Eff", "Other-cell")
# Other-cell is absent from the analysed set, so intersect rather than require
# equality
ct_levels <- intersect(ct_levels, unique(allc_supp$cell_type))
stopifnot(length(ct_levels) > 0,
          all(unique(allc_supp$cell_type) %in% ct_levels))

sid   <- allc_supp$Common_Minimal_Informative_ID
part  <- function(x, i) vapply(strsplit(x, "_", fixed = TRUE), function(p)
  if (length(p) >= i) p[i] else NA_character_, character(1))
cohort <- dplyr::case_when(
  startsWith(sid, "CoV_")  ~ "COVID-19",
  startsWith(sid, "Ctrl_") ~ "Control",
  startsWith(sid, "Flu_")  ~ "Influenza",
  startsWith(sid, "HIV_")  ~ "HIV",
  startsWith(sid, "OP_")   ~ "OP",
  TRUE ~ NA_character_
)
stopifnot(!any(is.na(cohort)))

allc_keyed <- allc_supp %>%
  dplyr::mutate(
    .sid    = sid,
    Cohort  = cohort,
    Subject = dplyr::case_when(
      Cohort == "COVID-19" ~ part(.sid, 3),
      Cohort == "Control"  ~ paste0("Ctrl", part(.sid, 2)),
      TRUE                 ~ part(.sid, 2)
    ),
    Timepoint = dplyr::case_when(
      Cohort == "COVID-19"  ~ part(.sid, 4),
      Cohort == "Influenza" ~ part(.sid, 3),
      Cohort == "HIV"       ~ dplyr::recode(part(.sid, 3),
                                Pre = "pre", Acu = "acute", Cro = "chronic"),
      Cohort == "OP"        ~ part(.sid, 4),
      TRUE                  ~ NA_character_
    ),
    `Exposure level` = dplyr::case_when(
      Cohort == "COVID-19" ~ dplyr::recode(part(.sid, 2), S = "severe", nS = "mild"),
      Cohort == "OP"       ~ part(.sid, 3),
      Cohort == "Control"  ~ "healthy",
      TRUE                 ~ NA_character_
    ),
    cell_type = factor(cell_type, levels = ct_levels)
  )

allc_meta <- allc_keyed %>%
  dplyr::group_by(`Sample ID` = .sid) %>%
  dplyr::summarise(
    Cohort           = dplyr::first(Cohort),
    Condition        = dplyr::first(condition),
    Subject          = dplyr::first(Subject),
    Timepoint        = dplyr::first(Timepoint),
    `Exposure level` = dplyr::first(`Exposure level`),
    Sex              = dplyr::first(stats::na.omit(Sex))[1],
    Age              = dplyr::first(stats::na.omit(age))[1],
    Libraries        = dplyr::n_distinct(Salk_ID),
    `Library IDs`    = paste(sort(unique(Salk_ID)), collapse = "; "),
    `Total cells`    = dplyr::n(),
    `Mean unique mapped reads` = round(mean(TotalUniqueMappedReads), 0),
    `Mean mapped ratio (%)`    = round(mean(TotalMappedRatio), 2),
    `Mean mCG rate`            = round(mean(CG_Rate), 4),
    `Mean mCH rate`            = round(mean(CH_Rate), 4),
    `Mean CCC rate`            = round(mean(CCC_Rate), 4),
    `Mean genome coverage`     = round(mean(genome_cov), 4),
    .groups = "drop"
  )

allc_counts <- allc_keyed %>%
  dplyr::count(`Sample ID` = .sid, cell_type, .drop = FALSE) %>%
  tidyr::pivot_wider(names_from = cell_type, values_from = n, values_fill = 0)

allc_sample_supp <- allc_meta %>%
  dplyr::left_join(allc_counts, by = "Sample ID") %>%
  dplyr::relocate(dplyr::all_of(ct_levels), .after = `Total cells`) %>%
  dplyr::arrange(Cohort, Subject, Timepoint)

# the join must not add rows, and the per-cell-type counts must sum to the total
stopifnot(nrow(allc_sample_supp) == dplyr::n_distinct(sid))
stopifnot(sum(allc_sample_supp$`Total cells`) == nrow(allc_supp))
stopifnot(all(rowSums(allc_sample_supp[, ct_levels]) ==
                allc_sample_supp$`Total cells`))
message("snmC-seq per-sample summary: ", nrow(allc_sample_supp), " samples, ",
        dplyr::n_distinct(allc_sample_supp$Subject), " subjects")

# Sex should not vary within a subject. Four do in the committed annotation
# (COVID-19 S11, S17; OP S4, S7); age is consistent. Warning, not an error.
sex_clash <- allc_sample_supp %>%
  dplyr::group_by(Cohort, Subject) %>%
  dplyr::filter(dplyr::n_distinct(stats::na.omit(Sex)) > 1) %>%
  dplyr::ungroup()
if (nrow(sex_clash)) {
  warning(nrow(sex_clash), " samples in ",
          dplyr::n_distinct(sex_clash$Subject),
          " subject(s) have inconsistent Sex across timepoints; resolve before ",
          "submission (see Reviewer 1, comment 1)")
  print(as.data.frame(sex_clash[, c("Sample ID", "Cohort", "Subject",
                                    "Timepoint", "Sex", "Age")]))
}

# (c) snmC-seq pseudobulk QC (src/meth/01b_meth_pseudobulk_qc.R)
pbqc_supp     <- read.csv(pbqc_csv, stringsAsFactors = FALSE, check.names = FALSE)
pbqc_grp_supp <- read.csv(pbqc_grp_csv, stringsAsFactors = FALSE, check.names = FALSE)
message("snmC-seq pseudobulk QC: ", nrow(pbqc_supp), " pseudobulks, ",
        nrow(pbqc_grp_supp), " exposure group x cell type combinations")

# (d) Covariate balance (07_6_covariate_balance.R) and the adjusted-vs-unadjusted
# DAR comparison (07_3_confounder_adjusted_DARs.R), as separate sheets.
bal_csv <- file.path(fig_dir, "covariate_balance_summary.csv")
adj_csv <- file.path(fig_dir, "confounder_adjusted_DAR_summary.csv")
stopifnot(file.exists(bal_csv), file.exists(adj_csv))

bal_supp <- read.csv(bal_csv, stringsAsFactors = FALSE) %>%
  dplyr::transmute(
    `QC metric`                              = covariate,
    `Comparisons tested`                     = n_comparisons,
    `Median |SMD|`                           = signif(median_abs_SMD, 3),
    `Median |SMD| expected under the null`    = signif(null_median_abs_SMD, 3),
    `Max |SMD|`                              = signif(max_abs_SMD, 3),
    `Comparisons imbalanced (permutation)`   = n_imbalanced_perm,
    `Smallest permutation p`                 = signif(min_p_perm, 3)
  )

adj_supp <- read.csv(adj_csv, stringsAsFactors = FALSE) %>%
  dplyr::transmute(
    `Cell type`                    = cell,
    `Comparison`                   = comparison,
    `Adjustment set`               = adj_set,
    `DARs, unadjusted`             = DARs_unadjusted,
    `DARs, adjusted`               = DARs_adjusted,
    `Shared`                       = shared,
    `Recovered (%)`                = signif(recovered_pct, 3),
    `Direction concordance (%)`    = signif(sign_concordance_pct, 3),
    `Fold-change correlation (r)`  = signif(lfc_pearson_all, 3),
    `Regions tested`               = n_regions_tested
  ) %>%
  dplyr::arrange(`Cell type`, `Comparison`, `Adjustment set`)

message("Covariate balance: ", nrow(bal_supp), " metrics, ",
        sum(bal_supp$`Comparisons imbalanced (permutation)`), " imbalanced comparisons in total")
message("Confounder-adjusted DARs: ", nrow(adj_supp), " cell x comparison x design rows")

ct_assoc_csv <- file.path(fig_dir, "dim_covariate_association_by_celltype.csv")
cc_varpart_csv <- file.path(repo_dir, "sample_annots/cvar_diffs",
  "varpart_within_celltype_by_cohort.csv")
read_optional <- function(path, label) {
  if (!file.exists(path)) {
    warning(label, " not found at ", path, "; that sheet will be skipped.")
    return(NULL)
  }
  read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
}
ct_assoc_supp   <- read_optional(ct_assoc_csv, "Cell-type association table")
cc_varpart_supp <- read_optional(cc_varpart_csv, "Per-cohort variance partition")

new_content <- list(
  "Table S3"   = ct_assoc_supp,
  "Table S8"   = comp_supp,
  "Table S11"  = cc_varpart_supp,
  "Table S17"  = bal_supp,
  "Table S18"  = adj_supp,
  "Table S21"  = allc_supp,
  "Table S22"  = allc_sample_supp,
  "Table S23"  = pbqc_supp,
  "Table S24"  = pbqc_grp_supp
)
new_content <- new_content[!vapply(new_content, is.null, logical(1))]

## ---- 3. Rename existing sheets ---------------------------------------------
wb <- openxlsx::loadWorkbook(supp_in)

if (all(names(new_content) %in% names(wb)) &&
    identical(dim(openxlsx::readWorkbook(wb, "Table S8")), dim(comp_supp))) {
  message("Workbook already assembled; nothing to do.")
} else {
  # temp prefix first: the old and new names overlap (old S4 -> S7 while an
  # old S7 still exists, old S10/S11 -> S20/S21), so nothing may be renamed
  # straight onto a name still in use
  renames <- final_layout %>%
    dplyr::filter(!is.na(old), old != new)

  tmp_map <- list()
  for (i in seq_len(nrow(renames))) {
    old_i <- renames$old[i]
    new_i <- renames$new[i]
    if (old_i %in% names(wb)) {
      tmp_i <- paste0("__tmp__", new_i)
      openxlsx::renameWorksheet(wb, old_i, tmp_i)
      tmp_map[[tmp_i]] <- new_i
    }
  }
  for (tmp_i in names(tmp_map)) {
    openxlsx::renameWorksheet(wb, tmp_i, tmp_map[[tmp_i]])
  }
  message("Renamed ", length(tmp_map), " sheet(s)")

  ## ---- 4. Insert the new tables ---------------------------------------------
  hdr <- openxlsx::createStyle(textDecoration = "bold")
  for (nm in names(new_content)) {
    if (nm %in% names(wb)) openxlsx::removeWorksheet(wb, nm)
    openxlsx::addWorksheet(wb, nm)
    dat <- new_content[[nm]]
    openxlsx::writeData(wb, nm, dat, withFilter = TRUE, headerStyle = hdr)
    # auto-width is slow on tens of thousands of rows
    if (nrow(dat) < 5000) {
      openxlsx::setColWidths(wb, nm, cols = seq_along(dat), widths = "auto")
    } else {
      openxlsx::setColWidths(wb, nm, cols = seq_along(dat), widths = 18)
    }
    message("Wrote ", nm, " (", nrow(dat), " rows)")
  }

  ## ---- 5. Order the sheets to match final_layout ----------------------------
  # one sheet per final_layout row, so the order is the layout order
  want <- final_layout$new[final_layout$new %in% names(wb)]
  other <- setdiff(names(wb), c(want, "Index"))
  ord_names <- c("Index"[("Index" %in% names(wb))], want, other)
  openxlsx::worksheetOrder(wb) <- match(ord_names, names(wb))

  ## ---- 6. Rebuild the Index sheet -------------------------------------------
  # built from final_layout, so it cannot drift from the sheet names
  idx <- data.frame(
    Table = final_layout$new,
    Description = final_layout$description,
    stringsAsFactors = FALSE
  )
  idx <- idx[idx$Table %in% names(wb), , drop = FALSE]

  if ("Index" %in% names(wb)) openxlsx::removeWorksheet(wb, "Index")
  openxlsx::addWorksheet(wb, "Index")
  openxlsx::writeData(wb, "Index", idx, headerStyle = hdr)
  openxlsx::setColWidths(wb, "Index", cols = 1:2, widths = c(14, 110))
  openxlsx::worksheetOrder(wb) <- c(
    which(names(wb) == "Index"),
    setdiff(seq_along(names(wb)), which(names(wb) == "Index"))
  )

  openxlsx::saveWorkbook(wb, supp_out, overwrite = TRUE)
  message("Wrote ", supp_out)

  ## ---- 7. Renumbering map for the manuscript text ---------------------------
  map_df <- final_layout %>%
    dplyr::filter(!is.na(old)) %>%
    dplyr::transmute(old, new, changed = old != new)
  message("\nAPPLY THIS RENUMBERING TO THE MANUSCRIPT TEXT ",
          "(work from the HIGHEST number down, so replacements do not collide):")
  chg <- as.data.frame(map_df[map_df$changed, c("old", "new")])
  print(chg[order(-as.numeric(gsub("[^0-9]", "", chg$old)), chg$old), ])
  message("\nNew tables: S4 = cell-type composition statistics, ",
          "S12/S13 = covariate balance and confounder-adjusted DARs, ",
          "S16 = snmC-seq per-cell annotation, ",
          "S17 = snmC-seq per-sample summary, ",
          "S18/S19 = snmC-seq pseudobulk QC.")
  message("NB the numbering is flat: the old lettered sub-tables are now ",
          "tables in their own right (S1B -> S2, S3A/S3B -> S5/S6, ",
          "S3C -> S15, and the new S12/S13, S17, S19).")
}

#####################################################################
