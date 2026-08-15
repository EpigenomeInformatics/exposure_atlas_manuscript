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
#         (snmC-seq per-cell annotation; committed in the repo)
#   - src/meth/01b_meth_pseudobulk_qc.R -> figures/meth_pseudobulk_qc.csv
#         and figures/meth_pseudobulk_qc_by_group.csv
#
# Several new tables are inserted here, so the numbering shifts more than once. Rather than
# chaining two rename passes (which is where off-by-one errors creep in), the
# final layout is declared once in `final_layout` below and every rename is
# derived from it. Edit that table if the order changes; nothing else needs
# touching.
#
# The script always reads the pre-shift workbook and writes a fresh final file,
# so it is idempotent (safe to re-run).
#
# FINAL NUMBERING
#   S1   Sample metadata of the scATAC dataset
#   S1B  Per-PC covariate association (before/after Harmony batch correction)
#   S2   snmC-seq per-cell annotation                                     [NEW]
#   S3   snmC-seq pseudobulk QC, per pseudobulk                           [NEW]
#   S3B  snmC-seq pseudobulk QC, per exposure group x cell type           [NEW]
#   S4   Cluster-to-cell-type annotation mapping                     (was S2)
#   S5   Cell-type composition statistics                                 [NEW]
#   S6   Variance partition of TF motif activity (S6A/B/C)      (was S3A/B/C)
#   S7   Pairwise differential TF motif activity across cell types   (was S4)
#   S8   Differentially accessible genes in CD8+ T-cell clusters    (was S5)
#   S9   Differentially accessible cCREs in CD8+ T-cell clusters    (was S6)
#   S10  Differential gene activity & expression, COVID-19 sev vs ctrl (was S7)
#   S11  Differentially accessible cCREs, COVID-19 severe vs control (was S8)
#   S12  Bulk RNA-seq differential expression, CD14+ monocytes      (was S9)
#   S13  Pairwise Wilcoxon (methylTFR vs chromVAR) across cell types (was S10)
#   S14  Differentially methylated cCREs, COVID-19 CD14+ monocytes  (was S11)
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble) # tribble() for the layout table
  library(openxlsx)
})

repo_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript"
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
# `old` is the sheet name in the incoming workbook, NA for the two new tables.
final_layout <- tibble::tribble(
  ~new,          ~old,          ~description,
  "Table S1",    "Table S1",    "Sample metadata of the scATAC dataset",
  "Table S1B",   "Table S1B",   paste(
    "Association of sample-level principal components with known covariates,",
    "before (IterativeLSI) and after (Harmony) batch correction"),
  "Table S2",    NA,            paste(
    "Per-cell annotation of the snmC-seq dataset: cell identifier, assigned cell",
    "type, donor sex and age, condition, mapping and coverage statistics, and",
    "global mCG / mCH / CCC methylation rates"),
  "Table S3",    NA,            paste(
    "Quality-control statistics for each snmC-seq pseudobulk (one per cell type",
    "x sample): donor, exposure group, cells pooled, per-cell sequencing depth",
    "and mapping metrics, CpGs called, CpG coverage and its distribution, global",
    "mCG level, and the CCC rate as the bisulfite non-conversion estimate"),
  "Table S3B",   NA,            paste(
    "snmC-seq pseudobulk quality control summarised per exposure group and cell",
    "type: number of pseudobulks and donors, cells pooled, CpGs called and global",
    "mCG level"),
  "Table S4",    "Table S2",    paste(
    "Cluster-to-cell-type annotation mapping: Jaccard similarity between each",
    "graph-based scATAC cluster and predicted scRNA-seq labels, with final",
    "assigned cell type and cluster size"),
  "Table S5",    NA,            paste(
    "Cell-type composition statistics: per-sample cell-type proportions compared",
    "between each exposure group and its matched within-cohort control",
    "(two-sided Wilcoxon rank-sum test, Bonferroni-adjusted)"),
  "Table S6",    "Table S3",    paste(
    "Variance partition of TF motif activity across cell type, exposure and donor",
    "components (S6A: exposure model, S6B: stage/severity model,",
    "S6C: COVID-19 severity)"),
  "Table S7",    "Table S4",    paste(
    "Pairwise differential TF motif activity (Wilcoxon test of chromVAR deviation",
    "scores) across different cell types"),
  "Table S8",    "Table S5",    "List of differentially accessible genes in different clusters within CD8+ T cells",
  "Table S9",    "Table S6",    "List of differentially accessible cCREs in different clusters within CD8+ T cells",
  "Table S10",   "Table S7",    paste(
    "Differential gene activity and gene expression table of protein-coding genes",
    "for COVID-19 severe vs control in CD14+ monocytes"),
  "Table S11",   "Table S8",    "Differentially accessible cCREs for COVID-19 severe vs control in CD14+ monocytes",
  "Table S12",   "Table S9",    paste(
    "Bulk RNA-seq differential expression results for CD14+ monocytes",
    "(filtered: adj p < 0.05 & |log2FC| > 0.5)"),
  "Table S13",   "Table S10",   paste(
    "Pairwise Wilcoxon test between one vs other cell type manner for methylTFR",
    "and chromVAR z-scores"),
  "Table S14",   "Table S11",   "Differentially methylated cCREs in CD14+ monocytes"
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

# (b) snmC-seq per-cell annotation.
# Two columns are deliberately dropped:
#   - the unnamed row-number column written by write.csv (no information)
#   - allC_FilePathfull, an absolute path under /icbb/projects/igunduz/... The
#     relative allC_FilePath is kept; publishing internal filesystem layout in a
#     supplementary table serves no reader and is easy to overlook.
allc_raw <- read.csv(allc_csv, stringsAsFactors = FALSE, check.names = FALSE)
drop_cols <- c("", "X", "allC_FilePathfull")
allc_supp <- allc_raw[, !colnames(allc_raw) %in% drop_cols, drop = FALSE]
message("snmC-seq annotation: ", nrow(allc_supp), " cells x ",
        ncol(allc_supp), " columns (dropped: ",
        paste(intersect(drop_cols, colnames(allc_raw)), collapse = ", "), ")")
if (nrow(allc_supp) > 1e5) {
  message("NB this sheet is large; Excel's row limit is 1,048,576 so it fits, ",
          "but the workbook will be slow to open.")
}

# (c) snmC-seq pseudobulk QC, written by src/meth/01b_meth_pseudobulk_qc.R
pbqc_supp     <- read.csv(pbqc_csv, stringsAsFactors = FALSE, check.names = FALSE)
pbqc_grp_supp <- read.csv(pbqc_grp_csv, stringsAsFactors = FALSE, check.names = FALSE)
message("snmC-seq pseudobulk QC: ", nrow(pbqc_supp), " pseudobulks, ",
        nrow(pbqc_grp_supp), " exposure group x cell type combinations")

new_content <- list(
  "Table S2"  = allc_supp,
  "Table S3"  = pbqc_supp,
  "Table S3B" = pbqc_grp_supp,
  "Table S5"  = comp_supp
)

## ---- 3. Rename existing sheets ---------------------------------------------
wb <- openxlsx::loadWorkbook(supp_in)

if (all(names(new_content) %in% names(wb)) &&
    identical(dim(openxlsx::readWorkbook(wb, "Table S5")), dim(comp_supp))) {
  message("Workbook already assembled; nothing to do.")
} else {
  # Rename via a temporary prefix so an old name can never collide with a new
  # one mid-pass (old "Table S2" -> new "Table S3" while an old "Table S3" still
  # exists). Two passes: to temp names, then to final names.
  renames <- final_layout %>%
    dplyr::filter(!is.na(old), old != new)

  tmp_map <- list()
  for (i in seq_len(nrow(renames))) {
    for (suf in c("", LETTERS[1:6])) {
      old_i <- paste0(renames$old[i], suf)
      new_i <- paste0(renames$new[i], suf)
      if (old_i %in% names(wb)) {
        tmp_i <- paste0("__tmp__", new_i)
        openxlsx::renameWorksheet(wb, old_i, tmp_i)
        tmp_map[[tmp_i]] <- new_i
      }
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
    # auto-width is very slow on a sheet with tens of thousands of rows
    if (nrow(dat) < 5000) {
      openxlsx::setColWidths(wb, nm, cols = seq_along(dat), widths = "auto")
    } else {
      openxlsx::setColWidths(wb, nm, cols = seq_along(dat), widths = 18)
    }
    message("Wrote ", nm, " (", nrow(dat), " rows)")
  }

  ## ---- 5. Order the sheets to match final_layout ----------------------------
  want <- unlist(lapply(final_layout$new, function(nm) {
    present <- c(nm, paste0(nm, LETTERS[1:6]))
    present[present %in% names(wb)]
  }))
  other <- setdiff(names(wb), c(want, "Index"))
  ord_names <- c("Index"[("Index" %in% names(wb))], want, other)
  openxlsx::worksheetOrder(wb) <- match(ord_names, names(wb))

  ## ---- 6. Rebuild the Index sheet -------------------------------------------
  # Built from final_layout rather than patched, so the Index cannot drift out of
  # step with the sheet names.
  idx <- data.frame(
    Table = final_layout$new,
    Description = final_layout$description,
    stringsAsFactors = FALSE
  )
  idx <- idx[idx$Table %in% names(wb) |
               paste0(idx$Table, "A") %in% names(wb), , drop = FALSE]

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
  print(as.data.frame(map_df[map_df$changed, c("old", "new")])[
    order(-as.numeric(sub("Table S", "", map_df$old[map_df$changed]))), ])
  message("\nNew tables: S2 = snmC-seq per-cell annotation, ",
          "S3/S3B = snmC-seq pseudobulk QC, ",
          "S5 = cell-type composition statistics.")
  message("NB sub-sheets move with their parent (old S3A/B/C -> S6A/B/C).")
}

#####################################################################
