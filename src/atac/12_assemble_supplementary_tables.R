#!/usr/bin/env Rscript

#####################################################################
# 12_assemble_supplementary_tables.R
# Single source of truth for the supplementary-table workbook + numbering.
#
# Run LAST, after the scripts that produce the individual tables:
#   - 01_v2_quality_control.R  -> All_Supplementary_Tables_updated.xlsx
#         (Table S1 with the recovered metadata columns, and Table S1B)
#   - 04_2_cellprops.R         -> figures/cell_prop_wilcox_bonferroni.csv
#         (cell-type composition statistics, added here as the new Table S3)
#
# This script reads the updated workbook (pre-shift) and the composition CSV,
# inserts the composition table as Table S3, shifts every table from S3 upward
# by one, rebuilds the Index, and writes All_Supplementary_Tables_final.xlsx.
# It always reads the pre-shift workbook and writes a fresh final file, so it is
# idempotent (safe to re-run).
#
# FINAL NUMBERING (edit here if the order changes):
#   S1   Sample metadata of the scATAC dataset
#   S1B  Per-PC covariate association (before/after Harmony batch correction)
#   S2   Cluster-to-cell-type annotation mapping
#   S3   Cell-type composition statistics (per-sample Wilcoxon, Bonferroni)   [NEW]
#   S4   Variance partition of TF motif activity (cell type / exposure / donor)   (was S3)
#   S5   Pairwise differential TF motif activity across cell types               (was S4)
#   S6   Differentially accessible genes in CD8+ T-cell clusters                 (was S5)
#   S7   Differentially accessible cCREs in CD8+ T-cell clusters                 (was S6)
#   S8   Differential gene activity & expression, COVID-19 severe vs control     (was S7)
#   S9   Differentially accessible cCREs, COVID-19 severe vs control             (was S8)
#   S10  Differential gene expression (scRNA-seq), COVID-19 CD14+ monocytes      (was S9)
#   S11  Pairwise Wilcoxon (methylTFR vs chromVAR z-score) across cell types     (was S10)
#   S12  Differentially methylated cCREs, COVID-19 CD14+ monocytes               (was S11)
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(openxlsx)
})

repo_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript"
fig_dir  <- file.path(repo_dir, "figures")
annot    <- file.path(repo_dir, "sample_annots")

# Inputs: the pre-shift workbook (S1/S1B added by 01_v2) and the composition CSV.
supp_updated <- file.path(annot, "All_Supplementary_Tables_updated.xlsx")
supp_master  <- file.path(annot, "All_Supplementary_Tables.xlsx")
supp_in      <- if (file.exists(supp_updated)) supp_updated else supp_master
supp_out     <- file.path(annot, "All_Supplementary_Tables_final.xlsx")
comp_csv     <- file.path(fig_dir, "cell_prop_wilcox_bonferroni.csv")

stopifnot(file.exists(supp_in), file.exists(comp_csv))
message("Reading workbook: ", supp_in)

# ---- publication-ready composition-statistics table (new Table S3) -----------
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

wb <- openxlsx::loadWorkbook(supp_in)

if ("Table S3" %in% names(wb) &&
    identical(dim(openxlsx::readWorkbook(wb, "Table S3")), dim(comp_supp))) {
  message("Composition table already inserted; nothing to do.")
} else {
  # ---- shift every numbered sheet from S3 upward by one -----------------------
  # Handles single sheets (Table S7) and lettered sub-sheets (Table S3A/B/C).
  # Rename from the highest number downwards so names never collide.
  for (n in 11:3) {
    for (suf in c("", LETTERS[1:6])) {
      old <- paste0("Table S", n, suf)
      if (old %in% names(wb)) {
        openxlsx::renameWorksheet(wb, old, paste0("Table S", n + 1, suf))
      }
    }
  }

  # ---- insert the new Table S3 and place it directly after Table S2 -----------
  openxlsx::addWorksheet(wb, "Table S3")
  openxlsx::writeData(wb, "Table S3", comp_supp,
    withFilter  = TRUE,
    headerStyle = openxlsx::createStyle(textDecoration = "bold")
  )
  openxlsx::setColWidths(wb, "Table S3", cols = seq_along(comp_supp), widths = "auto")

  nm        <- names(wb)
  pos_new   <- which(nm == "Table S3")
  pos_after <- which(nm == "Table S2")
  ord       <- append(seq_along(nm)[-pos_new], pos_new, after = pos_after)
  openxlsx::worksheetOrder(wb) <- ord

  # ---- rebuild the Index sheet ------------------------------------------------
  # Shift existing labels from S3 up (high -> low so nothing double-shifts),
  # then insert the new S3 row after S2.
  shift_label <- function(x) {
    for (n in 11:3) {
      x <- gsub(paste0("Table S", n, "([A-Z]?)"),
                paste0("Table S", n + 1, "\\1"), x)
    }
    x
  }
  if ("Index" %in% names(wb)) {
    idx <- openxlsx::readWorkbook(wb, sheet = "Index")
    idx[[1]] <- vapply(as.character(idx[[1]]), shift_label, character(1), USE.NAMES = FALSE)
    new_row  <- idx[1, ]
    new_row[1, ] <- NA
    new_row[[1]] <- "Table S3"
    new_row[[2]] <- paste(
      "Cell-type composition statistics: per-sample cell-type proportions compared",
      "between each exposure group and its matched within-cohort control",
      "(two-sided Wilcoxon rank-sum test, Bonferroni-adjusted)"
    )
    after <- which(idx[[1]] == "Table S2")
    idx   <- rbind(idx[seq_len(after), ], new_row, idx[-seq_len(after), ])

    openxlsx::removeWorksheet(wb, "Index")
    openxlsx::addWorksheet(wb, "Index")
    openxlsx::writeData(wb, "Index", idx,
      headerStyle = openxlsx::createStyle(textDecoration = "bold"))
    openxlsx::setColWidths(wb, "Index", cols = 1:2, widths = c(14, 110))
    openxlsx::worksheetOrder(wb) <- c(
      which(names(wb) == "Index"),
      setdiff(seq_along(names(wb)), which(names(wb) == "Index"))
    )
  } else {
    message("No 'Index' sheet found; skipping Index rebuild.")
  }

  openxlsx::saveWorkbook(wb, supp_out, overwrite = TRUE)
  message("Wrote ", supp_out)
  message("APPLY THIS RENUMBERING TO THE MANUSCRIPT TEXT:")
  print(data.frame(
    old = c("Table S3 (and S3A/B/C)", paste0("Table S", 4:11)),
    new = c("Table S4 (and S4A/B/C)", paste0("Table S", 5:12))
  ))
  message("New Table S3 = cell-type composition statistics.")
}
