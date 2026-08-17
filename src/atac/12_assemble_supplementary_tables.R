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
#         (snmC-seq per-cell annotation; committed in the repo, PRE-filter —
#          subset to the analysed cells via cellAnnot_meth.rds, see below)
#   - src/atac/10_4_tex_chromvar_volcano.R -> figures/tex_chromvar_volcano_altius.csv
#         (Tex vs other CD8+ consensus TFBS comparison)
#   - src/meth/01b_meth_pseudobulk_qc.R -> figures/meth_pseudobulk_qc.csv
#         and figures/meth_pseudobulk_qc_by_group.csv
#
# Several tables are inserted, so the numbering shifts more than once. The final
# layout is declared once in `final_layout` and every rename derives from it;
# edit that table if the order changes. Always reads the pre-shift workbook and
# writes a fresh final file, so re-running is safe.
#
# FINAL NUMBERING
# ATAC tables first, in the order the Results use them, then the snmC-seq block,
# then the tables that combine the two modalities.
#
#   S1    Sample metadata of the scATAC dataset
#   S1B   Per-PC covariate association (before/after Harmony batch correction)
#   S2    Cluster-to-cell-type annotation mapping                  (unchanged)
#   S3    Cell-type composition statistics                              [NEW]
#   S4    Variance partition of TF motif activity (S4A/B/C)    (was S3A/B/C)
#   S5    Pairwise differential TF motif activity across cell types (was S4)
#   S6    Differentially accessible genes in CD8+ T-cell clusters   (was S5)
#   S7    Differentially accessible cCREs in CD8+ T-cell clusters   (was S6)
#   S8    Differential gene activity & expression, COVID-19 sev vs ctrl (was S7)
#   S9    Differentially accessible cCREs, COVID-19 severe vs control (was S8)
#   S10   Bulk RNA-seq differential expression, CD14+ monocytes     (was S9)
#   S11   snmC-seq per-cell annotation                                  [NEW]
#   S11B  snmC-seq per-sample summary (donor x timepoint x cell type)   [NEW]
#   S12   snmC-seq pseudobulk QC, per pseudobulk                        [NEW]
#   S12B  snmC-seq pseudobulk QC, per exposure group x cell type        [NEW]
#   S13   Pairwise Wilcoxon (methylTFR vs chromVAR) across cell types (was S10)
#   S14   Differentially methylated cCREs, COVID-19 CD14+ monocytes  (was S11)
#   S15   Tex vs other CD8+ consensus TFBS comparison                     [NEW]
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)  # pivot_wider() for the snmC per-sample cell counts
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
# `old` is the incoming sheet name, NA for new tables
final_layout <- tibble::tribble(
  ~new,          ~old,          ~description,
  "Table S1",    "Table S1",    "Sample metadata of the scATAC dataset",
  "Table S1B",   "Table S1B",   paste(
    "Association of sample-level principal components with known covariates,",
    "before (IterativeLSI) and after (Harmony) batch correction"),
  "Table S2",    "Table S2",    paste(
    "Cluster-to-cell-type annotation mapping: Jaccard similarity between each",
    "graph-based scATAC cluster and predicted scRNA-seq labels, with final",
    "assigned cell type and cluster size"),
  "Table S3",    NA,            paste(
    "Cell-type composition statistics: per-sample cell-type proportions compared",
    "between each exposure group and its matched within-cohort control",
    "(two-sided Wilcoxon rank-sum test, Bonferroni-adjusted)"),
  "Table S4",    "Table S3",    paste(
    "Variance partition of TF motif activity across cell type, exposure and donor",
    "components (S4A: exposure model, S4B: stage/severity model,",
    "S4C: COVID-19 severity)"),
  "Table S5",    "Table S4",    paste(
    "Pairwise differential TF motif activity (Wilcoxon test of chromVAR deviation",
    "scores) across different cell types"),
  "Table S6",    "Table S5",    "List of differentially accessible genes in different clusters within CD8+ T cells",
  "Table S7",    "Table S6",    "List of differentially accessible cCREs in different clusters within CD8+ T cells",
  "Table S8",    "Table S7",    paste(
    "Differential gene activity and gene expression table of protein-coding genes",
    "for COVID-19 severe vs control in CD14+ monocytes"),
  "Table S9",    "Table S8",    "Differentially accessible cCREs for COVID-19 severe vs control in CD14+ monocytes",
  "Table S10",   "Table S9",    paste(
    "Bulk RNA-seq differential expression results for CD14+ monocytes",
    "(filtered: adj p < 0.05 & |log2FC| > 0.5)"),
  "Table S11",   NA,            paste(
    "Per-cell annotation of the snmC-seq dataset: cell identifier, assigned cell",
    "type, donor sex and age, condition, mapping and coverage statistics, and",
    "global mCG / mCH / CCC methylation rates"),
  "Table S11B",  NA,            paste(
    "Per-sample summary of the snmC-seq dataset: cohort, subject, timepoint,",
    "exposure level, donor sex and age, sequencing libraries, and the number of",
    "profiled nuclei in each FACS-sorted cell population, with per-sample means of",
    "the mapping, coverage and methylation-rate metrics"),
  "Table S12",   NA,            paste(
    "Quality-control statistics for each snmC-seq pseudobulk (one per cell type",
    "x sample): donor, exposure group, cells pooled, per-cell sequencing depth",
    "and mapping metrics, CpGs called, CpG coverage and its distribution, global",
    "mCG level, and the CCC rate as the bisulfite non-conversion estimate"),
  "Table S12B",  NA,            paste(
    "snmC-seq pseudobulk quality control summarised per exposure group and cell",
    "type: number of pseudobulks and donors, cells pooled, CpGs called and global",
    "mCG level"),
  "Table S13",   "Table S10",   paste(
    "Pairwise Wilcoxon test between one vs other cell type manner for methylTFR",
    "and chromVAR z-scores"),
  "Table S14",   "Table S11",   "Differentially methylated cCREs in CD14+ monocytes",
  "Table S15",   NA,            paste(
    "Differential TF motif activity between the exhausted (Tex) and the remaining",
    "CD8+ memory T-cell subclusters of the HIV cohort, across 286 consensus TFBS",
    "archetypes: difference in mean chromVAR deviation z-score, Wilcoxon p and",
    "Benjamini-Hochberg adjusted p, the number of subjects in which the change has",
    "the same direction, and whether the archetype was called")
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

# (b) snmC-seq per-cell annotation. Drops the write.csv row-number column and
# allC_FilePathfull (absolute path under /icbb/projects/igunduz/...); the
# relative allC_FilePath is kept.
allc_raw <- read.csv(allc_csv, stringsAsFactors = FALSE, check.names = FALSE)
drop_cols <- c("", "X", "allC_FilePathfull")
allc_supp <- allc_raw[, !colnames(allc_raw) %in% drop_cols, drop = FALSE]
message("snmC-seq annotation as read: ", nrow(allc_supp), " cells x ",
        ncol(allc_supp), " columns (dropped: ",
        paste(intersect(drop_cols, colnames(allc_raw)), collapse = ", "), ")")

# Restrict to the cells actually analysed. The CSV is the pre-filter annotation
# (115 samples, 81,907 cells); the Results report 22,488 cells from the samples
# that overlap the ATAC data, so S11 and S11B must match that set or they
# contradict the sentence that cites them.
#
# The filter chain lives in src/integration/:
#   01_prepare_sampleannot.R  exposure types COVID/FLU/HIV/OP        -> 105 samples
#                             intersect with ATAC sample IDs         ->  39 samples
#                             >= 200 ATAC and >= 50 snmC cells       ->  38 samples, 27,225 cells
#   (7 FACS populations, Other-cell dropped)                         ->  38 samples, 23,981 cells
#   03_quality_check.R        N_valid_sites in [5e5, 4e6]            ->  22,488 cells
#
# The last step cannot be reproduced here: the CSV has no N_valid_sites column.
# So prefer the saved post-QC annotation, and fall back to the reproducible
# subset with a loud warning if it is not reachable.
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
    warning("Expected ", n_cells_expected, " cells to match the Results, got ",
            nrow(allc_supp), ". Check that cellAnnot_meth.rds is the object the ",
            "manuscript numbers came from, and that the Results text matches.")
  }
} else {
  warning("cellAnnot_meth.rds not reachable at ", meth_qc_rds,
          ". Falling back to the filters reproducible from the repo, which stop ",
          "short of the N_valid_sites step, so the cell count will exceed the ",
          n_cells_expected, " reported in the Results. Do not submit this sheet.")
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

message("snmC-seq annotation for S11: ", nrow(allc_supp), " cells, ",
        length(unique(allc_supp$Common_Minimal_Informative_ID)), " samples")
if (nrow(allc_supp) > 1e5) {
  message("NB this sheet is large; Excel's row limit is 1,048,576 so it fits, ",
          "but the workbook will be slow to open.")
}

# (b2) snmC-seq per-sample summary. Reviewer 1 comment 3 asks for the breakdown
# of methylation cells across donors, timepoints and cell subsets; S2 holds one
# row per cell, which does not answer that directly. Aggregated from the same
# object so the two sheets cannot disagree.
#
# Subject / timepoint / exposure level are parsed from the sample identifier:
#   CoV_S_S11_D1          severity (S / nS), subject, day
#   Ctrl_10_M_White_39yo  index, sex, ethnicity, age
#   Flu_S10_D28           subject, day
#   HIV_S1_Pre            subject, stage (Pre / Acu / Cro)
#   OP_S12_High_D28       subject, exposure level, day (or NonFarm)
ct_levels <- c("B-cell", "Monocyte", "NK-cell",
               "Th-Naive", "Th-Mem", "Th-Eff",
               "Tc-Naive", "Tc-Mem", "Tc-Eff", "Other-cell")
# Other-cell is absent once the analysed set is used, so intersect rather than
# require equality
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

# Sex must not contradict itself within a subject. In the annotation as committed,
# four subjects carry two different values across their own timepoints
# (COVID-19 S11 and S17, OP S4 and S7); age is consistent throughout. Reviewer 1
# comment 1 is about exactly this metadata, so this is a warning, not a stop.
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

# (c) snmC-seq pseudobulk QC, written by src/meth/01b_meth_pseudobulk_qc.R
pbqc_supp     <- read.csv(pbqc_csv, stringsAsFactors = FALSE, check.names = FALSE)
pbqc_grp_supp <- read.csv(pbqc_grp_csv, stringsAsFactors = FALSE, check.names = FALSE)
message("snmC-seq pseudobulk QC: ", nrow(pbqc_supp), " pseudobulks, ",
        nrow(pbqc_grp_supp), " exposure group x cell type combinations")

# (d) Tex vs other CD8+ consensus TFBS comparison, written by
# src/atac/10_4_tex_chromvar_volcano.R and plotted in Supplementary Figure 3F.
# Selection is on effect size and agreement across subjects, not on padj, so both
# are reported and `Called` records the outcome.
tex_csv  <- file.path(fig_dir, "tex_chromvar_volcano_altius.csv")
stopifnot(file.exists(tex_csv))
tex_supp <- read.csv(tex_csv, stringsAsFactors = FALSE) %>%
  dplyr::transmute(
    `Consensus TFBS`                        = motif,
    `Delta z (Tex - other CD8+)`            = signif(zdiff, 3),
    `p (Wilcoxon)`                          = signif(p, 3),
    `p adjusted (Benjamini-Hochberg)`       = signif(padj, 3),
    `Subjects with same direction (of 4)`   = donors_consistent,
    `Called`                                = dplyr::recode(group,
                                                `Tex` = "higher in Tex",
                                                `Other` = "higher in other CD8+",
                                                `NO` = "not called")
  ) %>%
  dplyr::arrange(dplyr::desc(`Delta z (Tex - other CD8+)`))
message("Tex consensus TFBS comparison: ", nrow(tex_supp), " archetypes, ",
        sum(tex_supp$Called != "not called"), " called")

new_content <- list(
  "Table S3"   = comp_supp,
  "Table S11"  = allc_supp,
  "Table S11B" = allc_sample_supp,
  "Table S12"  = pbqc_supp,
  "Table S12B" = pbqc_grp_supp,
  "Table S15"  = tex_supp
)

## ---- 3. Rename existing sheets ---------------------------------------------
wb <- openxlsx::loadWorkbook(supp_in)

if (all(names(new_content) %in% names(wb)) &&
    identical(dim(openxlsx::readWorkbook(wb, "Table S3")), dim(comp_supp))) {
  message("Workbook already assembled; nothing to do.")
} else {
  # temp prefix first, so an old name cannot collide with a new one mid-pass
  # (old S3 -> new S4 while an old S4 still exists, and old S10/S11 -> S13/S14
  # must vacate before the new S11/S12 sheets are inserted)
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
    # auto-width is slow on tens of thousands of rows
    if (nrow(dat) < 5000) {
      openxlsx::setColWidths(wb, nm, cols = seq_along(dat), widths = "auto")
    } else {
      openxlsx::setColWidths(wb, nm, cols = seq_along(dat), widths = 18)
    }
    message("Wrote ", nm, " (", nrow(dat), " rows)")
  }

  ## ---- 5. Order the sheets to match final_layout ----------------------------
  # unique() matters: a sub-sheet is reached twice, once by the suffix expansion
  # of its parent (S11 -> S11B) and once by its own row in final_layout. Without
  # it, ord_names is longer than the workbook and worksheetOrder<- rejects it.
  want <- unique(unlist(lapply(final_layout$new, function(nm) {
    present <- c(nm, paste0(nm, LETTERS[1:6]))
    present[present %in% names(wb)]
  })))
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
  message("\nNew tables: S3 = cell-type composition statistics, ",
          "S11 = snmC-seq per-cell annotation, ",
          "S11B = snmC-seq per-sample summary, ",
          "S12/S12B = snmC-seq pseudobulk QC, ",
          "S15 = Tex vs other CD8+ consensus TFBS comparison.")
  message("NB sub-sheets move with their parent (old S3A/B/C -> S4A/B/C).")
}

#####################################################################
