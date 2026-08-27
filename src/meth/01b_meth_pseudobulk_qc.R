#!/usr/bin/env Rscript

#####################################################################
# 01b_meth_pseudobulk_qc.R
# created on 2026-08-12 by Irem B. Gunduz
# QC per snmC-seq pseudobulk (supplementary table)
#
# The methylation analyses all run on the pseudobulks built in
# 01_meth_pseudobulks.R (one per cell type x sample), but nothing reports how
# much data sits behind each one. Measured here from the bedGraphs.
#
# !! The bedGraph 5th column is not one convention. In aggregateALLCSamples()
# the chunked branch (>50 cells pooled) ends with mutate(cov = cov - mc), so it
# holds UNMETHYLATED counts; the unchunked branch (<=50) writes total coverage.
# Here that is 844 of 896 pseudobulks vs 52, i.e. both conventions in the same
# directory. This script derives which from the cell count and records it in
# `coverage_convention`. The real fix is in createPseudoBulks.R, but that means
# rebuilding the pseudobulks.
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})
set.seed(12)

repo_dir <- "/scratch/icbb/igunduz/irem_github/exposure_atlas_manuscript"
fig_dir  <- file.path(repo_dir, "figures")
annot    <- file.path(repo_dir, "sample_annots")
pb_dir   <- "/icbb/projects/igunduz/DARPA/data/pseudoBulks/perSample"

if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)
stopifnot(dir.exists(pb_dir))

# branch threshold in createPseudoBulks()
chunk_threshold <- 50

# reported as fraction of called CpGs reaching Nx
cov_thresholds <- c(1, 5, 10)

## ---- 1. Per-cell annotation -------------------------------------------------
# as in 01_meth_pseudobulks.R: Other-cell is not pseudobulked
sampleAnnot <- data.table::fread(
  file.path(annot, "allc_sample_annot_final.csv")
) %>%
  dplyr::select(-dplyr::any_of(c("V1", "X"))) %>%
  dplyr::filter(cell_type != "Other-cell")

message("Annotation: ", nrow(sampleAnnot), " cells, ",
  dplyr::n_distinct(sampleAnnot$cell_type), " cell types, ",
  dplyr::n_distinct(sampleAnnot$Common_Minimal_Informative_ID), " samples")

# map the methylation condition labels onto the ATAC exposure_group names
condition_to_group <- c(
  CommercialControl_healthy = "C19_ctrl",
  COVID_mild                = "C19_mild",
  COVID_severe              = "C19_sev",
  HIV_pre                   = "HIV_ctrl",
  HIV_acute                 = "HIV_acu",
  HIV_chronic               = "HIV_chr",
  FLU_healthy               = "Influenza_ctrl",
  FLU_day30                 = "Influenza_d28",
  OP_low                    = "OP_low",
  OP_medium                 = "OP_med",
  OP_high                   = "OP_high"
)
# "FLU_day30" is the day-28 mislabel also corrected in Table S1. Raw value kept.

cell_summary <- sampleAnnot %>%
  dplyr::group_by(cell_type, Common_Minimal_Informative_ID) %>%
  dplyr::summarise(
    donor_id            = paste(sort(unique(Salk_ID)), collapse = ";"),
    condition_raw       = paste(sort(unique(condition)), collapse = ";"),
    sex                 = paste(sort(unique(Sex)), collapse = ";"),
    age                 = paste(sort(unique(age)), collapse = ";"),
    n_cells             = dplyr::n(),
    median_unique_reads = stats::median(TotalUniqueMappedReads, na.rm = TRUE),
    min_unique_reads    = min(TotalUniqueMappedReads, na.rm = TRUE),
    mean_mapping_ratio  = mean(TotalMappedRatio, na.rm = TRUE),
    mean_genome_cov     = mean(genome_cov, na.rm = TRUE),
    mean_CG_rate        = mean(CG_Rate, na.rm = TRUE),
    mean_CH_rate        = mean(CH_Rate, na.rm = TRUE),
    mean_CCC_rate       = mean(CCC_Rate, na.rm = TRUE),
    max_CCC_rate        = max(CCC_Rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    exposure_group = unname(condition_to_group[condition_raw]),
    bedGraph = paste0(cell_type, "_", Common_Minimal_Informative_ID, ".bedGraph")
  )

if (any(is.na(cell_summary$exposure_group))) {
  message("Unmapped condition label(s): ",
    paste(unique(cell_summary$condition_raw[is.na(cell_summary$exposure_group)]),
      collapse = ", "))
}

## ---- 2. Measure each pseudobulk --------------------------------------------
# headerless from makeGRanges(): seqnames, start, strand, mc, cov
# cov is unmethylated or total depending on the branch (see header)
read_pseudobulk_qc <- function(path, n_cells) {
  if (!file.exists(path)) {
    return(data.frame(file_found = FALSE, stringsAsFactors = FALSE))
  }
  dt <- data.table::fread(path, header = FALSE,
    col.names = c("seqnames", "start", "strand", "mc", "cov"),
    showProgress = FALSE)
  if (nrow(dt) == 0) {
    return(data.frame(file_found = TRUE, n_CpG_called = 0L, stringsAsFactors = FALSE))
  }

  chunked <- n_cells > chunk_threshold
  total_cov <- if (chunked) dt$mc + dt$cov else dt$cov
  convention <- if (chunked) "mc + unmethylated" else "mc, total coverage"

  qc <- data.frame(
    file_found          = TRUE,
    coverage_convention = convention,
    n_CpG_called        = nrow(dt),
    total_CpG_coverage  = sum(total_cov, na.rm = TRUE),
    mean_CpG_coverage   = mean(total_cov, na.rm = TRUE),
    median_CpG_coverage = stats::median(total_cov, na.rm = TRUE),
    global_mCG          = sum(dt$mc, na.rm = TRUE) / sum(total_cov, na.rm = TRUE),
    n_chromosomes       = dplyr::n_distinct(dt$seqnames),
    stringsAsFactors = FALSE
  )
  for (thr in cov_thresholds) {
    qc[[paste0("frac_CpG_ge_", thr, "x")]] <- mean(total_cov >= thr, na.rm = TRUE)
  }
  qc
}

message("Measuring ", nrow(cell_summary), " pseudobulk file(s) in ", pb_dir)
qc_list <- vector("list", nrow(cell_summary))
for (i in seq_len(nrow(cell_summary))) {
  if (i %% 50 == 0) message("  ", i, "/", nrow(cell_summary))
  qc_list[[i]] <- read_pseudobulk_qc(
    file.path(pb_dir, cell_summary$bedGraph[i]), cell_summary$n_cells[i]
  )
}
qc_df <- dplyr::bind_rows(qc_list)
pb_qc <- dplyr::bind_cols(cell_summary, qc_df)

n_missing <- sum(!pb_qc$file_found)
if (n_missing > 0) {
  message("!! ", n_missing, " pseudobulk file(s) not found on disk:")
  print(head(pb_qc$bedGraph[!pb_qc$file_found], 20))
}
message("Convention split: ",
  paste(sprintf("%s = %d", names(table(pb_qc$coverage_convention)),
    as.integer(table(pb_qc$coverage_convention))), collapse = ", "))

## ---- 3. Consistency checks --------------------------------------------------
# these run every time
message("\nChecks:")
message("  cells in the pseudobulk table vs the annotation: ",
  sum(pb_qc$n_cells), " vs ", nrow(sampleAnnot),
  if (sum(pb_qc$n_cells) == nrow(sampleAnnot)) "  OK" else "  MISMATCH")

tiny <- pb_qc[pb_qc$n_cells < 10 & pb_qc$file_found, ]
message("  pseudobulks pooling fewer than 10 cells: ", nrow(tiny))
if (nrow(tiny) > 0) {
  print(as.data.frame(tiny[order(tiny$n_cells),
    c("cell_type", "Common_Minimal_Informative_ID", "exposure_group",
      "n_cells", "n_CpG_called")]))
}

bad_mcg <- pb_qc$file_found & !is.na(pb_qc$global_mCG) &
  (pb_qc$global_mCG < 0 | pb_qc$global_mCG > 1)
if (any(bad_mcg)) {
  message("  !! global mCG outside [0, 1] for ", sum(bad_mcg),
    " pseudobulk(s) -- the coverage convention for those files is wrong")
  print(as.data.frame(pb_qc[bad_mcg,
    c("bedGraph", "n_cells", "coverage_convention", "global_mCG")]))
} else {
  message("  global mCG within [0, 1] for every pseudobulk  OK")
}

## ---- 4. Publication-ready table ---------------------------------------------
pb_supp <- pb_qc %>%
  dplyr::filter(file_found) %>%
  dplyr::transmute(
    `Cell type`                      = cell_type,
    `Sample ID`                      = Common_Minimal_Informative_ID,
    `Donor ID`                       = donor_id,
    `Exposure group`                 = exposure_group,
    `Condition (raw label)`          = condition_raw,
    Sex                              = sex,
    Age                              = age,
    `Cells pooled`                   = n_cells,
    `Median unique mapped reads per cell` = round(median_unique_reads),
    `Minimum unique mapped reads per cell` = round(min_unique_reads),
    `Mean mapping ratio (%)`         = round(mean_mapping_ratio, 2),
    `Mean genome coverage per cell`  = signif(mean_genome_cov, 3),
    `CpGs called`                    = n_CpG_called,
    `Total CpG coverage`             = total_CpG_coverage,
    `Mean coverage per CpG`          = round(mean_CpG_coverage, 2),
    `Median coverage per CpG`        = median_CpG_coverage,
    `CpGs at >=1x (%)`               = round(100 * frac_CpG_ge_1x, 1),
    `CpGs at >=5x (%)`               = round(100 * frac_CpG_ge_5x, 1),
    `CpGs at >=10x (%)`              = round(100 * frac_CpG_ge_10x, 1),
    `Global mCG (%)`                 = round(100 * global_mCG, 2),
    `Mean mCH per cell (%)`          = round(100 * mean_CH_rate, 2),
    `Mean CCC rate per cell (%)`     = round(100 * mean_CCC_rate, 3),
    `Max CCC rate per cell (%)`      = round(100 * max_CCC_rate, 3),
    `Coverage convention`            = coverage_convention
  ) %>%
  dplyr::arrange(`Cell type`, `Exposure group`, `Sample ID`)

write.csv(pb_supp, file.path(fig_dir, "meth_pseudobulk_qc.csv"), row.names = FALSE)
message("\nWrote ", nrow(pb_supp), " rows to figures/meth_pseudobulk_qc.csv")

## ---- 5. Per exposure-group x cell-type summary ------------------------------
# the per-pseudobulk table is for checking; this is for the methods sentence
group_supp <- pb_qc %>%
  dplyr::filter(file_found) %>%
  dplyr::group_by(cell_type, exposure_group) %>%
  dplyr::summarise(
    `Pseudobulks`            = dplyr::n(),
    `Donors`                 = dplyr::n_distinct(donor_id),
    `Cells (total)`          = sum(n_cells),
    `Cells per pseudobulk (median)` = stats::median(n_cells),
    `Cells per pseudobulk (min)`    = min(n_cells),
    `Cells per pseudobulk (max)`    = max(n_cells),
    `CpGs called (median)`   = round(stats::median(n_CpG_called)),
    `CpGs called (min)`      = min(n_CpG_called),
    `Mean coverage per CpG (median)` = round(stats::median(mean_CpG_coverage), 2),
    `Global mCG (%) (median)` = round(100 * stats::median(global_mCG), 2),
    `Global mCG (%) (range)`  = sprintf("%.2f-%.2f",
      100 * min(global_mCG), 100 * max(global_mCG)),
    .groups = "drop"
  ) %>%
  dplyr::rename(`Cell type` = cell_type, `Exposure group` = exposure_group) %>%
  dplyr::arrange(`Cell type`, `Exposure group`)

write.csv(group_supp, file.path(fig_dir, "meth_pseudobulk_qc_by_group.csv"),
  row.names = FALSE)
message("Wrote ", nrow(group_supp),
  " rows to figures/meth_pseudobulk_qc_by_group.csv")
print(as.data.frame(group_supp))

message("\nBoth tables are picked up by src/atac/12_assemble_supplementary_tables.R")

#####################################################################
