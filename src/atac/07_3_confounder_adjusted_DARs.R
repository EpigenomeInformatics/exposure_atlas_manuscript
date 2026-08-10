#!/usr/bin/env Rscript

#####################################################################
# 07_3_confounder_adjusted_DARs.R
# created on 2026-08-08 by Irem B. Gunduz
#
# Reviewer response (R1.2 / R2 technical covariates; supervisor request):
# re-call differentially accessible regions with per-sample technical
# covariates included in the differential model, and show that the reported
# results are qualitatively unchanged.
#
# Approach
#  - Load the EXISTING processed ChrAccR DsATAC objects (no re-import, no
#    re-normalisation) so the adjusted and unadjusted calls differ ONLY in the
#    design matrix.
#  - Attach per-sample QC metrics (mean TSS enrichment, mean log10 fragments,
#    mean FRIP, cell number) computed from the ArchR project.
#  - Re-run ChrAccR's differential module with those metrics as adjustment
#    columns, then compare adjusted vs unadjusted DARs: counts, overlap,
#    direction concordance and effect-size correlation.
#
# NOTE ON POWER: these comparisons have few samples per group (e.g. n = 6
# severe vs 7 control). Every added covariate costs a degree of freedom, so we
# adjust for a SMALL set of metrics (TSS and FRIP by default) rather than all
# of them. The point of the analysis is to show the conclusions are robust,
# not to maximise the number of DARs.
#####################################################################

suppressPackageStartupMessages({
  library(ArchR)
  library(ChrAccR)
  library(dplyr)
  library(ggplot2)
})
set.seed(12)

repo_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript"
source(file.path(repo_dir, "utils/chraccr_plots.R")) # cutL0.5fc2Padj05
out_dir <- "/icbb/projects/igunduz/finalize_echo_050824/confounder_adjusted/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# adjustment columns to add on top of any existing ones
qc_adj_cols <- c("mean_TSS", "mean_FRIP")

# Normalise a sample identifier for joining ChrAccR sample ids to the ArchR QC
# table: drop any directory, the .tsv.gz extension, and a trailing "_fragments",
# so "hiv6_fragments.tsv.gz", "hiv6_fragments" and "hiv6" all collapse to "hiv6".
# (The old code stripped only ".tsv.gz", leaving "hiv6_fragments" != "hiv6",
# which is why every sample failed to match.)
norm_key <- function(x) {
  x <- basename(as.character(x))
  x <- sub("\\.tsv\\.gz.*$", "", x)   # drop .tsv.gz and anything after it
  x <- sub("_fragments$", "", x)      # drop a trailing _fragments
  x
}

## ---- 1. Per-sample QC metrics from the ArchR project ------------------------
archr_dir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/"
project <- ArchR::loadArchRProject(archr_dir, showLogo = FALSE)

qc_by_sample <- getCellColData(project,
  select = c("Sample", "TSSEnrichment", "nFrags", "FRIP")
) %>%
  as.data.frame() %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarise(
    n_cells = dplyr::n(),
    mean_TSS = mean(TSSEnrichment, na.rm = TRUE),
    mean_log10_nFrags = mean(log10(nFrags), na.rm = TRUE),
    mean_FRIP = mean(FRIP, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(sample_key = norm_key(Sample))
rm(project)
gc()

## ---- 2. Helper: run one comparison, adjusted vs unadjusted ------------------
# anaDir_existing : ChrAccR run directory containing <cell>/data/dsATAC_processed
# comp            : "GRP1 vs GRP2 [sample_exposure_group]"
# extra_adj       : adjustment columns already used in the original run
run_adjusted <- function(anaDir_existing, cell, comp, extra_adj = character(0)) {
  ds <- ChrAccR::loadDsAcc(file.path(anaDir_existing, cell, "data", "dsATAC_processed"))
  sa <- ChrAccR::getSampleAnnot(ds)

  # match ChrAccR sample ids to the ArchR QC table. ChrAccR sample ids derive
  # from the per-cell-type BED filenames; normalise both sides to a common key
  # before matching.
  sa_keys <- norm_key(rownames(sa))
  if (all(is.na(match(sa_keys, qc_by_sample$sample_key)))) {
    sa_keys <- norm_key(sa[[1]])   # fall back to the first annotation column
  }
  idx <- match(sa_keys, qc_by_sample$sample_key)
  # substring fallback for any residual mismatch (e.g. cohort prefixes)
  for (ii in which(is.na(idx))) {
    hit <- which(vapply(qc_by_sample$sample_key,
      function(k) grepl(k, sa_keys[ii], fixed = TRUE) ||
                  grepl(sa_keys[ii], k, fixed = TRUE), logical(1)))
    if (length(hit) == 1) idx[ii] <- hit
  }
  if (any(is.na(idx))) {
    warning(cell, " / ", comp, ": ", sum(is.na(idx)),
      " sample(s) could not be matched to QC metrics; they will be dropped by the model")
    message("  unmatched sample keys: ", paste(head(unique(sa_keys[is.na(idx)]), 10), collapse = ", "))
    message("  available QC keys: ", paste(head(qc_by_sample$sample_key, 20), collapse = ", "))
  } else {
    message("  matched all ", length(idx), " samples to QC metrics")
  }
  for (cc in c(qc_adj_cols, "mean_log10_nFrags", "n_cells")) {
    sa[[cc]] <- qc_by_sample[[cc]][idx]
  }
  # write the augmented annotation back onto the object
  ds@sampleAnnot <- sa

  # DESeq2-based differential needs raw integer counts. The processed DsATAC can
  # carry normalized/transformed (non-integer) counts, which makes
  # run_atac_differential fail with "some values in assay are not integers".
  # Loading a raw/pre-normalization object would be cleaner; as a safeguard we
  # round any non-integer count matrix back to integers so the model can run.
  ct_slots <- tryCatch(names(ds@counts), error = function(e) character(0))
  for (rt in ct_slots) {
    cm <- tryCatch(as.matrix(ds@counts[[rt]]), error = function(e) NULL)
    if (is.null(cm)) next
    if (any(abs(cm - round(cm)) > 1e-6, na.rm = TRUE)) {
      message("  region type '", rt, "': non-integer counts detected, rounding for DESeq2")
      ds@counts[[rt]] <- round(cm)
    }
  }

  # ---- adjusted run, into a fresh directory ----
  tag <- gsub("[^A-Za-z0-9]+", "_", sub(" \\[.*", "", comp))
  ana_adj <- file.path(out_dir, paste0(cell, "__", tag, "__adjusted"))
  if (!dir.exists(ana_adj)) dir.create(ana_adj, recursive = TRUE)

  setConfigElement("differentialColumns", c("sample_exposure_group"))
  setConfigElement("differentialCompNames", comp)
  setConfigElement("differentialAdjColumns", unique(c(extra_adj, qc_adj_cols)))
  setConfigElement("differentialCutoffL2FC", 0.5)
  setConfigElement("filteringSexChroms", TRUE)
  run_atac_differential(ds, ana_adj)

  # ---- read both tables and compare ----
  read_diff <- function(dir) {
    f <- list.files(dir, pattern = "diffTab.*archrPeaks.*\\.tsv$",
      recursive = TRUE, full.names = TRUE
    )
    if (length(f) == 0) {
      f <- list.files(dir, pattern = "diffTab.*\\.tsv$", recursive = TRUE, full.names = TRUE)
    }
    stopifnot(length(f) >= 1)
    dm <- read.delim(f[1])
    isDiff <- cutL0.5fc2Padj05(dm[, c("log2FoldChange", "padj")])
    isDiff[is.na(isDiff)] <- FALSE
    data.frame(
      id = paste0(dm$chrom, ":", dm$chromStart + 1, "-", dm$chromEnd),
      log2FC = dm$log2FoldChange, padj = dm$padj, isDiff = isDiff,
      stringsAsFactors = FALSE
    )
  }

  adj <- read_diff(ana_adj)
  # unadjusted: the comparison as originally run in the existing analysis dir
  unadj_all <- list.files(file.path(anaDir_existing, cell, "reports", "differential_data"),
    pattern = "diffTab.*archrPeaks.*\\.tsv$", full.names = TRUE
  )
  ct <- readRDS(file.path(anaDir_existing, cell, "reports", "differential_data", "comparisonTable.rds"))
  want <- sub(" \\[.*", "", comp)
  i <- which(paste0(ct$grp1, " vs ", ct$grp2) == want)
  stopifnot(length(i) == 1, length(unadj_all) >= i)
  dm <- read.delim(unadj_all[i])
  isD <- cutL0.5fc2Padj05(dm[, c("log2FoldChange", "padj")])
  isD[is.na(isD)] <- FALSE
  unadj <- data.frame(
    id = paste0(dm$chrom, ":", dm$chromStart + 1, "-", dm$chromEnd),
    log2FC = dm$log2FoldChange, padj = dm$padj, isDiff = isD, stringsAsFactors = FALSE
  )

  m <- merge(unadj, adj, by = "id", suffixes = c("_unadj", "_adj"))
  a <- m$id[m$isDiff_unadj]
  b <- m$id[m$isDiff_adj]
  shared <- intersect(a, b)
  conc <- if (length(shared)) {
    mean(sign(m$log2FC_unadj[match(shared, m$id)]) ==
      sign(m$log2FC_adj[match(shared, m$id)])) * 100
  } else NA_real_

  data.frame(
    cell = cell, comparison = want,
    DARs_unadjusted = length(a), DARs_adjusted = length(b),
    shared = length(shared),
    recovered_pct = ifelse(length(a) > 0, round(100 * length(shared) / length(a), 1), NA),
    sign_concordance_pct = round(conc, 1),
    lfc_pearson_all = round(cor(m$log2FC_unadj, m$log2FC_adj, use = "complete.obs"), 3),
    stringsAsFactors = FALSE
  )
}

## ---- 3. Run for the comparisons the manuscript's conclusions rest on --------
covid_dir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/"
other_dir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_2023-10-02/"

jobs <- list(
  list(dir = covid_dir, cell = "Mono_CD14",
       comp = "C19_sev vs C19_ctrl [sample_exposure_group]", adj = "processing_date"),
  list(dir = covid_dir, cell = "Mono_CD14",
       comp = "C19_mod vs C19_ctrl [sample_exposure_group]", adj = "processing_date"),
  list(dir = other_dir, cell = "T_mem_CD8",
       comp = "HIV_ctrl vs HIV_acu [sample_exposure_group]", adj = character(0)),
  list(dir = other_dir, cell = "T_mem_CD8",
       comp = "HIV_ctrl vs HIV_chr [sample_exposure_group]", adj = character(0))
)

res <- do.call(rbind, lapply(jobs, function(j) {
  message("=== ", j$cell, " | ", j$comp)
  tryCatch(run_adjusted(j$dir, j$cell, j$comp, j$adj),
    error = function(e) {
      message("  FAILED: ", conditionMessage(e))
      NULL
    }
  )
}))

print(res)
write.csv(res, file.path(out_dir, "confounder_adjusted_DAR_summary.csv"), row.names = FALSE)

## ---- 4. Summary panel -------------------------------------------------------
if (!is.null(res) && nrow(res) > 0) {
  plot_df <- res %>%
    dplyr::select(comparison, cell, DARs_unadjusted, DARs_adjusted) %>%
    tidyr::pivot_longer(c(DARs_unadjusted, DARs_adjusted),
      names_to = "model", values_to = "n_DARs"
    ) %>%
    dplyr::mutate(
      model = ifelse(model == "DARs_unadjusted", "Unadjusted", "QC-adjusted"),
      label = paste0(cell, "\n", comparison)
    )
  p <- ggplot(plot_df, aes(x = label, y = n_DARs, fill = model)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65, colour = "grey20") +
    geom_text(aes(label = n_DARs),
      position = position_dodge(width = 0.7), vjust = -0.3, size = 3
    ) +
    scale_fill_manual(values = c("Unadjusted" = "#3C5488", "QC-adjusted" = "#E64B35")) +
    labs(
      title = "Differential accessibility with and without technical-covariate adjustment",
      subtitle = paste0("Adjusted model adds ", paste(qc_adj_cols, collapse = " + "),
        " to the differential design"),
      x = NULL, y = "Number of DARs", fill = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top", plot.subtitle = element_text(size = 9))
  ggsave(file.path(out_dir, "confounder_adjusted_DARs.pdf"), p, width = 7.5, height = 5)
}

message("Done. Results in ", out_dir)

#####################################################################
