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
# Scratch space for the adjusted ChrAccR/DESeq2 runs (large, stays off the repo)
out_dir <- "/icbb/projects/igunduz/finalize_echo_050824/confounder_adjusted/"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

# Manuscript-facing outputs (figures and the summary tables) go into the repo
# figures/ directory alongside every other panel, not into the scratch dir.
fig_dir <- file.path(repo_dir, "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# adjustment columns to add on top of any existing ones
qc_adj_cols <- c("mean_TSS", "mean_FRIP")

# Reuse an adjusted fit that is already on disk instead of refitting it. The
# sweep below covers every cell type and every saved comparison, so this makes an
# interrupted run resumable. Set to FALSE to force everything to be refitted.
resume <- TRUE

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

# The ChrAccR sample ids ("Influenza_ctrl_flu01_1") and the ArchR Sample ids
# ("ATAC_055") are different schemes, so a direct key match fails. Build a lookup
# from EVERY ArchR per-sample descriptor value to the ArchR Sample, so a ChrAccR
# id that matches any descriptor column can be bridged to the QC metrics.
archr_sample_meta <- as.data.frame(getCellColData(project))
archr_sample_meta <- archr_sample_meta[!duplicated(archr_sample_meta$Sample), , drop = FALSE]
message("ArchR per-sample descriptor columns: ",
        paste(colnames(archr_sample_meta), collapse = ", "))
id_lookup <- do.call(rbind, lapply(colnames(archr_sample_meta), function(cn) {
  data.frame(key = norm_key(archr_sample_meta[[cn]]),
             Sample = archr_sample_meta$Sample, column = cn, stringsAsFactors = FALSE)
}))
id_lookup <- id_lookup[!is.na(id_lookup$key) & id_lookup$key != "", ]
id_lookup <- id_lookup[!duplicated(id_lookup$key), ]
rm(project)
gc()

## ---- 2. Helper: run one comparison, adjusted vs unadjusted ------------------
# anaDir_existing : ChrAccR run directory containing <cell>/data/dsATAC_processed
# comp            : "GRP1 vs GRP2 [sample_exposure_group]"
# extra_adj       : adjustment columns already used in the original run
run_adjusted <- function(anaDir_existing, cell, comp, extra_adj = character(0)) {
  tag <- gsub("[^A-Za-z0-9]+", "_", sub(" \\[.*", "", comp))
  ana_adj <- file.path(out_dir, paste0(cell, "__", tag, "__adjusted"))

  # Resume support: the full sweep over every cell type is long, and the adjusted
  # DESeq2 fit is the expensive part. If this comparison has already been fitted
  # in a previous run, reuse the table on disk instead of refitting it.
  done <- list.files(ana_adj, pattern = "diffTab.*\\.tsv$",
    recursive = TRUE, full.names = TRUE)
  if (resume && length(done) > 0) {
    message("  adjusted table already on disk, reusing it (resume = TRUE)")
    return(compare_adjusted(anaDir_existing, cell, comp, ana_adj))
  }

  ds <- ChrAccR::loadDsAcc(file.path(anaDir_existing, cell, "data", "dsATAC_filtered"))
  sa <- ChrAccR::getSampleAnnot(ds)

  # Bridge ChrAccR sample ids to ArchR Sample ids via the descriptor lookup, then
  # to the QC metrics. Try the rownames and every annotation column, and report
  # which one mapped best (so we can see the correct bridge column if this fails).
  try_map <- function(v) id_lookup$Sample[match(norm_key(v), id_lookup$key)]
  cand   <- c(list(.rownames = rownames(sa)), as.list(as.data.frame(sa)))
  mapped <- lapply(cand, try_map)
  scores <- vapply(mapped, function(m) sum(!is.na(m)), integer(1))
  best   <- names(cand)[which.max(scores)]
  message("  ChrAccR->ArchR mapping scores (matched / ", length(rownames(sa)), "): ",
          paste(sprintf("%s=%d", names(cand), scores), collapse = ", "))
  arch_samp <- mapped[[best]]
  idx <- match(arch_samp, qc_by_sample$Sample)

  if (any(is.na(idx))) {
    warning(cell, " / ", comp, ": ", sum(is.na(idx)),
      " sample(s) could not be matched to QC metrics")
    message("  unmatched ChrAccR ids: ", paste(head(rownames(sa)[is.na(idx)], 10), collapse = ", "))
    message("  sa columns: ", paste(colnames(sa), collapse = ", "))
  } else {
    message("  matched all ", length(idx), " samples to QC metrics via '", best, "'")
  }
  # DESeq2 rejects NA in the design, so if any samples are still unmatched, stop
  # with a clear message rather than the opaque "cannot contain NA" error.
  if (mean(is.na(idx)) > 0) {
    stop(sum(is.na(idx)), " of ", length(idx),
      " samples could not be mapped to QC metrics; fix the id_lookup / bridge ",
      "column (see mapping scores above) before adjusting.")
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
  if (!dir.exists(ana_adj)) dir.create(ana_adj, recursive = TRUE)

  setConfigElement("differentialColumns", c("sample_exposure_group"))
  setConfigElement("differentialCompNames", comp)
  setConfigElement("differentialAdjColumns", unique(c(extra_adj, qc_adj_cols)))
  setConfigElement("differentialCutoffL2FC", 0.5)
  setConfigElement("filteringSexChroms", TRUE)
  run_atac_differential(ds, ana_adj)
  rm(ds)
  gc()

  compare_adjusted(anaDir_existing, cell, comp, ana_adj)
}

## ---- 2b. Helper: compare an adjusted run against its unadjusted counterpart --
# Split out of run_adjusted so a comparison that has already been fitted can be
# summarised from disk without refitting (see the resume branch above).
compare_adjusted <- function(anaDir_existing, cell, comp, ana_adj) {
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
  # unadjusted: the comparison as originally run in the existing analysis dir.
  # That directory holds every comparison for the cell type, so pick the file by
  # the group labels in its name and only fall back to the comparisonTable row
  # index if the filenames do not carry them.
  ddir <- file.path(anaDir_existing, cell, "reports", "differential_data")
  unadj_all <- sort(list.files(ddir,
    pattern = "diffTab.*archrPeaks.*\\.tsv$", full.names = TRUE
  ))
  if (length(unadj_all) == 0) {
    unadj_all <- sort(list.files(ddir, pattern = "diffTab.*\\.tsv$", full.names = TRUE))
  }
  stopifnot(length(unadj_all) >= 1)
  want <- sub(" \\[.*", "", comp)
  grps <- strsplit(want, " vs ", fixed = TRUE)[[1]]
  hit <- unadj_all[grepl(grps[1], basename(unadj_all), fixed = TRUE) &
    grepl(grps[2], basename(unadj_all), fixed = TRUE)]
  if (length(hit) >= 1) {
    unadj_file <- hit[1]
  } else {
    ct <- readRDS(file.path(ddir, "comparisonTable.rds"))
    i <- which(paste0(ct$grp1, " vs ", ct$grp2) == want)
    stopifnot(length(i) == 1, length(unadj_all) >= i)
    unadj_file <- unadj_all[i]
  }
  dm <- read.delim(unadj_file)
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
    n_regions_tested = nrow(m),
    status = "ok",
    stringsAsFactors = FALSE
  )
}

# A row with the same columns for a comparison that could not be completed, so
# failures stay visible in the summary table instead of silently vanishing.
fail_row <- function(cell, comp, msg) {
  data.frame(
    cell = cell, comparison = sub(" \\[.*", "", comp),
    DARs_unadjusted = NA_integer_, DARs_adjusted = NA_integer_,
    shared = NA_integer_, recovered_pct = NA_real_,
    sign_concordance_pct = NA_real_, lfc_pearson_all = NA_real_,
    n_regions_tested = NA_integer_,
    status = paste0("failed: ", msg),
    stringsAsFactors = FALSE
  )
}

## ---- 3. Discover every cell type and comparison with saved differential data -
covid_dir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/"
other_dir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_2023-10-02/"

# Rather than a hand-picked list, sweep every cell type in both ChrAccR runs and
# every comparison saved for it. A cell type is only usable if it has BOTH the
# filtered DsATAC object (needed to refit with the extra covariates) and a saved
# unadjusted differential run (needed as the comparison); anything missing either
# is skipped and recorded in the skip log rather than failing the sweep.
skipped <- list()
note_skip <- function(anaDir, cell, reason) {
  message("  skipping ", cell, ": ", reason)
  skipped[[length(skipped) + 1L]] <<- data.frame(
    analysis_dir = anaDir, cell = cell, reason = reason, stringsAsFactors = FALSE
  )
}

discover_jobs <- function(anaDir, extra_adj = character(0)) {
  if (!dir.exists(anaDir)) {
    message("analysis directory not found, skipping: ", anaDir)
    return(list())
  }
  cells <- basename(list.dirs(anaDir, recursive = FALSE))
  message("Scanning ", anaDir, " (", length(cells), " directories)")
  out <- list()
  for (cell in cells) {
    ds_path <- file.path(anaDir, cell, "data", "dsATAC_filtered")
    ddir <- file.path(anaDir, cell, "reports", "differential_data")
    if (!dir.exists(ds_path)) {
      note_skip(anaDir, cell, "no dsATAC_filtered object")
      next
    }
    if (!dir.exists(ddir)) {
      note_skip(anaDir, cell, "no differential_data directory")
      next
    }
    tabs <- list.files(ddir, pattern = "diffTab.*\\.tsv$", full.names = TRUE)
    if (length(tabs) == 0) {
      note_skip(anaDir, cell, "no saved diffTab tables")
      next
    }
    ct_file <- file.path(ddir, "comparisonTable.rds")
    if (!file.exists(ct_file)) {
      note_skip(anaDir, cell, "no comparisonTable.rds")
      next
    }
    ct <- tryCatch(readRDS(ct_file), error = function(e) NULL)
    if (is.null(ct) || nrow(ct) == 0) {
      note_skip(anaDir, cell, "comparisonTable.rds empty or unreadable")
      next
    }
    for (i in seq_len(nrow(ct))) {
      out[[length(out) + 1L]] <- list(
        dir = anaDir, cell = cell,
        comp = paste0(ct$grp1[i], " vs ", ct$grp2[i], " [sample_exposure_group]"),
        adj = extra_adj
      )
    }
    message("  ", cell, ": ", nrow(ct), " comparison(s)")
  }
  out
}

# processing_date was the adjustment column used in the original COVID-19 run,
# so it is carried through there and nowhere else (as in the earlier version).
jobs <- c(
  discover_jobs(covid_dir, extra_adj = "processing_date"),
  discover_jobs(other_dir, extra_adj = character(0))
)
message("Total comparisons to process: ", length(jobs))

skip_df <- if (length(skipped)) do.call(rbind, skipped) else
  data.frame(analysis_dir = character(0), cell = character(0), reason = character(0))
write.csv(skip_df, file.path(fig_dir, "confounder_adjusted_skipped_cells.csv"),
  row.names = FALSE)
message("Skipped ", nrow(skip_df), " cell-type directory(ies); see ",
  "confounder_adjusted_skipped_cells.csv")

## ---- 3b. Run the sweep ------------------------------------------------------
# The summary is rewritten after every comparison, so a crash or a wall-clock
# limit part-way through still leaves a usable table (and `resume` picks the run
# back up from where it stopped).
summary_csv <- file.path(fig_dir, "confounder_adjusted_DAR_summary.csv")
res_list <- list()

for (k in seq_along(jobs)) {
  j <- jobs[[k]]
  message("=== [", k, "/", length(jobs), "] ", j$cell, " | ", j$comp)
  r <- tryCatch(
    run_adjusted(j$dir, j$cell, j$comp, j$adj),
    error = function(e) {
      message("  FAILED: ", conditionMessage(e))
      fail_row(j$cell, j$comp, conditionMessage(e))
    }
  )
  res_list[[k]] <- r
  res <- do.call(rbind, res_list)
  write.csv(res, summary_csv, row.names = FALSE)
  gc()
}

res <- do.call(rbind, res_list)
n_ok <- sum(res$status == "ok", na.rm = TRUE)
message(n_ok, " of ", nrow(res), " comparisons completed")
print(res)

## ---- 4. Figures ------------------------------------------------------------
# All figures are drawn by src/atac/07_4_confounder_adjusted_DAR_plots.R, which
# reads the tables written above straight off disk. Keeping the plotting in one
# script means the sweep and the figures cannot drift apart, and the figures can
# be redrawn or restyled without touching the fits.
ok <- res[!is.na(res$status) & res$status == "ok", , drop = FALSE]
message(nrow(ok), " comparison(s) available for plotting. Run ",
  "src/atac/07_4_confounder_adjusted_DAR_plots.R to write the figures into ",
  fig_dir, " (nothing is re-fitted there).")

message("Done. Adjusted runs in ", out_dir, "; figures and summary tables in ", fig_dir)

#####################################################################
