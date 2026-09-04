#!/usr/bin/env Rscript

#####################################################################
# 01_v2_quality_control.R -- Irem Gunduz, 2026-08-07 (rev. 2026-08-26)
# Association of the IterativeLSI / Harmony dimensions with known covariates.
# Each covariate is tested where it varies: donor level (cohort, age, sex),
# sample level with (1|donor) (sampling day, viral load, per-sample QC), and
# within-C19 for processing batch. Variance decomposition lives in 13b.
#####################################################################

## Load Libraries
suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(ArchR)
  library(dplyr)
  library(ggplot2)
  library(MASS)  # LDA for molecular sex classification
  library(lme4)  # mixed models for the nested / repeated-measures designs
})
set.seed(12)

## ---- Paths / parameters -----------------------------------------------------
# EXPOSURE_ATLAS_REPO / EXPOSURE_ATLAS_ARCHR override the paths below.
repo_candidates <- c(
  Sys.getenv("EXPOSURE_ATLAS_REPO"),
  "/scratch/icbb/igunduz/irem_github/exposure_atlas_manuscript",
  "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript"
)
repo_candidates <- repo_candidates[nzchar(repo_candidates)]
repo_dir <- repo_candidates[which(dir.exists(repo_candidates))[1L]]
if (is.na(repo_dir)) {
  stop("Could not locate the repository. Set EXPOSURE_ATLAS_REPO to its path.")
}
message("repo_dir: ", repo_dir)
suppTables <- file.path(repo_dir, "sample_annots/All_Supplementary_Tables.xlsx")
n_dims     <- 10                                     # top dimensions to test
embeddings_to_test <- c("IterativeLSI", "Harmony")   # pre- and post-correction
outputDir  <- Sys.getenv("EXPOSURE_ATLAS_ARCHR",
  "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/")

# corCutOff = 1 keeps dimensions ArchR's 0.75 default would drop for depth.
archr_cor_cutoff <- 0.75   # the default that 02_cluster_and_batch.R relies on

# Load the ArchR project (already subsetted to non-BA samples).
# loadArchRProject() remaps the stored GroupCoverages paths with a plain gsub of
# the old output directory and then asserts all of them exist:
#     zfiles <- gsub(outputDir, outputDirNew, zdata$File)
#     stopifnot(all(file.exists(zfiles)))
# That assertion is NOT covered by force = TRUE, so a moved project or deleted
# coverage files abort the load. Nothing in this script touches coverages, so
# fall back to a loader that repairs what it can and warns about the rest.
load_archr_lenient <- function(path) {
  p <- normalizePath(path, mustWork = TRUE)
  proj <- ArchR::recoverArchRProject(readRDS(file.path(p, "Save-ArchR-Project.rds")))

  arrows <- file.path(p, "ArrowFiles", basename(proj@sampleColData$ArrowFiles))
  if (!all(file.exists(arrows))) {
    stop("Missing arrow file(s) under ", file.path(p, "ArrowFiles"), ": ",
         paste(basename(arrows[!file.exists(arrows)]), collapse = ", "))
  }
  proj@sampleColData$ArrowFiles <- arrows

  gcov <- proj@projectMetadata$GroupCoverages
  if (length(gcov)) {
    for (z in seq_along(gcov)) {
      f   <- gcov[[z]]$coverageMetadata$File
      new <- file.path(p, "GroupCoverages", basename(dirname(f)), basename(f))
      hit <- file.exists(new)
      gcov[[z]]$coverageMetadata$File <- ifelse(hit, new, f)
      if (!all(hit)) {
        nm <- if (!is.null(names(gcov))) names(gcov)[z] else as.character(z)
        message("GroupCoverages '", nm, "': ", sum(!hit), "/", length(hit),
                " coverage file(s) not found. Not needed here, but 05_pseudobulk.R ",
                "would have to regenerate them.")
      }
    }
    proj@projectMetadata$GroupCoverages <- gcov
  }
  proj@projectMetadata$outputDirectory <- p
  proj   # nothing on disk is modified
}
project <- tryCatch(
  ArchR::loadArchRProject(outputDir, showLogo = FALSE),
  error = function(e) {
    message("loadArchRProject() failed: ", conditionMessage(e))
    message("Retrying without the GroupCoverages assertion.")
    load_archr_lenient(outputDir)
  }
)

## ---- Sample-level covariate table (Table S1) --------------------------------
# Join on arrow_name ("ATAC_055_fragments.tsv.gz"), NOT sampleId ("C19_mod_055").
covar <- readxl::read_excel(suppTables, sheet = "Table S1") %>%
  dplyr::select(sampleId, arrow_name, exposure_type, exposure_group,
    exposure_grouping, record_id, Age, Sex,
    supplier, processing_institution, tss_cutoff, log10frags_cutoff) %>%
  dplyr::mutate(
    # Cohort: C19 / HIV / Influenza / OP
    exposure_type = factor(exposure_type),
    # within-cohort condition arm, tested separately from cohort
    exposure_group = factor(exposure_group),
    # OP has no unexposed arm, so Control_status is partly nested in cohort.
    control_status = factor(
      ifelse(exposure_grouping == "healthy", "Control", "Exposed"),
      levels = c("Control", "Exposed")
    ),
    Sex = factor(Sex),
    Age = as.numeric(Age)
  )

# Age/sex recovered for C19; sampling day is relative to onset per cohort.

# ATAC label -> clinical Donor ID (C19). 555/564/66 are repeat draws.
atac_to_donor <- c(
  ATAC_055 = "55650-0055", ATAC_057 = "55650-0057", ATAC_132D0 = "55650-0132d0",
  ATAC_52 = "55650-0052", ATAC_555_1 = "28205-0555d0", ATAC_555_2 = "28205-0555d2",
  ATAC_556 = "28205-0556", ATAC_557 = "28205-0557", ATAC_558 = "28205-0558",
  ATAC_559 = "28205-0559", ATAC_560 = "28205-0560", ATAC_564A = "28205-0564d0",
  ATAC_564B = "28205-0564d2", ATAC_66D0 = "55650-0066d0", ATAC_66D7 = "55650-0066d7",
  ATAC_67 = "55650-0067", ATAC_83 = "55650-0083", ATAC_86 = "55650-0086",
  ATAC_EV08 = "EV08", ATAC_HIP02_frozen = "HIP002", ATAC_HIP023_frozen = "HIP023",
  ATAC_HIP043 = "HIP043", ATAC_HIP044 = "HIP044", ATAC_HIP045 = "HIP045",
  ATAC_HIP15_frozen = "HIP015"
)

strip_timepoint <- function(x) sub("d[0-9]+$", "", x)

# COVID-19 clinical metadata (one row per Donor; place the CSV in sample_annots/)
pm <- read.csv(file.path(repo_dir, "sample_annots/patient_metadata.csv"),
  stringsAsFactors = FALSE
)

# HIV longitudinal metadata: day relative to seroconversion (SC) + donor info
hiv_meta <- data.frame(
  stem = c("hiv6", "hiv12", "hiv9", "hiv8", "hiv4", "hiv1",
           "hiv2", "hiv7", "hiv3", "hiv11", "hiv10", "hiv5"),
  donor = c(rep("9313454", 3), rep("9313664", 3), rep("9320086", 3), rep("9320666", 3)),
  day_to_SC = c(-294, 0, 189, -232, 0, 265, -114, 0, 196, -364, 0, 253),
  viral_load = c(0, 1866453, 794, 0, 852580, 110691, 40, 415480, 244551, 40, 27250, 12341),
  age = c(rep(23, 3), rep(22, 3), rep(22, 3), rep(19, 3)),
  sex = "Male", stringsAsFactors = FALSE
)

# bin midpoint, so binned C19 ages can be used as a continuous covariate
bin_to_numeric <- function(x) {
  vapply(x, function(s) {
    if (is.na(s) || s == "" || s == "NA") return(NA_real_)
    if (grepl("\\+$", s)) return(as.numeric(sub("\\+", "", s)) + 5)
    p <- as.numeric(strsplit(s, "-")[[1]])
    if (length(p) == 2) mean(p) else suppressWarnings(as.numeric(s))
  }, numeric(1), USE.NAMES = FALSE)
}

label      <- sub("_fragments\\.tsv\\.gz$", "", covar$arrow_name)
sample_key <- unname(atac_to_donor[label])   # donor + timepoint
pm_i       <- match(sample_key, pm$Donor)
hiv_i      <- match(label, hiv_meta$stem)
# "day 30 after challenge" is a mislabel: the final timepoint is day 28.
flu_day <- c(
  "right before challenge" = -1, "day 3 after challenge" = 3,
  "day 6 after challenge" = 6, "day 30 after challenge" = 28,
  "day 28 after challenge" = 28
)

covar <- covar %>% dplyr::mutate(
  donor_id = strip_timepoint(dplyr::coalesce(sample_key, hiv_meta$donor[hiv_i])),
  # reported age: exact where known, bin label for COVID-19
  Age_reported = dplyr::case_when(
    !is.na(Age) ~ as.character(Age),
    !is.na(pm_i) ~ pm$Age[pm_i],
    !is.na(hiv_i) ~ as.character(hiv_meta$age[hiv_i]),
    TRUE ~ NA_character_
  ),
  # numeric age for the association test (bin midpoints for COVID-19)
  Age_numeric = dplyr::coalesce(
    Age, bin_to_numeric(ifelse(is.na(pm_i), NA_character_, pm$Age[pm_i])),
    as.numeric(hiv_meta$age[hiv_i])
  ),
  Sex = factor(dplyr::coalesce(
    as.character(Sex),
    c(M = "Male", F = "Female")[pm$Sex[pm_i]],
    hiv_meta$sex[hiv_i]
  ), levels = c("Female", "Male")),
  # sampling time relative to infection onset
  sampling_day = dplyr::case_when(
    exposure_type == "C19" & !is.na(pm_i) & pm$days_post_pos_test[pm_i] != "NA" ~
      suppressWarnings(as.numeric(pm$days_post_pos_test[pm_i])),
    exposure_type == "HIV" & !is.na(hiv_i) ~ as.numeric(hiv_meta$day_to_SC[hiv_i]),
    exposure_type == "Influenza" ~ unname(flu_day[as.character(record_id)]),
    TRUE ~ NA_real_
  ),
  onset_reference = dplyr::case_when(
    exposure_type == "C19" & !is.na(sampling_day) ~ "days after first positive SARS-CoV-2 test",
    exposure_type == "HIV" ~ "days relative to seroconversion",
    exposure_type == "Influenza" ~ "days after influenza challenge",
    exposure_type == "OP" ~ "chronic exposure; no datable onset",
    TRUE ~ NA_character_
  ),
  viral_load = ifelse(is.na(hiv_i), NA_real_, hiv_meta$viral_load[hiv_i])
)

# unresolved donor -> fall back to the sample, not a shared NA group
# influenza: ATAC_flu<ID>_<timepoint> -> ATAC_flu<ID>, 4 timepoints per subject
is_flu <- grepl("^ATAC_flu[0-9]+_[0-9]+$", label)
covar$donor_id[is.na(covar$donor_id) & is_flu] <-
  sub("^(ATAC_flu[0-9]+)_[0-9]+$", "\\1", label[is.na(covar$donor_id) & is_flu])
covar$donor_id[is.na(covar$donor_id)] <- covar$arrow_name[is.na(covar$donor_id)]

message("Recorded metadata coverage after recovery (n samples with a value):")
print(covar %>%
  dplyr::group_by(exposure_type) %>%
  dplyr::summarise(
    n = dplyr::n(), n_donors = dplyr::n_distinct(donor_id),
    sex = sum(!is.na(Sex)), age = sum(!is.na(Age_numeric)),
    sampling_day = sum(!is.na(sampling_day)), .groups = "drop"
  ))

rep_don <- covar %>% dplyr::count(donor_id) %>% dplyr::filter(n > 1)
message(sprintf(
  "Repeated sampling: %d sample(s) from %d donor(s) contribute more than one sample (max %d). The donor-level tests below account for this.",
  sum(rep_don$n), nrow(rep_don), if (nrow(rep_don)) max(rep_don$n) else 0L
))

# BATCH. processing_date is C19-only (710/720), so it is a within-C19 covariate.
# supplier is a relabelling of cohort and processing_institution is constant, so
# neither adds anything; tss_cutoff / log10frags_cutoff do vary within cohort.

atac_meta_covid_path <- file.path(repo_dir, "sample_annots/ATAC_metadata_covid.csv")
if (file.exists(atac_meta_covid_path)) {
  amc <- read.csv(atac_meta_covid_path, stringsAsFactors = FALSE)
  # Join on arrow_name: the record_id column here is truncated for some samples.
  stopifnot(!anyDuplicated(amc$arrow_name))
  bi <- match(label, amc$arrow_name)
  covar$processing_batch <- factor(ifelse(is.na(bi), NA_character_,
                                          as.character(amc$processing_date[bi])))
  message(sprintf(
    "Processing batch recovered for %d sample(s) (cohorts: %s); levels: %s",
    sum(!is.na(covar$processing_batch)),
    paste(sort(unique(as.character(covar$exposure_type[!is.na(covar$processing_batch)]))),
          collapse = ", "),
    paste(levels(covar$processing_batch), collapse = ", ")
  ))
} else {
  covar$processing_batch <- factor(NA_character_)
  warning("ATAC_metadata_covid.csv not found -- processing batch will not be tested.")
}

covar <- covar %>% dplyr::mutate(
  supplier = factor(supplier),
  processing_institution = factor(processing_institution),
  tss_cutoff = as.numeric(tss_cutoff),
  log10frags_cutoff = as.numeric(log10frags_cutoff)
)

## ---- Confound structure among the design and technical variables ------------
confound_tab <- function(a, b, name_a, name_b) {
  keep <- !is.na(a) & !is.na(b)
  if (!any(keep)) return(NULL)
  tb <- table(droplevels(factor(a[keep])), droplevels(factor(b[keep])))
  d <- as.data.frame(tb, stringsAsFactors = FALSE)
  names(d) <- c(name_a, name_b, "n")
  d[d$n > 0, , drop = FALSE]
}

message("\n--- Cohort vs supplier (expected to be aliased) ---")
print(with(covar, table(cohort = exposure_type, supplier = supplier)))
message("processing_institution levels: ",
  paste(levels(droplevels(covar$processing_institution)), collapse = ", "),
  " (constant -> not testable)")

# Is the C19 batch balanced across severity arms? If not, say so.
c19 <- covar[covar$exposure_type == "C19" & !is.na(covar$processing_batch), ]
if (nrow(c19) > 0) {
  bt <- table(group = droplevels(c19$exposure_group),
              batch = droplevels(c19$processing_batch))
  message("\n--- COVID-19 processing batch vs severity arm ---")
  print(bt)
  ft <- tryCatch(stats::fisher.test(bt, simulate.p.value = TRUE, B = 1e5),
                 error = function(e) NULL)
  if (!is.null(ft)) {
    message(sprintf(
      "Fisher exact test of batch vs severity arm within C19: p = %.4g%s",
      ft$p.value,
      if (ft$p.value < 0.05)
        "  <-- BATCH IS NOT BALANCED ACROSS THE C19 ARMS; any within-C19 contrast is partly confounded with it"
      else "  (balanced)"
    ))
  }
  # donor-level version, since repeat draws of one donor share a batch
  bt_don <- table(
    group = droplevels(c19$exposure_group[!duplicated(c19$donor_id)]),
    batch = droplevels(c19$processing_batch[!duplicated(c19$donor_id)])
  )
  message("--- same table at the donor level ---")
  print(bt_don)
}

confound_rows <- dplyr::bind_rows(
  confound_tab(covar$exposure_type, covar$supplier, "cohort", "supplier"),
  confound_tab(covar$exposure_group, covar$processing_batch, "exposure_group", "processing_batch"),
  confound_tab(covar$exposure_type, covar$tss_cutoff, "cohort", "tss_cutoff"),
  confound_tab(covar$exposure_type, covar$log10frags_cutoff, "cohort", "log10frags_cutoff")
)
write.csv(confound_rows,
  file = file.path(repo_dir, "figures/covariate_confound_structure.csv"),
  row.names = FALSE)


## ---- Per-sample QC metrics (descriptive; the real test is Part 3) -----------
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
  )
covar <- dplyr::left_join(covar, qc_by_sample, by = c("arrow_name" = "Sample"))

## ---- Cells behind each pseudobulk (sample x cell type) ----------------------
# Raw counts, including combinations below the >50-cell cutoff in 05_pseudobulk.R.
ct_col <- "ClusterCellTypes" # cell-type annotation used across the atlas

pb_ct <- as.data.frame(getCellColData(project, select = c("Sample", ct_col)))
names(pb_ct)[names(pb_ct) == ct_col] <- "CellType"
n_unannotated <- sum(is.na(pb_ct$CellType))
pb_ct <- pb_ct[!is.na(pb_ct$CellType), , drop = FALSE]

# one row per sample, one column per cell type (zeros where a sample has none)
pb_wide <- as.data.frame.matrix(table(pb_ct$Sample, pb_ct$CellType))
pb_wide <- pb_wide[, sort(colnames(pb_wide)), drop = FALSE]
pb_cellnames <- colnames(pb_wide)
colnames(pb_wide) <- paste0("n_cells_", gsub("[^A-Za-z0-9]+", "_", pb_cellnames))
pb_wide$arrow_name <- rownames(pb_wide)
rownames(pb_wide) <- NULL
pb_count_cols <- setdiff(colnames(pb_wide), "arrow_name")

# a duplicated key would expand covar and break the Table S1 append below
stopifnot(!anyDuplicated(pb_wide$arrow_name), !anyDuplicated(qc_by_sample$Sample))
n_covar_before <- nrow(covar)
covar <- dplyr::left_join(covar, pb_wide, by = "arrow_name")
stopifnot(nrow(covar) == n_covar_before)

# per-cell-type counts should sum to n_cells, minus unannotated cells
pb_mat <- as.matrix(as.data.frame(covar)[, pb_count_cols, drop = FALSE])
pb_rowsum <- rowSums(pb_mat, na.rm = TRUE)
n_match <- sum(pb_rowsum == covar$n_cells, na.rm = TRUE)
message(sprintf(
  "Pseudobulk cell counts: %d cell type(s), %d/%d samples where the per-cell-type counts sum to n_cells (%d cell(s) carry no cell-type label)",
  length(pb_count_cols), n_match, nrow(covar), n_unannotated
))
if (n_match < sum(!is.na(covar$n_cells))) {
  message("  samples where the sums differ (per-cell-type total vs n_cells):")
  d <- which(pb_rowsum != covar$n_cells)
  print(data.frame(
    arrow_name = covar$arrow_name[d], per_celltype_total = pb_rowsum[d],
    n_cells = covar$n_cells[d]
  ))
}

# how many pseudobulks clear the >50-cell threshold 05_pseudobulk.R uses
message(sprintf(
  "Pseudobulks with > 50 cells: %d of %d sample x cell-type combinations",
  sum(pb_mat > 50, na.rm = TRUE),
  nrow(pb_mat) * length(pb_count_cols)
))


## ---- Predict donor sex from chrY and XIST accessibility ---------------------
# Supervised: threshold XIST -> Sex_by_xist, train LDA on recorded sex where it
# agrees, predict all. chrY is held out of the training-set definition, so
# agreement with Sex_by_chrY is the only non-circular check.
xist_gr <- GenomicRanges::GRanges(
  "chrX", IRanges::IRanges(start = 73820651, end = 73852753)
) # XIST gene, hg38

# Reading chrX/chrY fragments from every arrow is the slow step, so it is cached.
# Only the raw counts are cached; the thresholds and the LDA below are recomputed
# every run, so changing the classifier does not need a re-extraction.
# Force a re-extraction with REDO_SEX=1.
sex_counts_csv <- file.path(repo_dir, "figures/sex_fragment_counts.csv")
count_cols <- c("Sample", "chrY", "chrX", "xist",
                "chrY_frac", "xist_frac", "chrY_chrX_ratio")

extract_sex_counts <- function() {
  arrowFiles      <- ArchR::getArrowFiles(project)
  cell_sample_all <- getCellColData(project, select = "Sample", drop = TRUE)
  all_cellNames   <- project$cellNames
  nfrags_all      <- getCellColData(project, select = "nFrags", drop = TRUE)
  total_by_sample <- tapply(nfrags_all, cell_sample_all, sum)

  out <- do.call(rbind, lapply(arrowFiles, function(af) {
    samp  <- gsub("\\.arrow$", "", basename(af))
    cells <- all_cellNames[cell_sample_all == samp]
    if (length(cells) == 0) {
      return(NULL)
    }
    # chrY may be absent from the arrow; fall back to zero rather than aborting
    fy <- tryCatch(
      ArchR::getFragmentsFromArrow(af, chr = "chrY", cellNames = cells, verbose = FALSE),
      error = function(e) {
        warning("no chrY fragments retrievable for ", samp, ": ", conditionMessage(e))
        GenomicRanges::GRanges()
      }
    )
    fx <- ArchR::getFragmentsFromArrow(af, chr = "chrX", cellNames = cells, verbose = FALSE)
    n_xist <- sum(IRanges::overlapsAny(fx, xist_gr))
    total  <- as.numeric(total_by_sample[samp])
    data.frame(
      Sample = samp,
      chrY = length(fy), chrX = length(fx), xist = n_xist,
      chrY_frac = length(fy) / total,
      xist_frac = n_xist / total,
      chrY_chrX_ratio = length(fy) / length(fx),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  out
}

project_samples <- unique(as.character(getCellColData(project, select = "Sample", drop = TRUE)))
sex_metrics <- NULL
if (file.exists(sex_counts_csv) && !nzchar(Sys.getenv("REDO_SEX"))) {
  cached <- read.csv(sex_counts_csv, stringsAsFactors = FALSE)
  missing_cols    <- setdiff(count_cols, names(cached))
  missing_samples <- setdiff(project_samples, cached$Sample)
  if (length(missing_cols)) {
    message("Sex-count cache is missing column(s): ",
      paste(missing_cols, collapse = ", "), " -- re-extracting.")
  } else if (length(missing_samples)) {
    message("Sex-count cache is missing ", length(missing_samples),
      " sample(s) present in the project -- re-extracting.")
  } else {
    sex_metrics <- cached[cached$Sample %in% project_samples, count_cols, drop = FALSE]
    rownames(sex_metrics) <- NULL
    message("Loaded cached chrX/chrY/XIST counts for ", nrow(sex_metrics),
      " sample(s) from ", sex_counts_csv, " (REDO_SEX=1 to re-extract).")
  }
}
if (is.null(sex_metrics)) {
  sex_metrics <- extract_sex_counts()
  write.csv(sex_metrics, sex_counts_csv, row.names = FALSE)
  message("Extracted and cached chrX/chrY/XIST counts to ", sex_counts_csv)
}

# attach recorded sex + cohort for calibration and diagnostics
sex_metrics <- dplyr::left_join(sex_metrics,
  dplyr::select(covar, arrow_name, exposure_type, Sex),
  by = c("Sample" = "arrow_name")
)

# data-driven pseudocount (half the smallest non-zero value), used throughout
half_min_nonzero <- function(x) min(x[x > 0], na.rm = TRUE) / 2
eps_chrY <- half_min_nonzero(sex_metrics$chrY_frac)
eps_xist <- half_min_nonzero(sex_metrics$xist_frac)
sex_metrics$log_chrY <- log10(sex_metrics$chrY_frac + eps_chrY)
sex_metrics$log_xist <- log10(sex_metrics$xist_frac + eps_xist)

# Do NOT cluster on chrY and XIST jointly: chrY carries cohort-specific offsets,
# so the split follows cohort, not sex. Threshold on XIST, cross-check with chrY.
gap_threshold <- function(x, eps, lo = 0.10, hi = 0.90) {
  lx   <- sort(log10(x + eps))
  n    <- length(lx)
  i_lo <- max(1L, floor(lo * n))
  i_hi <- min(n - 1L, ceiling(hi * n))
  idx  <- i_lo:i_hi
  i    <- idx[which.max(diff(lx)[idx])]
  10^mean(c(lx[i], lx[i + 1L])) - eps  # geometric midpoint of the largest gap
}

# Marker-based reference calls: high XIST = female, high chrY = male.
xist_thr <- gap_threshold(sex_metrics$xist_frac, eps_xist)
sex_metrics$Sex_by_xist <- factor(
  ifelse(sex_metrics$xist_frac >= xist_thr, "Female", "Male"),
  levels = c("Female", "Male")
)
chrY_thr <- gap_threshold(sex_metrics$chrY_frac, eps_chrY)
sex_metrics$Sex_by_chrY <- factor(
  ifelse(sex_metrics$chrY_frac >= chrY_thr, "Male", "Female"),
  levels = c("Female", "Male")
)
sex_metrics$features_agree <- sex_metrics$Sex_by_xist == sex_metrics$Sex_by_chrY
message(sprintf(
  "Marker thresholds: XIST %.3g, chrY %.3g; the two markers agree for %d/%d samples",
  xist_thr, chrY_thr, sum(sex_metrics$features_agree), nrow(sex_metrics)
))

## ---- LINEAR CLASSIFIER (LDA) on the two log features ------------------------
# Training on all recorded labels inverted the chrY -> male direction.
lab       <- !is.na(sex_metrics$Sex)
train_idx <- lab & sex_metrics$Sex == sex_metrics$Sex_by_xist
message(sprintf(
  "LDA training set: %d high-confidence samples (%d recorded labels excluded as inconsistent with XIST)",
  sum(train_idx), sum(lab) - sum(train_idx)
))
stopifnot(sum(train_idx) >= 6, dplyr::n_distinct(sex_metrics$Sex[train_idx]) == 2)

train_df <- data.frame(
  Sex = droplevels(factor(sex_metrics$Sex[train_idx], levels = c("Female", "Male"))),
  log_chrY = sex_metrics$log_chrY[train_idx],
  log_xist = sex_metrics$log_xist[train_idx]
)

# equal priors: the training subset is ~87% male and would bias calls otherwise
equal_prior <- c(Female = 0.5, Male = 0.5)

# LOO within a concordance-selected set: separability, not generalisation error
loo <- MASS::lda(Sex ~ log_chrY + log_xist,
  data = train_df, prior = equal_prior, CV = TRUE
)
loo_correct  <- sum(loo$class == train_df$Sex)
loo_n        <- nrow(train_df)
loo_accuracy <- loo_correct / loo_n
# Wilson interval, so the reported accuracy carries its uncertainty
loo_ci <- suppressWarnings(
  stats::prop.test(loo_correct, loo_n, correct = FALSE)$conf.int)
message(sprintf(
  "LDA leave-one-out accuracy on the (concordance-selected) training set: %.1f%% (%d/%d, 95%% CI %.1f-%.1f%%)",
  100 * loo_accuracy, loo_correct, loo_n, 100 * loo_ci[1], 100 * loo_ci[2]
))

# fit on the full training set and predict all samples
lda_fit <- MASS::lda(Sex ~ log_chrY + log_xist, data = train_df, prior = equal_prior)
print(lda_fit)
pred <- stats::predict(lda_fit, newdata = sex_metrics[, c("log_chrY", "log_xist")])
sex_metrics$Sex_predicted <- factor(as.character(pred$class), levels = c("Female", "Male"))
sex_metrics$p_male <- round(pred$posterior[, "Male"], 4)

# chrY was held out of the training-set definition; XIST agreement is circular.
agree_chrY <- mean(sex_metrics$Sex_predicted == sex_metrics$Sex_by_chrY)
agree_xist <- mean(sex_metrics$Sex_predicted == sex_metrics$Sex_by_xist)
message(sprintf(
  "LDA call vs chrY threshold (HELD OUT): %.1f%% agreement across all %d samples",
  100 * agree_chrY, nrow(sex_metrics)
))
message(sprintf(
  "LDA call vs XIST threshold (NOT independent -- used to select training set): %.1f%%",
  100 * agree_xist
))

message("Predicted sex by cohort:")
print(table(cohort = sex_metrics$exposure_type, predicted = sex_metrics$Sex_predicted))
low_conf <- sex_metrics$p_male > 0.05 & sex_metrics$p_male < 0.95
if (any(low_conf)) {
  message("Low-confidence calls (0.05 < P(male) < 0.95):")
  print(sex_metrics[low_conf,
    c("Sample", "exposure_type", "chrY_frac", "xist_frac", "p_male", "Sex_predicted")
  ])
}

# recorded sex used only to flag metadata errors, not as validation
val <- sex_metrics[lab, ]
message(sprintf(
  "Molecular sex vs recorded sex: %.1f%% concordant (%d labelled samples; not an independent check)",
  100 * mean(val$Sex_predicted == val$Sex), nrow(val)
))
print(table(recorded = val$Sex, predicted = val$Sex_predicted))
sex_metrics$sex_flag <- ifelse(
  lab & sex_metrics$Sex_predicted != sex_metrics$Sex, "check_metadata", ""
)
disc <- sex_metrics[sex_metrics$sex_flag == "check_metadata",
  c("Sample", "exposure_type", "Sex", "Sex_predicted", "chrY_frac", "xist_frac")
]
message("Recorded sex disagrees with the molecular call for ", nrow(disc),
  " sample(s) (verify metadata):"
)
print(disc)

## ---- Classifier performance, saved for the manuscript text ------------------
sex_lda_metrics <- data.frame(
  metric = c(
    "n_samples_total",
    "n_recorded_sex",
    "n_training",
    "n_excluded_from_training",
    "loo_correct",
    "loo_accuracy",
    "loo_accuracy_ci_low",
    "loo_accuracy_ci_high",
    "agreement_with_chrY_threshold_heldout",
    "agreement_with_xist_threshold_circular",
    "concordance_with_recorded_sex",
    "n_discordant_with_recorded_sex",
    "n_low_confidence_calls",
    "xist_threshold",
    "chrY_threshold",
    "prior_female",
    "prior_male"
  ),
  value = c(
    nrow(sex_metrics),
    sum(lab),
    loo_n,
    sum(lab) - loo_n,
    loo_correct,
    round(loo_accuracy, 4),
    round(loo_ci[1], 4),
    round(loo_ci[2], 4),
    round(agree_chrY, 4),
    round(agree_xist, 4),
    round(mean(val$Sex_predicted == val$Sex), 4),
    nrow(disc),
    sum(low_conf),
    signif(xist_thr, 4),
    signif(chrY_thr, 4),
    equal_prior[["Female"]],
    equal_prior[["Male"]]
  ),
  stringsAsFactors = FALSE
)
write.csv(sex_lda_metrics,
  file = file.path(repo_dir, "figures/sex_lda_classifier_metrics.csv"),
  row.names = FALSE)
message("LDA classifier metrics written to figures/sex_lda_classifier_metrics.csv:")
print(sex_lda_metrics)

# A ready-to-paste sentence, so the methods text and the file cannot drift apart
message(sprintf(
  "Methods sentence: LDA classifier trained on the %d samples whose recorded sex agreed with the XIST-based call (equal priors; leave-one-out accuracy %.1f%%, %d/%d); the resulting calls agreed with the held-out chrY threshold for %.1f%% of all %d samples.",
  loo_n, 100 * loo_accuracy, loo_correct, loo_n, 100 * agree_chrY, nrow(sex_metrics)
))

## ---- Diagnostic figure: molecular sex classification ------------------------
sex_grid <- expand.grid(
  log_chrY = seq(min(sex_metrics$log_chrY), max(sex_metrics$log_chrY), length.out = 200),
  log_xist = seq(min(sex_metrics$log_xist), max(sex_metrics$log_xist), length.out = 200)
)
sex_grid$p_male <- stats::predict(lda_fit, newdata = sex_grid)$posterior[, "Male"]

sex_metrics$Recorded <- ifelse(is.na(sex_metrics$Sex), "not recorded",
  as.character(sex_metrics$Sex)
)
sex_metrics$Discordant <- sex_metrics$sex_flag == "check_metadata"

base_sex_plot <- function(dat) {
  ggplot(dat, aes(x = log_chrY, y = log_xist)) +
    geom_contour(
      data = sex_grid, aes(z = p_male), breaks = 0.5,
      colour = "grey35", linetype = "dashed", linewidth = 0.5
    ) +
    geom_point(aes(colour = Sex_predicted, shape = Recorded), size = 2.8, alpha = 0.9) +
    scale_colour_manual(values = c(Female = "#E64B35", Male = "#3C5488"), name = "Predicted") +
    scale_shape_manual(values = c(Female = 17, Male = 15, `not recorded` = 1), name = "Recorded") +
    labs(
      x = expression(log[10] ~ "chrY fragment fraction"),
      y = expression(log[10] ~ "XIST accessibility fraction")
    ) +
    theme_classic(base_size = 12)
}

p_sex <- base_sex_plot(sex_metrics) +
  geom_point(
    data = subset(sex_metrics, Discordant), shape = 21, size = 5,
    stroke = 0.8, colour = "black", fill = NA
  ) +
  labs(
    title = "Molecular sex inference from chrY and XIST accessibility",
    subtitle = paste0(
      "Dashed line, LDA boundary. Black rings, recorded sex disagrees with the\n",
      "molecular call (n = ", sum(sex_metrics$Discordant), ")."
    )
  ) +
  theme(plot.subtitle = element_text(size = 9))
ggsave(file.path(repo_dir, "figures/sex_prediction_scatter.pdf"), p_sex,
  width = 8.5, height = 5.5
)

# same view split by cohort: shows the separation is cohort-independent
p_sex_cohort <- base_sex_plot(sex_metrics) +
  facet_wrap(~exposure_type) +
  labs(title = "Molecular sex inference by cohort")
ggsave(file.path(repo_dir, "figures/sex_prediction_by_cohort.pdf"), p_sex_cohort,
  width = 9, height = 6.5
)
message("Wrote sex-prediction figures to ", file.path(repo_dir, "figures/"))

# per-sample evidence behind the sex call; kept out of Table S1
write.csv(
  sex_metrics[, c(
    "Sample", "exposure_type", "Sex", "Sex_predicted", "p_male",
    "Sex_by_xist", "Sex_by_chrY", "features_agree",
    "chrY_frac", "xist_frac", "chrY_chrX_ratio", "sex_flag"
  )],
  file = file.path(repo_dir, "figures/sex_prediction_metrics.csv"), row.names = FALSE
)

# carry the raw chrX / chrY / XIST signals so the call can be checked
stopifnot(!anyDuplicated(sex_metrics$Sample))
n_covar_before <- nrow(covar)
covar <- dplyr::left_join(covar,
  dplyr::select(
    sex_metrics, Sample, chrY, chrX, xist, chrY_frac, xist_frac,
    chrY_chrX_ratio, p_male,
    Sex_predicted, Sex_by_xist, Sex_by_chrY, sex_flag
  ),
  by = c("arrow_name" = "Sample")
)
stopifnot(nrow(covar) == n_covar_before)


#####################################################################
# ASSOCIATION TESTING MACHINERY
#####################################################################

## ---- Is this covariate constant within every group? -------------------------
is_group_constant <- function(x, group) {
  ok <- !is.na(x) & !is.na(group)
  if (!any(ok)) return(TRUE)
  v <- as.character(x[ok]); g <- as.character(group[ok])
  all(tapply(v, g, function(z) length(unique(z))) == 1L)
}

## ---- Collapse y / x / cohort to one row per group ---------------------------
collapse_to_group <- function(y, x, group, cohort = NULL) {
  g  <- droplevels(factor(group))
  ym <- tapply(y, g, mean)
  if (is.factor(x)) {
    xv <- tapply(as.character(x), g, function(z) z[1L])
    xm <- droplevels(factor(unname(xv), levels = levels(x)))
  } else {
    xm <- as.numeric(tapply(as.numeric(x), g, mean))
  }
  cm <- NULL
  if (!is.null(cohort)) {
    cv <- tapply(as.character(cohort), g, function(z) z[1L])
    cm <- droplevels(factor(unname(cv), levels = levels(cohort)))
  }
  list(y = as.numeric(ym), x = xm, cohort = cm, n = nlevels(g))
}

## ---- One association test, at the right level -------------------------------
# group = donor (sample-level) or sample (pseudobulk-level). With `cohort`, the
# test is partial. Returns raw + adjusted R^2, p, n, and the unit tested.
assoc_test <- function(y, x, group, cohort = NULL,
                       collapsed_unit = "per donor",
                       obs_unit = "per sample, donor-adjusted") {
  keep <- !is.na(y) & !is.na(x) & !is.na(group)
  if (!is.null(cohort)) keep <- keep & !is.na(cohort)
  fail <- list(r2 = NA_real_, adj_r2 = NA_real_, p = NA_real_,
               n = sum(keep), n_group = NA_integer_, unit = NA_character_)
  if (sum(keep) < 4L) return(fail)

  y  <- y[keep]
  g  <- droplevels(factor(as.character(group[keep])))
  x  <- if (is.factor(x)) droplevels(x[keep]) else as.numeric(x[keep])
  ch <- if (is.null(cohort)) NULL else droplevels(cohort[keep])

  usable <- function(v) if (is.factor(v)) nlevels(v) >= 2L else stats::sd(v) > 0
  if (!usable(x)) return(fail)

  if (is_group_constant(x, g)) {
    ## ---- group (donor) level ordinary least squares -------------------------
    cl <- collapse_to_group(y, x, g, ch)
    if (!usable(cl$x) || cl$n < 4L) return(fail)
    dd <- data.frame(y = cl$y, x = cl$x)
    if (!is.null(cl$cohort) && nlevels(cl$cohort) >= 2L) {
      dd$ch <- cl$cohort
      f0 <- y ~ ch
    } else {
      f0 <- y ~ 1
    }
    f1 <- stats::update(f0, . ~ . + x)
    fit0 <- try(stats::lm(f0, data = dd), silent = TRUE)
    fit1 <- try(stats::lm(f1, data = dd), silent = TRUE)
    if (inherits(fit0, "try-error") || inherits(fit1, "try-error")) return(fail)
    an <- stats::anova(fit0, fit1)
    if (is.na(an[["Df"]][2L]) || an[["Df"]][2L] <= 0) return(fail)
    rss0 <- sum(stats::resid(fit0)^2); df0 <- an[["Res.Df"]][1L]
    rss1 <- sum(stats::resid(fit1)^2); df1 <- an[["Res.Df"]][2L]
    if (rss0 <= 0 || df1 <= 0) return(fail)
    list(
      r2      = (rss0 - rss1) / rss0,
      # adjusted: a many-level factor scores a large raw R^2 on ~40 donors
      adj_r2  = 1 - (rss1 / df1) / (rss0 / df0),
      p       = unname(an[["Pr(>F)"]][2L]),
      n       = cl$n, n_group = cl$n, unit = collapsed_unit
    )
  } else {
    ## ---- observation level, group random intercept --------------------------
    if (nlevels(g) < 3L) return(fail)
    dd <- data.frame(y = y, x = x, g = g)
    if (!is.null(ch) && nlevels(ch) >= 2L) {
      dd$ch <- ch
      f0 <- y ~ ch + (1 | g); f1 <- y ~ ch + x + (1 | g)
    } else {
      f0 <- y ~ 1 + (1 | g);  f1 <- y ~ x + (1 | g)
    }
    # singular fits expected: most donors have one sample, so it degrades to OLS
    fits <- try(suppressMessages(suppressWarnings({
      ctrl <- lme4::lmerControl(check.conv.singular = "ignore")
      m0 <- lme4::lmer(f0, data = dd, REML = FALSE, control = ctrl)
      m1 <- lme4::lmer(f1, data = dd, REML = FALSE, control = ctrl)
      list(m0 = m0, m1 = m1, an = stats::anova(m0, m1))
    })), silent = TRUE)
    if (inherits(fits, "try-error")) return(fail)
    # Nakagawa marginal R^2: fixed-effect variance over total variance. The
    # conditional-residual version does not track what the LRT tests, because
    # adding x can move the variance components without moving the residuals.
    r2_marginal <- function(m) {
      vf <- stats::var(as.vector(lme4::getME(m, "X") %*% lme4::fixef(m)))
      vr <- sum(vapply(lme4::VarCorr(m),
                       function(v) sum(diag(as.matrix(v))), numeric(1)))
      ve <- stats::sigma(m)^2
      den <- vf + vr + ve
      if (!is.finite(den) || den <= 0) NA_real_ else vf / den
    }
    r2 <- r2_marginal(fits$m1) - r2_marginal(fits$m0)
    k <- length(lme4::fixef(fits$m1)) - length(lme4::fixef(fits$m0))
    p_col <- grep("^Pr\\(>Chisq\\)$", colnames(fits$an), value = TRUE)
    if (!length(p_col) || k <= 0) return(fail)
    dfr <- length(y) - nlevels(g) - k
    if (dfr < 5 || !is.finite(r2)) return(fail)
    list(
      r2      = r2,
      adj_r2  = NA_real_,
      p       = unname(fits$an[[p_col[1L]]][2L]),
      n       = length(y), n_group = nlevels(g), unit = obs_unit
    )
  }
}


## ---- Small summary helpers that tolerate all-NA groups ----------------------
first_unit <- function(u) { u <- u[!is.na(u)]; if (length(u)) u[1L] else NA_character_ }
safe_max   <- function(v) { v <- v[!is.na(v)]; if (length(v)) max(v) else NA_real_ }
safe_at_max <- function(lab, v) {
  ok <- !is.na(v); if (!any(ok)) return(NA_character_)
  as.character(lab[ok][which.max(v[ok])])
}


#####################################################################
# PART 1 -- DONOR- AND SAMPLE-LEVEL ASSOCIATION WITH THE EMBEDDING DIMS
#####################################################################

# depth correlation per dimension -- the quantity ArchR's corCutOff acts on
depth_cor_tbl <- list()
assoc_results <- list()

for (emb in embeddings_to_test) {

  # corCutOff = 1: retain everything the 0.75 default would have removed
  reddim <- getReducedDims(project, reducedDims = emb, corCutOff = 1)

  cell_sample <- factor(getCellColData(project, select = "Sample", drop = TRUE))
  cell_qc <- as.data.frame(getCellColData(project,
    select = c("TSSEnrichment", "nFrags", "FRIP")))
  cell_qc$log10_nFrags <- log10(cell_qc$nFrags)
  cell_depth <- cell_qc$log10_nFrags

  # cell-level: each QC metric varies per cell, so it is correlated per cell
  abscor <- function(v, w) abs(stats::cor(v, w, use = "complete.obs"))
  dcor      <- apply(reddim, 2, abscor, w = cell_depth)
  tsscor    <- apply(reddim, 2, abscor, w = cell_qc$TSSEnrichment)
  fripcor   <- apply(reddim, 2, abscor, w = cell_qc$FRIP)
  # the same correlations after removing the sample means, i.e. purely
  # within-sample, which is where the per-sample QC means cannot look
  demean <- function(v) v - ave(v, cell_sample)
  wdepth <- demean(cell_depth); wtss <- demean(cell_qc$TSSEnrichment)
  wfrip  <- demean(cell_qc$FRIP)
  dcor_w    <- apply(reddim, 2, function(v) abscor(demean(v), wdepth))
  tsscor_w  <- apply(reddim, 2, function(v) abscor(demean(v), wtss))
  fripcor_w <- apply(reddim, 2, function(v) abscor(demean(v), wfrip))

  depth_cor_tbl[[emb]] <- data.frame(
    embedding = emb, dim = colnames(reddim), dim_index = seq_along(dcor),
    abs_cor_with_log10_nFrags = round(unname(dcor), 3),
    abs_cor_with_TSSEnrichment = round(unname(tsscor), 3),
    abs_cor_with_FRIP = round(unname(fripcor), 3),
    within_sample_abs_cor_log10_nFrags = round(unname(dcor_w), 3),
    within_sample_abs_cor_TSSEnrichment = round(unname(tsscor_w), 3),
    within_sample_abs_cor_FRIP = round(unname(fripcor_w), 3),
    excluded_by_default_cutoff = unname(dcor) > archr_cor_cutoff,
    stringsAsFactors = FALSE
  )
  message(sprintf(
    "%s cell-level QC correlation, max over dimensions: nFrags %.2f, TSS %.2f, FRIP %.2f (within-sample: %.2f / %.2f / %.2f)",
    emb, max(dcor), max(tsscor), max(fripcor),
    max(dcor_w), max(tsscor_w), max(fripcor_w)
  ))
  message(sprintf(
    "%s: %d/%d dimensions exceed the default corCutOff of %.2f with log10(nFrags) (%s)",
    emb, sum(dcor > archr_cor_cutoff), length(dcor), archr_cor_cutoff,
    if (any(dcor > archr_cor_cutoff))
      paste(colnames(reddim)[dcor > archr_cor_cutoff], collapse = ", ") else "none"
  ))

  # ---- per-sample mean of each dimension (these ARE the LSI/Harmony dims) ----
  sums         <- rowsum(reddim, group = cell_sample)
  n_per_sample <- as.numeric(table(cell_sample)[rownames(sums)])
  sample_emb   <- sums / n_per_sample                    # samples x dims
  k            <- min(n_dims, ncol(sample_emb))
  dims         <- sample_emb[, seq_len(k), drop = FALSE]
  colnames(dims) <- paste0("Dim", seq_len(k))

  # weight = each dimension's share of between-sample variance
  var_between <- apply(dims, 2, stats::var)
  var_share   <- var_between / sum(var_between)
  var_cells   <- apply(reddim[, seq_len(k), drop = FALSE], 2, stats::var)
  var_share_cells <- var_cells / sum(apply(reddim, 2, stats::var))

  # align covariates to the embedding sample order (join on arrow_name)
  cv <- covar[match(rownames(dims), covar$arrow_name), ]
  if (anyNA(cv$arrow_name)) {
    warning(sprintf(
      "%s: %d sample(s) in the embedding not found in Table S1",
      emb, sum(is.na(cv$arrow_name))
    ))
  }

  covariate_cols <- list(
    Cohort            = cv$exposure_type,      # donor-level
    Control_status    = cv$control_status,     # donor-level
    Exposure_group    = cv$exposure_group,     # donor-level
    Age               = cv$Age_numeric,        # donor-level
    Sex_observed      = cv$Sex,                # donor-level
    Sex_predicted     = cv$Sex_predicted,      # donor-level (chrY/XIST call)
    Days_from_onset   = cv$sampling_day,       # varies within donor
    QC_nCells         = cv$n_cells,            # varies within donor
    QC_meanTSS        = cv$mean_TSS,
    QC_meanLog10Frags = cv$mean_log10_nFrags,
    QC_meanFRIP       = cv$mean_FRIP
  )

  res <- do.call(rbind, lapply(seq_len(k), function(i) {
    do.call(rbind, lapply(names(covariate_cols), function(cn) {
      a <- assoc_test(dims[, i], covariate_cols[[cn]], cv$donor_id)
      adj <- if (cn == "Cohort") {
        list(r2 = NA_real_, adj_r2 = NA_real_, p = NA_real_,
             n = a$n, n_group = a$n_group, unit = NA_character_)
      } else {
        assoc_test(dims[, i], covariate_cols[[cn]], cv$donor_id, cv$exposure_type)
      }
      data.frame(
        embedding = emb,
        Dim = paste0("Dim", i),
        dim_index = i,
        var_share_between_samples = var_share[i],
        var_share_cells = var_share_cells[i],
        depth_cor = unname(dcor[i]),
        covariate = cn,
        r2 = a$r2, adj_r2 = a$adj_r2, p = a$p,
        r2_partial_adjCohort = adj$r2, adj_r2_partial_adjCohort = adj$adj_r2,
        p_adjCohort = adj$p,
        n = a$n, n_group = a$n_group, unit = a$unit,
        stringsAsFactors = FALSE
      )
    }))
  }))
  assoc_results[[emb]] <- res
}

depth_cor_df <- dplyr::bind_rows(depth_cor_tbl)
write.csv(depth_cor_df,
  file = file.path(repo_dir, "figures/dim_depth_correlation.csv"), row.names = FALSE)

assoc_df <- dplyr::bind_rows(assoc_results)

# BH-adjust within each embedding x covariate family
assoc_df <- assoc_df %>%
  dplyr::group_by(embedding, covariate) %>%
  dplyr::mutate(
    p_adj = p.adjust(p, method = "BH"),
    p_adjCohort_bh = p.adjust(p_adjCohort, method = "BH")
  ) %>%
  dplyr::ungroup()

## ---- Variance-weighted summary ----------------------------------------------
# sum over dimensions of (adjusted R^2) x (share of between-sample variance)
var_summary <- assoc_df %>%
  dplyr::group_by(embedding, covariate) %>%
  dplyr::summarise(
    unit = first_unit(unit),
    total_var_marginal = sum(pmax(r2, 0) * var_share_between_samples, na.rm = TRUE),
    total_var_adjCohort = sum(pmax(r2_partial_adjCohort, 0) * var_share_between_samples,
                              na.rm = TRUE),
    max_r2 = safe_max(r2),
    max_r2_Dim = safe_at_max(Dim, r2),
    .groups = "drop"
  ) %>%
  dplyr::arrange(embedding, dplyr::desc(total_var_marginal))
message("Total between-sample variance associated with each covariate:")
print(as.data.frame(var_summary))
write.csv(var_summary,
  file = file.path(repo_dir, "figures/dim_covariate_variance_summary.csv"),
  row.names = FALSE
)

write.csv(assoc_df,
  file = file.path(repo_dir, "figures/dim_covariate_association.csv"),
  row.names = FALSE
)


#####################################################################
# PART 1B -- THE SAME TESTS RUN AT THE CELL LEVEL, FOR COMPARISON
#####################################################################
# Repeats every donor-level test with the sample value copied onto each of its
# cells, which is the alternative that gets proposed for cell-level embeddings.
# The point is the comparison: the effect size barely moves while the p-value
# collapses, because the extra rows are copies rather than observations.
# Set SKIP_CELL_LEVEL_DEMO=1 to skip it.
do_cell_demo <- !nzchar(Sys.getenv("SKIP_CELL_LEVEL_DEMO"))

if (do_cell_demo) {
  demo_rows <- list()
  for (emb in embeddings_to_test) {
    reddim      <- getReducedDims(project, reducedDims = emb, corCutOff = 1)
    cell_sample <- as.character(getCellColData(project, select = "Sample", drop = TRUE))
    ci          <- match(cell_sample, covar$arrow_name)
    kk          <- min(n_dims, ncol(reddim))

    demo_covs <- list(
      Cohort          = covar$exposure_type[ci],
      Control_status  = covar$control_status[ci],
      Exposure_group  = covar$exposure_group[ci],
      Age             = covar$Age_numeric[ci],
      Sex_observed    = covar$Sex[ci],
      Sex_predicted   = covar$Sex_predicted[ci],
      Days_from_onset = covar$sampling_day[ci]
    )

    for (i in seq_len(kk)) {
      y <- reddim[, i]
      for (cn in names(demo_covs)) {
        x  <- demo_covs[[cn]]
        ok <- !is.na(y) & !is.na(x)
        if (sum(ok) < 10) next
        xx <- if (is.factor(x)) droplevels(x[ok]) else as.numeric(x[ok])
        if (is.factor(xx) && nlevels(xx) < 2) next
        if (!is.factor(xx) && stats::sd(xx) == 0) next
        fit <- stats::lm(y[ok] ~ xx)
        fs  <- summary(fit)$fstatistic
        demo_rows[[length(demo_rows) + 1L]] <- data.frame(
          embedding = emb, Dim = paste0("Dim", i), covariate = cn,
          r2_cell = summary(fit)$r.squared,
          p_cell  = if (is.null(fs)) NA_real_ else
            stats::pf(fs[1L], fs[2L], fs[3L], lower.tail = FALSE),
          n_cells = sum(ok), stringsAsFactors = FALSE
        )
      }
    }
    message(sprintf("%s: cell-level comparison fitted on %d cells", emb, length(y)))
  }

  demo_df <- dplyr::bind_rows(demo_rows) %>%
    dplyr::group_by(embedding, covariate) %>%
    dplyr::mutate(p_cell_adj = p.adjust(p_cell, method = "BH")) %>%
    dplyr::ungroup()

  cmp <- assoc_df %>%
    dplyr::select(embedding, Dim, covariate, r2, p_adj, n, n_group, unit) %>%
    dplyr::inner_join(demo_df, by = c("embedding", "Dim", "covariate")) %>%
    dplyr::rename(r2_grouped = r2, p_grouped = p_adj,
                  n_grouped = n, n_groups = n_group, unit_grouped = unit)

  write.csv(cmp,
    file = file.path(repo_dir, "figures/celllevel_vs_grouped_comparison.csv"),
    row.names = FALSE)

  message("Cell-level vs grouped testing, median over dimensions and covariates:")
  print(cmp %>%
    dplyr::group_by(embedding) %>%
    dplyr::summarise(
      median_r2_grouped   = round(median(r2_grouped, na.rm = TRUE), 3),
      median_r2_cell      = round(median(r2_cell, na.rm = TRUE), 3),
      median_p_grouped    = signif(median(p_grouped, na.rm = TRUE), 3),
      median_p_cell       = signif(median(p_cell_adj, na.rm = TRUE), 3),
      sig_grouped         = sum(p_grouped < 0.05, na.rm = TRUE),
      sig_cell            = sum(p_cell_adj < 0.05, na.rm = TRUE),
      n_tests             = dplyr::n(), .groups = "drop") %>%
    as.data.frame())

  # Paired plot: one segment per covariate x dimension joining the two schemes.
  # Segments are near-vertical (the effect size hardly moves) and very long
  # (the p-value does), which is the whole argument in one picture.
  demo_cap <- 300
  nl10 <- function(p) pmin(-log10(pmax(p, .Machine$double.xmin)), demo_cap)
  pair_id <- paste(cmp$embedding, cmp$Dim, cmp$covariate, sep = "|")
  cmp_long <- rbind(
    data.frame(pair = pair_id, embedding = cmp$embedding, covariate = cmp$covariate,
      scheme = "grouped (donor / sample)", r2 = cmp$r2_grouped,
      negp = nl10(cmp$p_grouped), stringsAsFactors = FALSE),
    data.frame(pair = pair_id, embedding = cmp$embedding, covariate = cmp$covariate,
      scheme = "cell level (sample value repeated per cell)", r2 = cmp$r2_cell,
      negp = nl10(cmp$p_cell_adj), stringsAsFactors = FALSE)
  )
  cmp_long$scheme <- factor(cmp_long$scheme,
    levels = c("grouped (donor / sample)", "cell level (sample value repeated per cell)"))

  p_demo <- ggplot(cmp_long, aes(x = r2, y = negp)) +
    geom_line(aes(group = pair), colour = "grey80", linewidth = 0.3) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed",
               colour = "grey40", linewidth = 0.4) +
    geom_point(aes(colour = scheme), size = 1.7, alpha = 0.85) +
    scale_colour_manual(values = c("grouped (donor / sample)" = "#3C5488",
      "cell level (sample value repeated per cell)" = "#E64B35"), name = NULL) +
    facet_wrap(~embedding) +
    theme_classic(base_size = 12) +
    theme(legend.position = "top", plot.subtitle = element_text(size = 8)) +
    labs(
      title = "Testing the same covariates at the cell level inflates significance without changing effect size",
      subtitle = paste0(
        "One grey segment per covariate and dimension, joining the two ways of testing it. ",
        "Dashed line, adjusted p = 0.05.\n",
        "Segments are near-vertical: repeating a sample's value across its cells leaves the variance ",
        "explained almost unchanged\nwhile moving the p-value by tens of orders of magnitude, ",
        "because the added rows are copies rather than independent observations."
      ),
      x = "Variance explained (R\u00b2)",
      y = expression(-log[10] ~ "(BH-adjusted p)")
    )
  ggsave(p_demo,
    filename = file.path(repo_dir, "figures/celllevel_vs_grouped_comparison.pdf"),
    width = 10, height = 5.5, units = "in", dpi = 300)
}


#####################################################################
# PART 2 -- PSEUDOBULK (SAMPLE x CELL TYPE) ASSOCIATION
#####################################################################
# QC is recomputed per pseudobulk and the grouping is the SAMPLE, so
# cohort/age/sex collapse back instead of repeating across ~10 pseudobulks.

min_pb_cells <- 20
assoc_ct_results <- list()

for (emb in embeddings_to_test) {
  reddim <- getReducedDims(project, reducedDims = emb, corCutOff = 1)
  cd <- as.data.frame(getCellColData(project,
    select = c("Sample", ct_col, "TSSEnrichment", "nFrags", "FRIP")))
  names(cd)[names(cd) == ct_col] <- "CellType"
  cd$log10_nFrags <- log10(cd$nFrags)

  ok  <- !is.na(cd$CellType)
  grp <- paste(cd$Sample, cd$CellType, sep = "||")

  sums   <- rowsum(reddim[ok, , drop = FALSE], group = grp[ok])
  n_grp  <- as.numeric(table(grp[ok])[rownames(sums)])
  pb_emb <- sums / n_grp
  keep   <- n_grp >= min_pb_cells
  pb_emb <- pb_emb[keep, , drop = FALSE]
  pb_ncells <- n_grp[keep]

  # QC summarised at the SAME unit as the embedding rows
  qc_pb <- cd[ok, ] %>%
    dplyr::mutate(.grp = grp[ok]) %>%
    dplyr::group_by(.grp) %>%
    dplyr::summarise(
      pb_TSS = mean(TSSEnrichment, na.rm = TRUE),
      pb_log10_nFrags = mean(log10_nFrags, na.rm = TRUE),
      pb_FRIP = mean(FRIP, na.rm = TRUE), .groups = "drop"
    )
  qc_pb <- qc_pb[match(rownames(pb_emb), qc_pb$.grp), ]
  stopifnot(identical(qc_pb$.grp, rownames(pb_emb)))

  kk   <- min(n_dims, ncol(pb_emb))
  dims <- pb_emb[, seq_len(kk), drop = FALSE]
  colnames(dims) <- paste0("Dim", seq_len(kk))
  var_between <- apply(dims, 2, stats::var)
  var_share   <- var_between / sum(var_between)

  key       <- do.call(rbind, strsplit(rownames(dims), "\\|\\|"))
  pb_sample <- key[, 1]
  pb_celltype <- key[, 2]
  cvp       <- covar[match(pb_sample, covar$arrow_name), ]

  covariate_cols <- list(
    CellType          = factor(pb_celltype),   # varies within sample
    Cohort            = cvp$exposure_type,     # constant within sample
    Control_status    = cvp$control_status,
    Exposure_group    = cvp$exposure_group,
    Age               = cvp$Age_numeric,
    Sex_observed      = cvp$Sex,
    Sex_predicted     = cvp$Sex_predicted,
    Days_from_onset   = cvp$sampling_day,
    QC_nCells         = pb_ncells,             # per pseudobulk
    QC_meanTSS        = qc_pb$pb_TSS,          # per pseudobulk
    QC_meanLog10Frags = qc_pb$pb_log10_nFrags, # per pseudobulk
    QC_meanFRIP       = qc_pb$pb_FRIP          # per pseudobulk
  )

  res <- do.call(rbind, lapply(seq_len(kk), function(i) {
    do.call(rbind, lapply(names(covariate_cols), function(cn) {
      # group = sample: sample-constant covariates collapse, others get (1|sample)
      a <- assoc_test(dims[, i], covariate_cols[[cn]], pb_sample,
                      collapsed_unit = "per sample",
                      obs_unit = "per pseudobulk, sample-adjusted")
      data.frame(
        embedding = emb, Dim = paste0("Dim", i), dim_index = i,
        var_share_between_pseudobulks = var_share[i], covariate = cn,
        r2 = a$r2, adj_r2 = a$adj_r2, p = a$p,
        n = a$n, n_group = a$n_group, unit = a$unit,
        stringsAsFactors = FALSE
      )
    }))
  }))
  assoc_ct_results[[emb]] <- res
}

assoc_ct_df <- dplyr::bind_rows(assoc_ct_results) %>%
  dplyr::group_by(embedding, covariate) %>%
  dplyr::mutate(p_adj = p.adjust(p, method = "BH")) %>%
  dplyr::ungroup()

write.csv(assoc_ct_df,
  file = file.path(repo_dir, "figures/dim_covariate_association_by_celltype.csv"),
  row.names = FALSE
)

var_summary_ct <- assoc_ct_df %>%
  dplyr::group_by(embedding, covariate) %>%
  dplyr::summarise(
    unit = first_unit(unit),
    total_var = sum(pmax(r2, 0) * var_share_between_pseudobulks, na.rm = TRUE),
    max_r2 = safe_max(r2),
    max_r2_Dim = safe_at_max(Dim, r2),
    .groups = "drop"
  ) %>%
  dplyr::arrange(embedding, dplyr::desc(total_var))
message("Total pseudobulk-level variance associated with each covariate:")
print(as.data.frame(var_summary_ct))
write.csv(var_summary_ct,
  file = file.path(repo_dir, "figures/dim_covariate_variance_summary_by_celltype.csv"),
  row.names = FALSE
)




#####################################################################
# FIGURES -- association heatmaps
#####################################################################
# p can underflow to 0 and -log10(0) is Inf, which ggplot paints as missing.
p_cap <- 50
neglog10 <- function(p, cap = p_cap) {
  pmin(-log10(pmax(p, .Machine$double.xmin)), cap)
}
# design covariates first, then donor, then technical
covariate_order <- c(
  "CellType", "Cohort", "Control_status", "Exposure_group",
  "Age", "Sex_observed", "Sex_predicted", "Days_from_onset",
  "QC_nCells", "QC_meanTSS", "QC_meanLog10Frags", "QC_meanFRIP"
)

# label each row with the unit the test was run on
add_row_label <- function(d) {
  u <- d %>%
    dplyr::group_by(covariate) %>%
    dplyr::summarise(unit = first_unit(unit), .groups = "drop")
  d <- dplyr::left_join(d, dplyr::rename(u, unit_lab = unit), by = "covariate")
  d$row_label <- ifelse(is.na(d$unit_lab), d$covariate,
                        paste0(d$covariate, "  [", d$unit_lab, "]"))
  d
}
order_rows <- function(d) {
  lv <- d %>% dplyr::distinct(covariate, row_label)
  lv <- lv[match(intersect(covariate_order, lv$covariate), lv$covariate), ]
  rev(lv$row_label)
}

assoc_plot_df <- add_row_label(assoc_df)

p_assoc <- ggplot(assoc_plot_df, aes(x = Dim, y = row_label, fill = neglog10(p_adj))) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = ifelse(is.na(r2), "", sprintf("%.2f", pmax(r2, 0)))), size = 3) +
  facet_wrap(~embedding, ncol = 1) +
  scale_fill_gradient(
    low = "white", high = "#006400", na.value = "grey95",
    name = paste0("-log10(adj. p)\n(capped at ", p_cap, ")")
  ) +
  scale_x_discrete(limits = paste0("Dim", seq_len(n_dims))) +
  scale_y_discrete(limits = order_rows(assoc_plot_df)) +
  labs(
    title = paste0("Association of ", paste(embeddings_to_test, collapse = " / "),
                   " dimensions with known covariates"),
    subtitle = paste0(
      "Columns are the LSI / Harmony dimensions themselves (per-sample means), not principal components of them.\n",
      "Tile label = variance explained (partial R^2). Shading = BH-adjusted significance, which already\n",
      "accounts for the degrees of freedom each covariate uses. Adjusted R^2 is in Table S1B.\n",
      "Bracketed label = one observation in that test. 'per donor': the covariate is constant within a\n",
      "donor, so repeat draws are collapsed. 'per sample, donor-adjusted': the covariate varies within a\n",
      "donor, so samples are kept separate and a donor random intercept absorbs the repeated sampling."
    ),
    x = NULL, y = NULL
  ) +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p_assoc,
  filename = file.path(repo_dir, "figures/dim_covariate_association.pdf"),
  width = 10,
  height = max(7.5, 0.5 * length(unique(assoc_df$covariate)) + 3),
  units = "in", dpi = 300, limitsize = FALSE
)

assoc_ct_plot_df <- add_row_label(assoc_ct_df)

p_assoc_ct <- ggplot(assoc_ct_plot_df,
  aes(x = Dim, y = row_label, fill = neglog10(p_adj))) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = ifelse(is.na(r2), "", sprintf("%.2f", pmax(r2, 0)))), size = 3) +
  facet_wrap(~embedding, ncol = 1) +
  scale_fill_gradient(low = "white", high = "#006400", na.value = "grey95",
    name = paste0("-log10(adj. p)\n(capped at ", p_cap, ")")) +
  scale_x_discrete(limits = paste0("Dim", seq_len(n_dims))) +
  scale_y_discrete(limits = order_rows(assoc_ct_plot_df)) +
  labs(
    title = "Association of pseudobulk (sample x cell type) dimensions with covariates",
    subtitle = paste0(
      "Pseudobulks of >= ", min_pb_cells, " cells. Tile label = variance explained (partial R^2); shading = BH-adjusted significance.\n",
      "Cell type and the QC metrics vary within a sample and are tested per pseudobulk with a sample random intercept;\n",
      "cohort, age and sex are constant within a sample and are collapsed back, so they are not counted once per pseudobulk."
    ),
    x = NULL, y = NULL
  ) +
  theme_classic() +
  theme(plot.subtitle = element_text(size = 8),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p_assoc_ct,
  filename = file.path(repo_dir, "figures/dim_covariate_association_by_celltype.pdf"),
  width = 10, height = 8, units = "in", dpi = 300
)


#####################################################################
# SUPPLEMENTARY TABLE OUTPUT
#####################################################################
assoc_supp <- assoc_df %>%
  dplyr::transmute(
    Embedding = embedding,
    Dimension = Dim,
    `Share of between-sample variance (%)` = round(100 * var_share_between_samples, 2),
    `|cor| with log10(nFrags), cell level` = round(depth_cor, 3),
    Covariate = covariate,
    `Unit of observation` = unit,
    `R2 (marginal)` = round(r2, 3),
    `Adjusted R2 (marginal)` = round(adj_r2, 3),
    `p (marginal, BH)` = signif(p_adj, 3),
    `Partial R2 (cohort-adjusted)` = round(r2_partial_adjCohort, 3),
    `Adjusted partial R2 (cohort-adjusted)` = round(adj_r2_partial_adjCohort, 3),
    `p (cohort-adjusted, BH)` = signif(p_adjCohort_bh, 3),
    `n observations` = n,
    `n groups` = n_group
  ) %>%
  dplyr::arrange(
    Embedding, factor(Dimension, levels = paste0("Dim", seq_len(n_dims))), Covariate
  )

wb <- openxlsx::loadWorkbook(suppTables) # keeps all existing sheets intact

# fix the influenza timepoint mislabel: the final timepoint is day 28, not 30
s1_now  <- openxlsx::readWorkbook(wb, sheet = "Table S1")
rec_col <- which(names(s1_now) == "record_id")
if (length(rec_col) == 1) {
  fixed_rec <- sub("day 30 after challenge", "day 28 after challenge",
    as.character(s1_now$record_id)
  )
  n_fixed <- sum(fixed_rec != as.character(s1_now$record_id), na.rm = TRUE)
  openxlsx::writeData(wb, "Table S1", x = fixed_rec, startCol = rec_col, startRow = 2)
  message("Corrected 'day 30' -> 'day 28' in Table S1 record_id for ", n_fixed, " row(s)")
}

new_cols <- data.frame(
  # Donor_ID: HIV is longitudinal and some C19 donors were sampled twice
  Donor_ID = covar$donor_id,
  Age_reported = covar$Age_reported,
  Sex_reported = as.character(covar$Sex),
  Sampling_day_rel_onset = covar$sampling_day,
  Onset_reference = covar$onset_reference,
  Viral_load = covar$viral_load,
  # processing batch, C19 only (from sample_annots/ATAC_metadata_covid.csv)
  Processing_batch = as.character(covar$processing_batch),
  # full sex-call detail stays in figures/sex_prediction_metrics.csv
  Sex_predicted = as.character(covar$Sex_predicted),
  P_male = covar$p_male,
  chrY_fragments = covar$chrY,
  chrX_fragments = covar$chrX,
  XIST_fragments = covar$xist,
  chrY_fraction = signif(covar$chrY_frac, 4),
  XIST_fraction = signif(covar$xist_frac, 4),
  chrY_chrX_ratio = signif(covar$chrY_chrX_ratio, 4),
  Sex_metadata_flag = covar$sex_flag,
  stringsAsFactors = FALSE
)
# cells behind each sample x cell-type pseudobulk, one column per cell type
new_cols <- cbind(new_cols, as.data.frame(covar)[, pb_count_cols, drop = FALSE])

# GUARD: written at a column OFFSET, so covar must stay row-aligned with the sheet
if (nrow(covar) != nrow(s1_now)) {
  stop(sprintf(
    "Table S1 has %d rows but covar has %d: a join above changed the row count. Not writing.",
    nrow(s1_now), nrow(covar)
  ))
}
if ("arrow_name" %in% names(s1_now) &&
    !identical(as.character(covar$arrow_name), as.character(s1_now$arrow_name))) {
  stop("covar is no longer in Table S1 row order; the positional append would misalign. Not writing.")
}
# refuse rather than appending a second copy of the same columns
clash <- intersect(names(new_cols), names(s1_now))
if (length(clash)) {
  stop("Table S1 already contains: ", paste(clash, collapse = ", "),
       ". Start from the pristine All_Supplementary_Tables.xlsx rather than an ",
       "already-updated copy, or remove those columns first.")
}

s1_ncol <- ncol(s1_now)
openxlsx::writeData(wb, "Table S1", x = new_cols, startCol = s1_ncol + 1, startRow = 1)
message("Added ", ncol(new_cols), " column(s) to Table S1 (",
  length(pb_count_cols), " of them per-cell-type pseudobulk counts)")

write_sheet <- function(wb, name, df) {
  if (name %in% names(wb)) openxlsx::removeWorksheet(wb, name)
  openxlsx::addWorksheet(wb, name)
  openxlsx::writeData(wb, name, df, withFilter = TRUE,
    headerStyle = openxlsx::createStyle(textDecoration = "bold"))
  openxlsx::setColWidths(wb, name, cols = seq_along(df), widths = "auto")
}

write_sheet(wb, "Table S1B", assoc_supp)          # per-dimension associations
write_sheet(wb, "Table S1C", as.data.frame(depth_cor_df))  # depth correlations
if (exists("confound_rows") && !is.null(confound_rows)) {
  write_sheet(wb, "Table S1D", as.data.frame(confound_rows))               # confound structure
}
write_sheet(wb, "Table S1E", sex_lda_metrics)                              # sex-classifier performance

out_xlsx <- file.path(repo_dir, "sample_annots/All_Supplementary_Tables_updated.xlsx")
openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message("Wrote association results to: ", out_xlsx)


#####################################################################
# CONSOLE SUMMARY
#####################################################################
message("Significant MARGINAL dimension-covariate associations (adj. p < 0.05):")
print(subset(assoc_df, p_adj < 0.05,
  select = c(embedding, Dim, covariate, adj_r2, p_adj, n, n_group, unit)))

# does any covariate survive once cohort is controlled for?
message("Significant COHORT-ADJUSTED associations (non-cohort covariates, adj. p < 0.05):")
print(subset(assoc_df, covariate != "Cohort" & p_adjCohort_bh < 0.05,
  select = c(embedding, Dim, covariate, adj_r2_partial_adjCohort, p_adjCohort_bh, n, unit)))

message("Cohort association before vs after Harmony (total between-sample variance):")
print(as.data.frame(subset(var_summary, covariate == "Cohort",
  select = c(embedding, covariate, total_var_marginal, max_r2, max_r2_Dim))))

#####################################################################
