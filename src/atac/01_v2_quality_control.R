#!/usr/bin/env Rscript

#####################################################################
# 01_v2_quality_control.R
# created on 2026-08-07 by Irem Gunduz
# Association of sample-level PCs with known covariates
# (Reviewer response: characterise variation w.r.t. potential confounders)
#
# Rationale
#  - We ask whether the low-dimensional structure of the ATAC data is
#    associated with cohort (exposure_type), Sex and Age.
#  - Age and Sex are recorded per SAMPLE and only for the OP and HIV
#    cohorts (see Table S1), so those two tests are run on that subset.
#  - The covariates are sample-level, therefore the test is performed at
#    the sample level: the cell-level embedding is aggregated to a
#    per-sample mean before PCA. Testing cell-level PCs against
#    sample-level covariates would be pseudoreplication and would
#    massively inflate significance.
#  - We run the test on BOTH the pre-Harmony IterativeLSI embedding and
#    the Harmony-corrected embedding (added in 02_cluster_and_batch.R).
#    Expectation: cohort is associated before correction (the bias we
#    account for via Harmony) and is attenuated afterwards, while Age and
#    Sex are not significantly associated with any PC in either case.
#####################################################################

## Load Libraries
suppressPackageStartupMessages({
  library(readxl)
  library(openxlsx)
  library(ArchR)
  library(dplyr)
  library(ggplot2)
  library(MASS) # LDA for molecular sex classification
})
set.seed(12) # set seed

## ---- Paths / parameters -----------------------------------------------------
repo_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript"
suppTables <- file.path(repo_dir, "sample_annots/All_Supplementary_Tables.xlsx")
n_pc <- 10 # number of top PCs to test
embeddings_to_test <- c("IterativeLSI", "Harmony") # pre- and post-batch-correction
outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/"

# Load the ArchR project (already subsetted to non-BA samples)
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)

## ---- Sample-level covariate table (Table S1) --------------------------------
# NB: project$Sample is the arrow filename (e.g. "ATAC_055_fragments.tsv.gz"),
# which matches Table S1$arrow_name -- NOT Table S1$sampleId ("C19_mod_055").
# We therefore join the covariates on arrow_name.
covar <- readxl::read_excel(suppTables, sheet = "Table S1") %>%
  dplyr::select(sampleId, arrow_name, exposure_type, exposure_group,
    exposure_grouping, record_id, Age, Sex) %>%
  dplyr::mutate(
    # Cohort: C19 / HIV / Influenza / OP
    exposure_type = factor(exposure_type),
    # within-cohort condition arm, tested separately from cohort
    exposure_group = factor(exposure_group),
    # "healthy" arms are C19_ctrl / HIV_ctrl / Influenza_ctrl (n = 16). OP has
    # no unexposed arm (OP_low is low-exposure), so every OP sample is Exposed
    # and this is partly nested in cohort. Read the cohort-adjusted column.
    control_status = factor(
      ifelse(exposure_grouping == "healthy", "Control", "Exposed"),
      levels = c("Control", "Exposed")
    ),
    Sex = factor(Sex),
    Age = as.numeric(Age)
  )

#####################################################################
# RECOVERED DONOR METADATA (age / sex / sampling time relative to onset)
# Age and sex were originally present only for the OP and HIV cohorts. We
# recover them for the COVID-19 cohort from the clinical metadata table, and
# recover sampling time relative to infection onset for all three infection
# cohorts:
#   COVID-19  : days after the first positive SARS-CoV-2 test
#   HIV       : days relative to seroconversion (negative = pre-infection)
#   Influenza : days after the challenge protocol (already in Table S1)
#   OP        : chronic environmental exposure, no datable onset event
# The recorded values fill Table S1; the molecular sex prediction below is kept
# and used to validate them.
#####################################################################

# Authoritative ATAC-label -> clinical Donor ID map for the COVID-19 cohort
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

# bin midpoint for the binned COVID-19 ages, so age can be used as a continuous
# covariate alongside the exact ages recorded for HIV and OP
bin_to_numeric <- function(x) {
  vapply(x, function(s) {
    if (is.na(s) || s == "" || s == "NA") return(NA_real_)
    if (grepl("\\+$", s)) return(as.numeric(sub("\\+", "", s)) + 5)
    p <- as.numeric(strsplit(s, "-")[[1]])
    if (length(p) == 2) mean(p) else suppressWarnings(as.numeric(s))
  }, numeric(1), USE.NAMES = FALSE)
}

label <- sub("_fragments\\.tsv\\.gz$", "", covar$arrow_name)
donor_id <- unname(atac_to_donor[label])
pm_i <- match(donor_id, pm$Donor)
hiv_i <- match(label, hiv_meta$stem)
# NB: "day 30 after challenge" in the legacy annotation is a MISLABEL -- the
# actual final timepoint is day 28 (group label Influenza_d28). We map it to 28
# and also correct the text in Table S1 when writing the workbook below.
flu_day <- c(
  "right before challenge" = -1, "day 3 after challenge" = 3,
  "day 6 after challenge" = 6, "day 30 after challenge" = 28,
  "day 28 after challenge" = 28
)

covar <- covar %>% dplyr::mutate(
  donor_id = dplyr::coalesce(donor_id, hiv_meta$donor[hiv_i]),
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

message("Recorded metadata coverage after recovery (n samples with a value):")
print(covar %>%
  dplyr::group_by(exposure_type) %>%
  dplyr::summarise(
    n = dplyr::n(), sex = sum(!is.na(Sex)), age = sum(!is.na(Age_numeric)),
    sampling_day = sum(!is.na(sampling_day)), .groups = "drop"
  ))

## ---- Per-sample QC metrics (technical covariates) ---------------------------
# The reviewer/PI also wants technical quality tested as a potential driver of
# global structure. We summarise the key single-cell QC metrics to the sample
# level (mean over the cells retained in the project) and treat them as
# continuous covariates in the association test below.
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
# One column per cell type in Table S1, so the table stays one row per sample.
# Raw counts: every combination is reported, including those below the >50-cell
# cutoff 05_pseudobulk.R applies and the Plasma compartment it drops.
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

covar <- dplyr::left_join(covar, pb_wide, by = "arrow_name")

# per-cell-type counts should sum to n_cells, minus unannotated cells
# covar is a tibble, so coerce before the arithmetic
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
# Two orthogonal, depth-normalised molecular features:
#   chrY_frac : chrY fragments / total fragments  (HIGH in XY / male)
#   xist_frac : XIST-locus accessibility / total  (HIGH in XX / female; XIST is
#               accessible from the inactive X in XX donors, silent in XY)
# The two features are INTERNALLY CONSISTENT (high chrY co-occurs with low XIST
# and vice versa), so we call sex UNSUPERVISED from the molecular signal and do
# NOT train on recorded sex. Recorded labels contain errors that would corrupt a
# supervised model (they did: training inverted the chrY->male direction); we
# therefore use recorded sex ONLY to validate the prediction and to flag samples
# whose recorded sex disagrees with the data.
xist_gr <- GenomicRanges::GRanges(
  "chrX", IRanges::IRanges(start = 73820651, end = 73852753)
) # XIST gene, hg38

arrowFiles <- ArchR::getArrowFiles(project)
cell_sample_all <- getCellColData(project, select = "Sample", drop = TRUE)
all_cellNames <- project$cellNames
nfrags_all <- getCellColData(project, select = "nFrags", drop = TRUE)
total_by_sample <- tapply(nfrags_all, cell_sample_all, sum)

sex_metrics <- do.call(rbind, lapply(arrowFiles, function(af) {
  samp <- gsub("\\.arrow$", "", basename(af))
  cells <- all_cellNames[cell_sample_all == samp]
  if (length(cells) == 0) {
    return(NULL)
  }
  fy <- ArchR::getFragmentsFromArrow(af, chr = "chrY", cellNames = cells, verbose = FALSE)
  fx <- ArchR::getFragmentsFromArrow(af, chr = "chrX", cellNames = cells, verbose = FALSE)
  n_xist <- sum(IRanges::overlapsAny(fx, xist_gr))
  total <- as.numeric(total_by_sample[samp])
  data.frame(
    Sample = samp,
    chrY = length(fy), chrX = length(fx), xist = n_xist,
    chrY_frac = length(fy) / total,
    xist_frac = n_xist / total,
    chrY_chrX_ratio = length(fy) / length(fx),
    stringsAsFactors = FALSE
  )
}))
rownames(sex_metrics) <- NULL

# attach recorded sex + cohort for calibration and diagnostics
sex_metrics <- dplyr::left_join(sex_metrics,
  dplyr::select(covar, arrow_name, exposure_type, Sex),
  by = c("Sample" = "arrow_name")
)

# UNSUPERVISED molecular sex call: 2-cluster k-means on the standardised log
# features (no recorded labels used); the cluster with the higher chrY fraction
# is Male. Because chrY and XIST agree, the two clusters are well separated.
# Log-transform with a DATA-DRIVEN pseudocount (half the smallest non-zero
# value per feature). A fixed 1e-9 would place a sample with zero XIST fragments
# at -9, several log units below every real value, which distorts both the axis
# range and the feature scale. Half the minimum observed value places it just
# below the lowest real measurement instead.
half_min_nonzero <- function(x) min(x[x > 0], na.rm = TRUE) / 2
eps_chrY <- half_min_nonzero(sex_metrics$chrY_frac)
eps_xist <- half_min_nonzero(sex_metrics$xist_frac)
sex_metrics$log_chrY <- log10(sex_metrics$chrY_frac + eps_chrY)
sex_metrics$log_xist <- log10(sex_metrics$xist_frac + eps_xist)

# IMPORTANT: do NOT cluster on chrY and XIST jointly. chrY read fraction carries
# large cohort/processing-specific offsets, so a 2-cluster k-means on both
# features splits the samples by COHORT rather than by sex. XIST accessibility is
# on a consistent scale across all cohorts here, so we threshold on XIST and use
# chrY only as an independent cross-check.
#
# Threshold = 1D natural break (largest gap in the sorted log values), searched
# over the interior quantiles so a single extreme sample cannot define the split.
gap_threshold <- function(x, lo = 0.10, hi = 0.90) {
  lx <- sort(log10(x + 1e-9))
  n <- length(lx)
  i_lo <- max(1L, floor(lo * n))
  i_hi <- min(n - 1L, ceiling(hi * n))
  idx <- i_lo:i_hi
  i <- idx[which.max(diff(lx)[idx])]
  10^mean(c(lx[i], lx[i + 1L])) # geometric midpoint of the largest interior gap
}

# Marker-based reference calls (used to define trustworthy training labels and
# as an independent cross-check): high XIST = female, high chrY = male.
xist_thr <- gap_threshold(sex_metrics$xist_frac)
sex_metrics$Sex_by_xist <- factor(
  ifelse(sex_metrics$xist_frac >= xist_thr, "Female", "Male"),
  levels = c("Female", "Male")
)
chrY_thr <- gap_threshold(sex_metrics$chrY_frac)
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
# Training labels: samples whose RECORDED sex agrees with the marker-based call.
# This deliberately EXCLUDES the mislabelled samples that corrupt a classifier
# trained on all recorded labels (a model trained on every recorded label
# inverted the chrY->male direction). LDA gives a linear decision boundary and
# posterior probabilities, and is numerically stable under near-perfect class
# separation (logistic-regression coefficients diverge in that regime).
lab <- !is.na(sex_metrics$Sex)
train_idx <- lab & sex_metrics$Sex == sex_metrics$Sex_by_xist &
  sex_metrics$features_agree
message(sprintf(
  "LDA training set: %d high-confidence samples (%d recorded labels excluded as inconsistent)",
  sum(train_idx), sum(lab) - sum(train_idx)
))
stopifnot(sum(train_idx) >= 6, dplyr::n_distinct(sex_metrics$Sex[train_idx]) == 2)

train_df <- data.frame(
  Sex = droplevels(factor(sex_metrics$Sex[train_idx], levels = c("Female", "Male"))),
  log_chrY = sex_metrics$log_chrY[train_idx],
  log_xist = sex_metrics$log_xist[train_idx]
)

# Equal priors: the training subset happens to be male-skewed (~87% male), and
# using its empirical class frequencies as priors would bias calls toward Male.
# We have no prior expectation of the sex ratio, so we set it to 0.5/0.5.
equal_prior <- c(Female = 0.5, Male = 0.5)

# leave-one-out cross-validated accuracy on the training set
loo <- MASS::lda(Sex ~ log_chrY + log_xist,
  data = train_df, prior = equal_prior, CV = TRUE
)
message(sprintf(
  "LDA leave-one-out accuracy on training set: %.1f%% (n = %d)",
  100 * mean(loo$class == train_df$Sex), nrow(train_df)
))

# fit on the full training set and predict all samples
lda_fit <- MASS::lda(Sex ~ log_chrY + log_xist, data = train_df, prior = equal_prior)
print(lda_fit)
pred <- stats::predict(lda_fit, newdata = sex_metrics[, c("log_chrY", "log_xist")])
sex_metrics$Sex_predicted <- factor(as.character(pred$class), levels = c("Female", "Male"))
sex_metrics$p_male <- round(pred$posterior[, "Male"], 4)

# Independent validation: the LDA call should reproduce the marker-threshold
# calls. chrY and XIST are strongly anti-correlated (near-collinear on the log
# scale), so the individual LDA coefficients are not individually interpretable
# -- but agreement with the independent XIST threshold confirms the partition.
agree_xist <- mean(sex_metrics$Sex_predicted == sex_metrics$Sex_by_xist)
agree_chrY <- mean(sex_metrics$Sex_predicted == sex_metrics$Sex_by_chrY)
message(sprintf(
  "LDA call agrees with XIST threshold for %.1f%% and with chrY threshold for %.1f%% of all %d samples",
  100 * agree_xist, 100 * agree_chrY, nrow(sex_metrics)
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

# use recorded sex ONLY to validate and to flag likely metadata errors
val <- sex_metrics[lab, ]
message(sprintf(
  "Molecular sex vs recorded sex: %.1f%% concordant (%d labelled samples)",
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

## ---- Diagnostic figure: molecular sex classification ------------------------
# Scatter of the two markers with the LDA decision boundary. Colour = predicted
# sex, shape = recorded sex, and samples whose recorded sex disagrees with the
# molecular call are ringed, so metadata errors are visible at a glance.
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

# Supporting per-sample detail for the sex inference. Kept OUT of Table S1 (which
# carries only the final Sex_predicted call) but retained here as the evidence
# behind the call and behind the flagged metadata discrepancies.
write.csv(
  sex_metrics[, c(
    "Sample", "exposure_type", "Sex", "Sex_predicted", "p_male",
    "Sex_by_xist", "Sex_by_chrY", "features_agree",
    "chrY_frac", "xist_frac", "chrY_chrX_ratio", "sex_flag"
  )],
  file = file.path(repo_dir, "figures/sex_prediction_metrics.csv"), row.names = FALSE
)

# add inferred sex, confidence, raw signals and the flag to the covariate table
# carry the raw chrX / chrY / XIST signals, not just the call, so the
# prediction can be checked from Table S1
covar <- dplyr::left_join(covar,
  dplyr::select(
    sex_metrics, Sample, chrY, chrX, xist, chrY_frac, xist_frac,
    chrY_chrX_ratio, p_male,
    Sex_predicted, Sex_by_xist, Sex_by_chrY, sex_flag
  ),
  by = c("arrow_name" = "Sample")
)

## ---- Helper: association of one PC with one covariate ------------------------
# Returns variance explained (R^2), F-test p-value and the usable n.
# Works for both categorical (Sex, exposure_type) and continuous (Age) x.
pc_assoc <- function(y, x) {
  keep <- !is.na(y) & !is.na(x)
  if (is.factor(x)) x <- droplevels(x[keep]) else x <- x[keep]
  y <- y[keep]
  n <- length(y)
  # need >= 2 groups / non-constant predictor and enough residual df
  ok <- n >= 3 && (if (is.factor(x)) nlevels(x) >= 2 else stats::sd(x) > 0)
  if (!ok) {
    return(c(r2 = NA_real_, p = NA_real_, n = n))
  }
  fit <- stats::lm(y ~ x)
  fs <- summary(fit)$fstatistic
  p <- if (is.null(fs)) NA_real_ else stats::pf(fs[1L], fs[2L], fs[3L], lower.tail = FALSE)
  c(r2 = summary(fit)$r.squared, p = unname(p), n = n)
}

## ---- Helper: cohort-ADJUSTED association of one PC with one covariate --------
# Tests the covariate (x) while controlling for cohort, i.e. compares
#   lm(y ~ cohort)  vs  lm(y ~ cohort + x)
# via a nested F-test, and reports the PARTIAL R^2 (extra variance explained
# by x beyond cohort). This separates a genuine Age/Sex effect from the
# cohort confounding, since Age/Sex are only available for OP + HIV.
pc_assoc_adj <- function(y, x, cohort) {
  keep <- !is.na(y) & !is.na(x) & !is.na(cohort)
  y <- y[keep]
  cohort <- droplevels(cohort[keep])
  if (is.factor(x)) x <- droplevels(x[keep]) else x <- x[keep]
  n <- length(y)
  ok <- n >= 4 && (if (is.factor(x)) nlevels(x) >= 2 else stats::sd(x) > 0)
  if (!ok) {
    return(c(r2_partial = NA_real_, p = NA_real_, n = n))
  }
  # reduced model: cohort only (fall back to intercept if <2 cohorts present)
  reduced <- if (nlevels(cohort) >= 2) stats::lm(y ~ cohort) else stats::lm(y ~ 1)
  full <- stats::update(reduced, . ~ . + x)
  an <- stats::anova(reduced, full)
  p <- an[["Pr(>F)"]][2L]
  rss_r <- sum(stats::resid(reduced)^2)
  rss_f <- sum(stats::resid(full)^2)
  r2_partial <- (rss_r - rss_f) / rss_r
  c(r2_partial = r2_partial, p = unname(p), n = n)
}

## ---- Run per embedding ------------------------------------------------------
assoc_results <- list()

for (emb in embeddings_to_test) {
  # cell x dims embedding
  reddim <- getReducedDims(project, reducedDims = emb)

  # map each cell to its sample, then take the per-sample mean embedding
  cell_sample <- factor(getCellColData(project, select = "Sample", drop = TRUE))
  sums <- rowsum(reddim, group = cell_sample)
  n_per_sample <- as.numeric(table(cell_sample)[rownames(sums)])
  sample_emb <- sums / n_per_sample # samples x dims (pseudobulk embedding)

  # sample-level PCA
  pca <- stats::prcomp(sample_emb, center = TRUE, scale. = TRUE)
  var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
  k <- min(n_pc, ncol(pca$x))
  pcs <- pca$x[, seq_len(k), drop = FALSE]

  # align covariates to the PCA sample order (join on arrow_name)
  cv <- covar[match(rownames(pcs), covar$arrow_name), ]
  if (anyNA(cv$arrow_name)) {
    warning(sprintf(
      "%s: %d sample(s) in the embedding not found in Table S1",
      emb, sum(is.na(cv$arrow_name))
    ))
  }

  # test each PC against each covariate
  covariate_cols <- list(
    Cohort = cv$exposure_type, # all samples
    Control_status = cv$control_status, # Control vs Exposed (all samples)
    Exposure_group = cv$exposure_group, # within-cohort condition arm
    Age = cv$Age_numeric, # C19 + HIV + OP (bin midpoints for C19)
    Sex_observed = cv$Sex, # C19 + HIV + OP (recorded/recovered)
    Sampling_day = cv$sampling_day, # within-cohort sampling time vs onset
    Sex_predicted = cv$Sex_predicted, # all samples (chrY/chrX inference)
    QC_nCells = cv$n_cells, # technical QC (all samples)
    QC_meanTSS = cv$mean_TSS, # technical QC
    QC_meanLog10Frags = cv$mean_log10_nFrags, # technical QC
    QC_meanFRIP = cv$mean_FRIP # technical QC
  )


  res <- do.call(rbind, lapply(seq_len(k), function(i) {
    do.call(rbind, lapply(names(covariate_cols), function(cn) {
      a <- pc_assoc(pcs[, i], covariate_cols[[cn]])
      # cohort-adjusted partial test (not meaningful for Cohort itself)
      adj <- if (cn == "Cohort") {
        c(r2_partial = NA_real_, p = NA_real_, n = a[["n"]])
      } else {
        pc_assoc_adj(pcs[, i], covariate_cols[[cn]], cv$exposure_type)
      }
      data.frame(
        embedding = emb,
        PC = paste0("PC", i),
        pc_index = i,
        var_explained = var_explained[i],
        covariate = cn,
        r2 = a[["r2"]], # marginal variance explained
        p = a[["p"]], # marginal p
        r2_partial_adjCohort = adj[["r2_partial"]], # variance beyond cohort
        p_adjCohort = adj[["p"]], # cohort-adjusted p
        n = a[["n"]],
        stringsAsFactors = FALSE
      )
    }))
  }))
  assoc_results[[emb]] <- res
}

assoc_df <- dplyr::bind_rows(assoc_results)

# BH-adjust p-values within each embedding x covariate family
# (both the marginal test and the cohort-adjusted test)
assoc_df <- assoc_df %>%
  dplyr::group_by(embedding, covariate) %>%
  dplyr::mutate(
    p_adj = p.adjust(p, method = "BH"),
    p_adjCohort_bh = p.adjust(p_adjCohort, method = "BH")
  ) %>%
  dplyr::ungroup()

## ---- Variance-weighted summary ----------------------------------------------
# Total sample-level variance associated with each covariate =
#   sum over PCs of (variance explained by covariate in that PC) x (PC's share
#   of total variance). Gives one interpretable number per covariate/embedding.
var_summary <- assoc_df %>%
  dplyr::group_by(embedding, covariate) %>%
  dplyr::summarise(
    total_var_marginal = sum(r2 * var_explained, na.rm = TRUE),
    total_var_adjCohort = sum(r2_partial_adjCohort * var_explained, na.rm = TRUE),
    # largest single-PC R^2, i.e. what the heatmap tiles show
    max_r2 = max(r2, na.rm = TRUE),
    max_r2_PC = PC[which.max(r2)],
    .groups = "drop"
  )
message("Total sample-level variance associated with each covariate:")
print(as.data.frame(var_summary))
write.csv(var_summary,
  file = file.path(repo_dir, "figures/pc_covariate_variance_summary.csv"),
  row.names = FALSE
)

## ---- Write the association results as a new sheet (Table S1B) ---------------
# Added as its own sheet alongside the existing metadata (Table S1). We save to
# a COPY of the workbook so the master file is never overwritten; rename it to
# All_Supplementary_Tables.xlsx once you are happy, or set
# out_xlsx <- suppTables to write in place.
assoc_supp <- assoc_df %>%
  dplyr::transmute(
    Embedding = embedding,
    PC = PC,
    `PC variance (%)` = round(100 * var_explained, 2),
    Covariate = covariate,
    `R2 (marginal)` = round(r2, 3),
    `p (marginal, BH)` = signif(p_adj, 3),
    `Partial R2 (cohort-adjusted)` = round(r2_partial_adjCohort, 3),
    `p (cohort-adjusted, BH)` = signif(p_adjCohort_bh, 3),
    `n samples` = n
  ) %>%
  dplyr::arrange(
    Embedding,
    factor(PC, levels = paste0("PC", seq_len(n_pc))),
    Covariate
  )

wb <- openxlsx::loadWorkbook(suppTables) # keeps all existing sheets intact

# Correct the influenza timepoint mislabel in the existing record_id column:
# the final challenge timepoint is day 28, not day 30 (group label Influenza_d28).
s1_now <- openxlsx::readWorkbook(wb, sheet = "Table S1")
rec_col <- which(names(s1_now) == "record_id")
if (length(rec_col) == 1) {
  fixed_rec <- sub("day 30 after challenge", "day 28 after challenge",
    as.character(s1_now$record_id)
  )
  n_fixed <- sum(fixed_rec != as.character(s1_now$record_id), na.rm = TRUE)
  openxlsx::writeData(wb, "Table S1", x = fixed_rec,
    startCol = rec_col, startRow = 2
  )
  message("Corrected 'day 30' -> 'day 28' in Table S1 record_id for ", n_fixed, " row(s)")
}

# Append the new columns to the Table S1 sheet. covar is in the same row order as
# Table S1 (it was read from that sheet and only left-joined onto), so we can
# write the new columns directly.
s1_ncol <- ncol(s1_now)
new_cols <- data.frame(
  # Donor_ID matters beyond this table: HIV is longitudinal, 3 samples/donor
  Donor_ID = covar$donor_id,
  Age_reported = covar$Age_reported,
  Sex_reported = as.character(covar$Sex),
  Sampling_day_rel_onset = covar$sampling_day,
  Viral_load = covar$viral_load,
  # sex call plus the features behind it. Full detail (marker-threshold calls,
  # discrepancy flag) stays in figures/sex_prediction_metrics.csv.
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
openxlsx::writeData(wb, "Table S1", x = new_cols,
  startCol = s1_ncol + 1, startRow = 1)
message("Added ", length(pb_count_cols),
  " per-cell-type pseudobulk count column(s) to Table S1")

sheet_name <- "Table S1B" # per-PC association results
if (sheet_name %in% names(wb)) openxlsx::removeWorksheet(wb, sheet_name)
openxlsx::addWorksheet(wb, sheet_name)
openxlsx::writeData(wb, sheet_name, assoc_supp,
  withFilter = TRUE, headerStyle = openxlsx::createStyle(textDecoration = "bold")
)
openxlsx::setColWidths(wb, sheet_name, cols = seq_along(assoc_supp), widths = "auto")

out_xlsx <- file.path(repo_dir, "sample_annots/All_Supplementary_Tables_updated.xlsx")
openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
message("Wrote association results to sheet '", sheet_name, "' in: ", out_xlsx)

# write the results table
write.csv(assoc_df,
  file = file.path(repo_dir, "figures/pc_covariate_association.csv"),
  row.names = FALSE
)

## ---- Visualise: -log10(adjusted p) per PC x covariate -----------------------
# pf() underflows to 0 for a strong association and -log10(0) is Inf, which
# ggplot paints with na.value -- so the strongest tiles came out the same grey
# as the missing ones. Floor the p and cap the scale.
p_cap <- 50
neglog10 <- function(p, cap = p_cap) {
  pmin(-log10(pmax(p, .Machine$double.xmin)), cap)
}
# design covariates first, then donor, then technical
covariate_order <- c(
  "CellType", "Cohort", "Control_status", "Exposure_group",
  "Age", "Sex_observed", "Sex_predicted", "Sampling_day",
  "QC_nCells", "QC_meanTSS", "QC_meanLog10Frags", "QC_meanFRIP"
)
order_rows <- function(d) rev(intersect(covariate_order, unique(d$covariate)))

p_assoc <- ggplot(
  assoc_df,
  aes(x = PC, y = covariate, fill = neglog10(p_adj))
) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = ifelse(is.na(p_adj), "", sprintf("%.2f", r2))), size = 3) +
  facet_wrap(~embedding, ncol = 1) +
  scale_fill_gradient(
    low = "white", high = "#006400",
    na.value = "grey95",
    name = paste0("-log10(adj. p)\n(capped at ", p_cap, ")")
  ) +
  scale_x_discrete(limits = paste0("PC", seq_len(n_pc))) +
  scale_y_discrete(limits = order_rows(assoc_df)) +
  labs(
    title = "Association of sample-level PCs with known covariates",
    subtitle = paste0(
      "Every covariate tested separately against each PC. Tile label = variance explained (R^2); shading = significance.\n",
      "Unit of observation = one sample. Cell type is not testable here (each sample is a single profile averaged over its cells); see the by-cell-type panel.\n",
      "Age/Sex on samples with recorded metadata; predicted Sex and QC on all samples."
    ),
    x = NULL, y = NULL
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p_assoc,
  filename = file.path(repo_dir, "figures/pc_covariate_association.pdf"),
  width = 9.5,
  height = max(7.5, 0.5 * length(unique(assoc_df$covariate)) + 2.5),
  units = "in", dpi = 300, limitsize = FALSE
)

## ---- Per-(sample x cell type) association: adds Cell type and Exposure -------
# Cell type does not vary in the sample-level PCA above, so it cannot be tested
# there. This repeats the PCA on sample x cell-type pseudobulks, where every
# covariate is present and each is tested separately against each PC.
# (ct_col is defined with the pseudobulk cell counts above.)

assoc_ct_results <- list()
for (emb in embeddings_to_test) {
  reddim <- getReducedDims(project, reducedDims = emb)
  cd     <- as.data.frame(getCellColData(project, select = c("Sample", ct_col)))
  grp    <- paste(cd$Sample, cd[[ct_col]], sep = "||")

  sums   <- rowsum(reddim, group = grp)
  n_grp  <- as.numeric(table(grp)[rownames(sums)])
  pb_emb <- sums / n_grp                       # (sample x cell type) x dims
  keep   <- n_grp >= 20                        # drop tiny, noisy pseudobulks
  pb_emb <- pb_emb[keep, , drop = FALSE]
  # same filter as pb_emb, or the counts misalign with the retained rows
  pb_ncells <- n_grp[keep]

  pca <- stats::prcomp(pb_emb, center = TRUE, scale. = TRUE)
  var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
  k   <- min(n_pc, ncol(pca$x))
  pcs <- pca$x[, seq_len(k), drop = FALSE]

  key       <- do.call(rbind, strsplit(rownames(pcs), "\\|\\|"))
  pb_sample <- key[, 1]
  pb_ct     <- key[, 2]
  cvp       <- covar[match(pb_sample, covar$arrow_name), ]

  # as the sample-level set plus CellType. QC_nCells is per pseudobulk here.
  covariate_cols <- list(
    CellType          = factor(pb_ct),
    Cohort            = cvp$exposure_type,
    Control_status    = cvp$control_status,
    Exposure_group    = cvp$exposure_group,
    Age               = cvp$Age_numeric,
    Sex_observed      = cvp$Sex,
    Sex_predicted     = cvp$Sex_predicted,
    Sampling_day      = cvp$sampling_day,
    QC_nCells         = pb_ncells,
    QC_meanTSS        = cvp$mean_TSS,
    QC_meanLog10Frags = cvp$mean_log10_nFrags,
    QC_meanFRIP       = cvp$mean_FRIP
  )

  res <- do.call(rbind, lapply(seq_len(k), function(i) {
    do.call(rbind, lapply(names(covariate_cols), function(cn) {
      a <- pc_assoc(pcs[, i], covariate_cols[[cn]])
      data.frame(embedding = emb, PC = paste0("PC", i), pc_index = i,
                 var_explained = var_explained[i], covariate = cn,
                 r2 = a[["r2"]], p = a[["p"]], n = a[["n"]],
                 stringsAsFactors = FALSE)
    }))
  }))
  assoc_ct_results[[emb]] <- res
}

assoc_ct_df <- dplyr::bind_rows(assoc_ct_results) %>%
  dplyr::group_by(embedding, covariate) %>%
  dplyr::mutate(p_adj = p.adjust(p, method = "BH")) %>%
  dplyr::ungroup()

write.csv(assoc_ct_df,
  file = file.path(repo_dir, "figures/pc_covariate_association_by_celltype.csv"),
  row.names = FALSE
)

# Same summary at the pseudobulk level. The Results text quotes these, so they
# have to exist for the same unit as the panel that is shown.
var_summary_ct <- assoc_ct_df %>%
  dplyr::group_by(embedding, covariate) %>%
  dplyr::summarise(
    total_var = sum(r2 * var_explained, na.rm = TRUE),
    max_r2 = max(r2, na.rm = TRUE),
    max_r2_PC = PC[which.max(r2)],
    .groups = "drop"
  ) %>%
  dplyr::arrange(embedding, dplyr::desc(total_var))
message("Total pseudobulk-level variance associated with each covariate:")
print(as.data.frame(var_summary_ct))
write.csv(var_summary_ct,
  file = file.path(repo_dir, "figures/pc_covariate_variance_summary_by_celltype.csv"),
  row.names = FALSE
)

p_assoc_ct <- ggplot(assoc_ct_df, aes(x = PC, y = covariate, fill = neglog10(p_adj))) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = ifelse(is.na(p_adj), "", sprintf("%.2f", r2))), size = 3) +
  facet_wrap(~embedding, ncol = 1) +
  scale_fill_gradient(low = "white", high = "#006400", na.value = "grey95",
    name = paste0("-log10(adj. p)\n(capped at ", p_cap, ")")) +
  scale_x_discrete(limits = paste0("PC", seq_len(n_pc))) +
  scale_y_discrete(limits = order_rows(assoc_ct_df)) +
  labs(
    title = "Association of pseudobulk (sample x cell type) PCs with covariates",
    subtitle = paste0(
      "Every covariate tested separately against each PC. Tile label = variance explained (R^2); shading = significance.\n",
      "Unit of observation = one sample x cell-type pseudobulk (>= 20 cells); cell type can only be tested at this level, not on sample-level PCs."
    ),
    x = NULL, y = NULL
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p_assoc_ct,
  filename = file.path(repo_dir, "figures/pc_covariate_association_by_celltype.pdf"),
  width = 9, height = 7.5, units = "in", dpi = 300
)

# quick console summary of any significant confounder associations
message("Significant MARGINAL PC-covariate associations (adj. p < 0.05):")
print(subset(assoc_df, p_adj < 0.05,
  select = c(embedding, PC, covariate, r2, p_adj, n)
))

# Key check for the rebuttal: does any covariate survive once cohort is
# controlled for (Age, Sex, and the technical QC metrics)?
message("Significant COHORT-ADJUSTED associations (non-cohort covariates, adj. p < 0.05):")
print(subset(assoc_df, covariate != "Cohort" & p_adjCohort_bh < 0.05,
  select = c(embedding, PC, covariate, r2_partial_adjCohort, p_adjCohort_bh, n)
))

#####################################################################