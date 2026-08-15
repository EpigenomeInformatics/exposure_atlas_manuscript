#!/usr/bin/env Rscript

#####################################################################
# 07_6_covariate_balance.R
# created on 2026-08-12 by Irem B. Gunduz
# Are the technical covariates balanced between the compared groups?
#
# Different question from 07_3/07_4:
#   07_3/07_4  if we adjust, do the DAR calls change?   (sensitivity)
#   here       is there anything to adjust FOR?         (confounding)
#
# Adjustment only removes confounding if the covariate differs between the two
# groups being compared. TSS enrichment differs a lot BETWEEN cohorts (Supp Fig
# 1A) but every comparison is WITHIN a cohort, so that variation never enters
# the model.
#
# Reading the output:
#  - ignore the rank-sum p. At 4-7 per group it has no power, and n = 4 vs 4
#    cannot go below 0.029.
#  - |SMD| < 0.25 does NOT apply here. That convention is for large matched
#    studies. Under the null SMD has SE ~ sqrt(2/n), so median |SMD| lands near
#    0.43 at n = 5 and 0.65 at n = 3. An earlier version used 0.25 and called
#    203/300 tests imbalanced, which was the threshold, not the data.
#  - use p_perm / p_perm_BH: permute group labels within the comparison and see
#    where the observed |SMD| falls. Exact at these group sizes.
#  - permutation is per comparison. A covariate shifted the same way across many
#    comparisons sharing a control arm is one fact repeated, so read the
#    direction summary at the end too.
#####################################################################

suppressPackageStartupMessages({
  library(ArchR)
  library(ChrAccR)
  library(dplyr)
  library(readxl)
  library(ggplot2)
})
set.seed(12)

repo_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript"
fig_dir  <- file.path(repo_dir, "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

covid_dir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/"
other_dir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_2023-10-02/"
archr_dir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/"

grp_col <- "sample_exposure_group"

# groups are 3-8 samples, so the label space is small; this is plenty
n_perm <- 5000

# only drawn on the figure for reference; nothing is classified with it
smd_cut_convention <- 0.25

# as in 07_3, so the ChrAccR -> ArchR bridge behaves the same
norm_key <- function(x) {
  x <- basename(as.character(x))
  x <- sub("\\.tsv\\.gz.*$", "", x)
  x <- sub("_fragments$", "", x)
  x
}

## ---- 1. Per-sample covariates ----------------------------------------------
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
  )

archr_sample_meta <- as.data.frame(getCellColData(project))
archr_sample_meta <- archr_sample_meta[!duplicated(archr_sample_meta$Sample), , drop = FALSE]
id_lookup <- do.call(rbind, lapply(colnames(archr_sample_meta), function(cn) {
  data.frame(key = norm_key(archr_sample_meta[[cn]]),
    Sample = archr_sample_meta$Sample, stringsAsFactors = FALSE)
}))
id_lookup <- id_lookup[!is.na(id_lookup$key) & id_lookup$key != "", ]
id_lookup <- id_lookup[!duplicated(id_lookup$key), ]
rm(project); gc()

# donor and predicted sex from Table S1 (written by 01_v2_quality_control.R)
suppTables <- file.path(repo_dir, "sample_annots/All_Supplementary_Tables_updated.xlsx")
if (!file.exists(suppTables)) {
  suppTables <- file.path(repo_dir, "sample_annots/All_Supplementary_Tables.xlsx")
}
s1 <- readxl::read_excel(suppTables, sheet = "Table S1")
donor_sex <- data.frame(
  Sample = as.character(s1$arrow_name),
  Donor_ID = if ("Donor_ID" %in% names(s1)) as.character(s1$Donor_ID) else NA_character_,
  Sex_predicted = if ("Sex_predicted" %in% names(s1)) as.character(s1$Sex_predicted) else NA_character_,
  stringsAsFactors = FALSE
)

continuous_covs <- c("mean_TSS", "mean_FRIP", "mean_log10_nFrags", "n_cells")
categorical_covs <- c("Sex_predicted", "Donor_ID")

## ---- 2. Balance for one comparison -----------------------------------------
# pooled SD. Both groups constant -> 0 if means match, Inf otherwise.
smd <- function(x1, x2) {
  s1 <- stats::sd(x1, na.rm = TRUE); s2 <- stats::sd(x2, na.rm = TRUE)
  pooled <- sqrt((s1^2 + s2^2) / 2)
  d <- mean(x1, na.rm = TRUE) - mean(x2, na.rm = TRUE)
  if (!is.finite(pooled) || pooled == 0) return(ifelse(d == 0, 0, Inf))
  d / pooled
}

balance_one <- function(sa, cell, grp1, grp2) {
  keep <- as.character(sa[[grp_col]]) %in% c(grp1, grp2)
  sub <- sa[keep, , drop = FALSE]
  g <- as.character(sub[[grp_col]])
  i1 <- g == grp1; i2 <- g == grp2
  if (sum(i1) < 2 || sum(i2) < 2) {
    message("  too few samples (", sum(i1), " vs ", sum(i2), "), skipping")
    return(NULL)
  }

  # p floor for these group sizes
  p_floor <- tryCatch(
    stats::wilcox.test(seq_len(sum(i1)), sum(i1) + seq_len(sum(i2)))$p.value,
    error = function(e) NA_real_
  )

  cont <- do.call(rbind, lapply(continuous_covs, function(cc) {
    x1 <- as.numeric(sub[[cc]][i1]); x2 <- as.numeric(sub[[cc]][i2])
    p <- tryCatch(stats::wilcox.test(x1, x2, exact = FALSE)$p.value,
      error = function(e) NA_real_)
    obs <- smd(x1, x2)

    # null: same values, same group sizes, shuffled labels
    xall <- c(x1, x2)
    n1 <- length(x1)
    null_smd <- replicate(n_perm, {
      idx <- sample.int(length(xall))
      smd(xall[idx[seq_len(n1)]], xall[idx[-seq_len(n1)]])
    })
    null_smd <- null_smd[is.finite(null_smd)]
    # +1 so p is never exactly 0
    p_perm <- if (length(null_smd) == 0 || !is.finite(obs)) NA_real_ else
      (1 + sum(abs(null_smd) >= abs(obs))) / (1 + length(null_smd))

    data.frame(
      cell = cell, comparison = paste(grp1, "vs", grp2), covariate = cc,
      type = "continuous",
      n_grp1 = sum(i1), n_grp2 = sum(i2),
      mean_grp1 = mean(x1, na.rm = TRUE), mean_grp2 = mean(x2, na.rm = TRUE),
      sd_grp1 = stats::sd(x1, na.rm = TRUE), sd_grp2 = stats::sd(x2, na.rm = TRUE),
      SMD = obs,
      # what balance looks like at these group sizes
      null_median_abs_SMD = if (length(null_smd)) stats::median(abs(null_smd)) else NA_real_,
      null_q95_abs_SMD = if (length(null_smd)) as.numeric(stats::quantile(abs(null_smd), 0.95)) else NA_real_,
      p_perm = p_perm,
      p_wilcoxon = p, p_floor = p_floor,
      stringsAsFactors = FALSE
    )
  }))

  cat_rows <- lapply(categorical_covs, function(cc) {
    v <- as.character(sub[[cc]])
    if (all(is.na(v)) || dplyr::n_distinct(v[!is.na(v)]) < 2) return(NULL)
    tb <- table(v, g)
    p <- tryCatch(stats::fisher.test(tb, simulate.p.value = TRUE, B = 20000)$p.value,
      error = function(e) NA_real_)
    data.frame(
      cell = cell, comparison = paste(grp1, "vs", grp2), covariate = cc,
      type = "categorical",
      n_grp1 = sum(i1), n_grp2 = sum(i2),
      mean_grp1 = NA_real_, mean_grp2 = NA_real_,
      sd_grp1 = NA_real_, sd_grp2 = NA_real_,
      SMD = NA_real_,
      null_median_abs_SMD = NA_real_, null_q95_abs_SMD = NA_real_,
      p_perm = p, p_wilcoxon = p, p_floor = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(cont, dplyr::bind_rows(cat_rows))
}

## ---- 3. Sweep every cell type and comparison -------------------------------
# One loadDsAcc per CELL TYPE (not per comparison): the sample annotation is
# shared across the comparisons of a cell type, so loading it once is enough.
collect_balance <- function(anaDir) {
  if (!dir.exists(anaDir)) {
    message("analysis directory not found: ", anaDir)
    return(NULL)
  }
  cells <- basename(list.dirs(anaDir, recursive = FALSE))
  out <- list()
  for (cell in cells) {
    ds_path <- file.path(anaDir, cell, "data", "dsATAC_filtered")
    ct_file <- file.path(anaDir, cell, "reports", "differential_data",
      "comparisonTable.rds")
    if (!dir.exists(ds_path) || !file.exists(ct_file)) next

    sa <- tryCatch({
      ds <- ChrAccR::loadDsAcc(ds_path)
      a <- ChrAccR::getSampleAnnot(ds)
      rm(ds); gc()
      a
    }, error = function(e) {
      message("  ", cell, ": could not load DsATAC (", conditionMessage(e), ")")
      NULL
    })
    if (is.null(sa)) next
    if (!grp_col %in% colnames(sa)) {
      message("  ", cell, ": no '", grp_col, "' column, skipping")
      next
    }

    # bridge to the ArchR QC metrics exactly as 07_3 does
    try_map <- function(v) id_lookup$Sample[match(norm_key(v), id_lookup$key)]
    cand <- c(list(.rownames = rownames(sa)), as.list(as.data.frame(sa)))
    mapped <- lapply(cand, try_map)
    best <- names(cand)[which.max(vapply(mapped, function(m) sum(!is.na(m)), integer(1)))]
    arch_samp <- mapped[[best]]
    idx <- match(arch_samp, qc_by_sample$Sample)
    if (mean(is.na(idx)) > 0) {
      message("  ", cell, ": ", sum(is.na(idx)), "/", length(idx),
        " samples unmatched to QC metrics, skipping")
      next
    }
    for (cc in continuous_covs) sa[[cc]] <- qc_by_sample[[cc]][idx]
    ds_i <- match(arch_samp, donor_sex$Sample)
    sa[["Donor_ID"]] <- donor_sex$Donor_ID[ds_i]
    sa[["Sex_predicted"]] <- donor_sex$Sex_predicted[ds_i]

    ct <- readRDS(ct_file)
    message("=== ", cell, ": ", nrow(ct), " comparison(s)")
    for (i in seq_len(nrow(ct))) {
      message("  ", ct$grp1[i], " vs ", ct$grp2[i])
      out[[length(out) + 1L]] <- balance_one(sa, cell, ct$grp1[i], ct$grp2[i])
    }
  }
  dplyr::bind_rows(out)
}

bal <- dplyr::bind_rows(collect_balance(covid_dir), collect_balance(other_dir))
if (is.null(bal) || nrow(bal) == 0) {
  stop("No comparison could be assessed. Check the ChrAccR run directories.")
}

bal <- bal %>%
  dplyr::mutate(
    # BH over all continuous tests (conservative; per-covariate also defensible)
    p_perm_BH = ifelse(type == "continuous",
      stats::p.adjust(ifelse(type == "continuous", p_perm, NA_real_), method = "BH"),
      NA_real_),
    imbalanced = ifelse(type == "continuous" & !is.na(p_perm_BH) & p_perm_BH < 0.05,
      "imbalanced", "not distinguishable from chance"),
    # observed effect relative to chance at this n
    SMD_vs_null = ifelse(is.finite(SMD) & !is.na(null_median_abs_SMD) &
        null_median_abs_SMD > 0, abs(SMD) / null_median_abs_SMD, NA_real_),
    p_uninformative = !is.na(p_floor) & p_floor > 0.05
  )

write.csv(bal %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 4))),
  file.path(fig_dir, "covariate_balance_by_comparison.csv"), row.names = FALSE)

## ---- 4. Headline summary ----------------------------------------------------
cont <- bal[bal$type == "continuous" & is.finite(bal$SMD), ]
n_imb <- sum(cont$imbalanced == "imbalanced", na.rm = TRUE)
message("\n", nrow(cont), " covariate x comparison tests (continuous covariates)")
message(n_imb, " imbalanced beyond chance (permutation p, BH < 0.05); ",
  nrow(cont) - n_imb, " not distinguishable from chance")
message("For reference, ", sum(abs(cont$SMD) >= smd_cut_convention, na.rm = TRUE),
  " exceed the |SMD| > ", smd_cut_convention, " convention -- which is NOT a ",
  "meaningful count at these group sizes (see header note 2)")
if (any(cont$p_uninformative)) {
  message(sum(cont$p_uninformative), " test(s) have a rank-sum p-value floor above 0.05")
}

message("\nStrongest imbalances (by permutation p, then effect relative to null):")
print(as.data.frame(
  cont %>%
    dplyr::arrange(p_perm, dplyr::desc(SMD_vs_null)) %>%
    dplyr::transmute(cell, comparison, covariate, n_grp1, n_grp2,
      mean_grp1 = signif(mean_grp1, 4), mean_grp2 = signif(mean_grp2, 4),
      SMD = round(SMD, 3),
      null_median = round(null_median_abs_SMD, 3),
      SMD_vs_null = round(SMD_vs_null, 2),
      p_perm = signif(p_perm, 3), p_perm_BH = signif(p_perm_BH, 3)) %>%
    utils::head(20)
))

per_cov <- cont %>%
  dplyr::group_by(covariate) %>%
  dplyr::summarise(
    n_comparisons     = dplyr::n(),
    median_abs_SMD    = round(stats::median(abs(SMD)), 3),
    null_median_abs_SMD = round(stats::median(null_median_abs_SMD, na.rm = TRUE), 3),
    max_abs_SMD       = round(max(abs(SMD)), 3),
    n_imbalanced_perm = sum(imbalanced == "imbalanced", na.rm = TRUE),
    min_p_perm        = signif(min(p_perm, na.rm = TRUE), 3),
    # consistent direction across comparisons, even if none is significant
    # on its own. This is what caught the COVID-19 cell-number difference.
    pct_higher_in_grp1 = round(100 * mean(SMD > 0), 1),
    p_sign_test       = signif(stats::binom.test(sum(SMD > 0), dplyr::n())$p.value, 3),
    .groups = "drop"
  )
message("\nPer covariate (median_abs_SMD is only interpretable against null_median_abs_SMD):")
print(as.data.frame(per_cov))
write.csv(per_cov, file.path(fig_dir, "covariate_balance_summary.csv"),
  row.names = FALSE)

## ---- 4b. Direction consistency within a cohort ------------------------------
# Comparisons sharing a control arm are not independent, so this is
# descriptive. It exists because per-comparison tests cannot see a consistent
# shift across comparisons.
cohort_of <- function(x) sub("_.*", "", x)
direction <- cont %>%
  dplyr::mutate(cohort = cohort_of(comparison)) %>%
  dplyr::group_by(cohort, covariate) %>%
  dplyr::summarise(
    n_comparisons = dplyr::n(),
    median_SMD = round(stats::median(SMD), 3),
    pct_same_direction = round(100 * max(mean(SMD > 0), mean(SMD < 0)), 1),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(abs(median_SMD)))
message("\nDirection consistency by cohort (descriptive; comparisons sharing a control arm are not independent):")
print(as.data.frame(direction))
write.csv(direction, file.path(fig_dir, "covariate_balance_direction.csv"),
  row.names = FALSE)

## ---- 5. Figure --------------------------------------------------------------
# Signed SMD, cohort on the y axis: is a covariate pushed to one side in a
# given cohort? Shaded band is the chance envelope from the permutations.
# (An earlier version plotted |SMD| against the null and was unreadable --
# nothing was significant, so the colour carried nothing, and absolute values
# hid the direction, which is the actual finding.)
# canonical cohort colours (same values as comp_group_colors in 14_l2fc.R)
cohort_colors <- c(
  "C19" = "#238B45", "HIV" = "#88419D",
  "Influenza" = "#D95F0E", "OP" = "#084594"
)

plot_bal <- cont %>%
  dplyr::mutate(cohort = cohort_of(comparison)) %>%
  dplyr::filter(is.finite(SMD))
plot_bal$cohort <- factor(plot_bal$cohort,
  levels = rev(intersect(names(cohort_colors), unique(plot_bal$cohort))))

# one band per covariate; null median varies with group size, so take the median
band <- plot_bal %>%
  dplyr::group_by(covariate) %>%
  dplyr::summarise(hi = stats::median(null_median_abs_SMD, na.rm = TRUE),
    .groups = "drop") %>%
  dplyr::mutate(lo = -hi)

# covariate x cohort combinations outside the band, i.e. worth a sentence
callout <- plot_bal %>%
  dplyr::group_by(covariate, cohort) %>%
  dplyr::summarise(median_SMD = stats::median(SMD), n = dplyr::n(), .groups = "drop") %>%
  dplyr::left_join(band, by = "covariate") %>%
  dplyr::filter(abs(median_SMD) > hi)

p_bal <- ggplot(plot_bal, aes(x = SMD, y = cohort)) +
  geom_rect(data = band, inherit.aes = FALSE,
    aes(xmin = lo, xmax = hi, ymin = -Inf, ymax = Inf),
    fill = "grey88", alpha = 0.7) +
  geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.3) +
  geom_point(aes(colour = cohort), size = 1.9, alpha = 0.75,
    position = position_jitter(height = 0.18, width = 0)) +
  # cohort median, so a consistent shift shows even when points scatter
  stat_summary(fun = median, geom = "point", shape = 124, size = 5,
    colour = "grey15") +
  facet_wrap(~covariate) +
  scale_colour_manual(values = cohort_colors, guide = "none") +
  labs(
    title = "Technical covariate balance within each differential comparison",
    subtitle = paste0(
      "Each point is one comparison; positive values mean the covariate is higher in the first group of the contrast.\n",
      "Grey band: the range that label permutation alone produces at these group sizes (n = 3-8, ", n_perm,
      " permutations). Vertical tick: cohort median.\n",
      "Points inside the band are indistinguishable from chance. A cohort whose points sit consistently to one side is systematically shifted."),
    x = "standardised mean difference between the compared groups", y = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(plot.subtitle = element_text(size = 8),
    strip.background = element_blank(), strip.text = element_text(face = "bold"),
    panel.grid.major.y = element_line(colour = "grey94", linewidth = 0.3))

ggsave(file.path(fig_dir, "covariate_balance.pdf"), p_bal,
  width = 9.5, height = 5.5)

if (nrow(callout) > 0) {
  message("\nCovariate x cohort combinations whose median sits outside the chance band:")
  print(as.data.frame(callout %>%
    dplyr::transmute(covariate, cohort, n_comparisons = n,
      median_SMD = round(median_SMD, 3), chance_band = round(hi, 3))))
} else {
  message("\nNo covariate x cohort combination sits outside the chance band.")
}

message("\nDone. Tables and figure in ", fig_dir)
message("Read alongside confounder_adjusted_DARs.pdf: this script asks whether ",
  "there is anything to adjust for, that one asks whether adjusting changes the answer.")
message("Report p_perm_BH, not the |SMD| convention. If a covariate is flagged ",
  "imbalanced, or shows a consistent direction in covariate_balance_direction.csv, ",
  "name it in the limitations rather than claiming balance.")

#####################################################################
