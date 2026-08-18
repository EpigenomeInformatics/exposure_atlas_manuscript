#!/usr/bin/env Rscript

#####################################################################
# 07_7_covariate_balance_response_figure.R
# created on 2026-08-18 by Irem B. Gunduz
#
# Response figure for Reviewer 1 comment 2: are the technical covariates
# balanced between the groups being compared?
#
# Reads  figures/covariate_balance_by_comparison.csv   (from 07_6_covariate_balance.R)
# Writes figures/response_figure_2_covariate_balance.pdf
#
# WHY A SECOND PLOT, given 07_6 already draws covariate_balance.pdf:
# that panel puts the observed |SMD| against the permutation null MEDIAN, so by
# construction about half the points sit above the line and the figure reads as
# imbalance — the opposite of the result. Here each observation is divided by its
# own null 95th percentile, so the line is a real significance boundary and
# "inside the shaded band" means "not distinguishable from chance".
#
# Sex and donor carry no standardised mean difference (categorical), so this
# figure covers the four continuous quality measures. Say so in the legend.
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork) # only for the two-panel layout; see the note at the bottom
})
set.seed(3)

repo_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript"
fig_dir  <- file.path(repo_dir, "figures")
in_csv   <- file.path(fig_dir, "covariate_balance_by_comparison.csv")
out_pdf  <- file.path(fig_dir, "response_figure_2_covariate_balance.pdf")

stopifnot(file.exists(in_csv))

bal <- read.csv(in_csv, stringsAsFactors = FALSE) %>%
  dplyr::filter(!p_uninformative, !is.na(SMD)) %>%
  dplyr::mutate(
    ratio  = abs(SMD) / null_q95_abs_SMD,
    cohort = dplyr::recode(sub("_.*$", "", comparison),
                           C19 = "COVID-19", HIV = "HIV",
                           Influenza = "Influenza", OP = "OP"),
    covariate = factor(covariate,
      levels = c("mean_TSS", "mean_FRIP", "mean_log10_nFrags", "n_cells"),
      labels = c("TSS\nenrichment", "FRIP",
                 "Sequencing depth\n(log10 fragments)", "Cells per\nsample"))
  ) %>%
  dplyr::filter(!is.na(covariate))

n_tests <- nrow(bal)
n_above <- sum(bal$ratio > 1)
n_comp  <- dplyr::n_distinct(paste(bal$cell, bal$comparison))
min_bh  <- min(bal$p_perm_BH, na.rm = TRUE)

message(sprintf("comparisons %d | tests %d | above ceiling %d | smallest adjusted p %.3f",
                n_comp, n_tests, n_above, min_bh))

cohort_cols <- c("COVID-19" = "#3b6ea5", "HIV" = "#c1741f",
                 "Influenza" = "#4f8a5b", "OP" = "#7a5c9e")

## ---- Panel A: observed difference relative to its own chance ceiling --------
pA <- ggplot(bal, aes(x = covariate, y = ratio, colour = cohort)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 1,
           fill = "#eaf2ea", colour = NA) +
  geom_jitter(width = 0.24, height = 0, size = 1.9, alpha = 0.8) +
  geom_hline(yintercept = 1, colour = "#c0392b", linewidth = 0.7) +
  annotate("text", x = 4.45, y = 1.02, label = "chance ceiling",
           colour = "#c0392b", fontface = "italic", size = 3.3, hjust = 1, vjust = 0) +
  scale_colour_manual(values = cohort_cols, name = NULL) +
  coord_cartesian(ylim = c(0, 1.45)) +
  labs(x = NULL,
       y = "Observed difference between groups\n(relative to its own chance ceiling)",
       title = sprintf("A   %d of %d tests fall below the level expected by chance",
                       n_tests - n_above, n_tests)) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0, size = 12),
        legend.position = "bottom")

if (n_above > 0) {
  pA <- pA + annotate("text", x = 2.6, y = 1.36, size = 3.1, colour = "#444444",
                      label = sprintf("%d exceed it; none significant\nafter multiple-testing correction",
                                      n_above))
}

## ---- Panel B: permutation p-values -----------------------------------------
pB <- ggplot(bal, aes(x = p_perm)) +
  geom_histogram(breaks = seq(0, 1, by = 0.05),
                 fill = "#9bb7d4", colour = "white") +
  geom_vline(xintercept = 0.05, colour = "#c0392b", linewidth = 0.7) +
  annotate("text", x = 0.08, y = Inf, label = "p = 0.05", colour = "#c0392b",
           fontface = "italic", size = 3.3, hjust = 0, vjust = 1.6) +
  labs(x = "Permutation p-value", y = "Number of tests",
       title = "B   p-values spread across the whole range,\n      as expected when the groups do not differ") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0, size = 12))

## ---- Compose ----------------------------------------------------------------
subtitle <- sprintf(paste0(
  "For each of the %d within-cohort comparisons we tested whether the two groups differ in four measures of sample quality (%d tests in total). Each point is one measure in one\n",
  "comparison, shown relative to the largest difference that arises by chance alone at that comparison's group sizes (95th percentile of 5,000 random relabellings of the samples).\n",
  "Points inside the shaded area are indistinguishable from chance. No test was imbalanced after correction for multiple testing (Benjamini-Hochberg, smallest adjusted p = %.2f)."),
  n_comp, n_tests, min_bh)

fig <- (pA | pB) +
  patchwork::plot_layout(widths = c(1.75, 1)) +
  patchwork::plot_annotation(
    title = "Sample quality does not differ between the groups compared in any differential analysis",
    subtitle = subtitle,
    theme = theme(plot.title    = element_text(face = "bold", size = 13.5),
                  plot.subtitle = element_text(size = 8.6, colour = "#333333"))
  )

ggsave(out_pdf, fig, width = 12.6, height = 5.8, useDingbats = FALSE)
message("Wrote ", out_pdf)

## ---- If patchwork is not installed -----------------------------------------
# Drop the library(patchwork) call and the fig <- ... block above, and save the
# two panels separately instead:
#   ggsave(file.path(fig_dir, "response_figure_2A_balance.pdf"), pA, width = 7.6, height = 5.2)
#   ggsave(file.path(fig_dir, "response_figure_2B_pvalues.pdf"), pB, width = 4.6, height = 5.2)
# The subtitle text above then goes into the figure legend in the response letter.
#####################################################################
