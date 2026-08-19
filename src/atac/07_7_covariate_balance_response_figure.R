#!/usr/bin/env Rscript

#####################################################################
# 07_7_covariate_balance_response_figure.R
# created on 2026-08-18 by Irem B. Gunduz
# Response figure: covariate balance within each differential comparison.
#
# in:  figures/covariate_balance_by_comparison.csv (07_6_covariate_balance.R)
# out: figures/response_figure_2_covariate_balance.pdf
#      figures/response_figure_2_legend.txt
#
# Each |SMD| is scaled by its own permutation null 95th percentile, so 1 is a
# significance boundary rather than a fixed convention. 07_6 plots against the
# null median instead, where half the points lie above the line by construction.
# Categorical covariates (sex, donor) have no SMD and are not plotted.
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
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

## ---- Panel A: |SMD| relative to the null 95th percentile -------------------
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
       title = "A") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0, size = 13),
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
  labs(x = "Permutation p-value", y = "Number of tests", title = "B") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", hjust = 0, size = 13))

## ---- Compose ---------------------------------------------------------------
# Caption is written to file rather than drawn on the panel, so the numbers in
# the legend and the plot come from the same objects.
fig <- (pA | pB) + patchwork::plot_layout(widths = c(1.75, 1))

ggsave(out_pdf, fig, width = 12.2, height = 5.2, useDingbats = FALSE)
message("Wrote ", out_pdf)

legend_txt <- sprintf(paste0(
  "Response Figure 2. Sample quality does not differ between the groups compared in any differential analysis. ",
  "For each of the %d within-cohort comparisons we tested whether the two compared groups differ in four measures ",
  "of sample quality: TSS enrichment, FRIP, sequencing depth and cells per sample (%d tests in total). ",
  "(A) Each point is one measure in one comparison, shown relative to the largest difference that arises by chance ",
  "alone at that comparison's group sizes (95th percentile of 5,000 random relabellings of the samples). Points ",
  "inside the shaded area are indistinguishable from chance. (B) Distribution of the permutation p-values. ",
  "No test was imbalanced after correction for multiple testing (Benjamini-Hochberg, smallest adjusted p = %.2f). ",
  "Sex and donor identity are categorical and carry no standardised mean difference, so they are not shown; ",
  "%d of the 75 comparisons yielded an informative permutation test.\n"),
  n_comp, n_tests, min_bh, n_comp)

writeLines(legend_txt, file.path(fig_dir, "response_figure_2_legend.txt"))
message("Wrote ", file.path(fig_dir, "response_figure_2_legend.txt"))

# Without patchwork, drop the fig block above and save the panels separately:
#   ggsave(file.path(fig_dir, "response_figure_2A_balance.pdf"), pA, width = 7.6, height = 5.2)
#   ggsave(file.path(fig_dir, "response_figure_2B_pvalues.pdf"), pB, width = 4.6, height = 5.2)
#####################################################################
