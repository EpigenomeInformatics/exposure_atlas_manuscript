#!/usr/bin/env Rscript

#####################################################################
# 07_7_covariate_balance_response_figure.R
# created on 2026-08-18 by Irem B. Gunduz
# Response figure: covariate balance, and the stability of the effect sizes
# and TF enrichments under QC adjustment.
#
# in:  figures/covariate_balance_by_comparison.csv   (07_6_covariate_balance.R)
#      figures/confounder_adjusted_DAR_summary.csv   (07_4_confounder_adjusted_DAR_plots.R)
#      figures/confounder_adjusted_lola_summary.csv  (07_5_confounder_adjusted_lola.R)
#      figures/confounder_adjusted_lola_terms.csv    (07_5_confounder_adjusted_lola.R)
#      <out_dir>/confounder_adjusted_DAR_regions.csv (07_4_confounder_adjusted_DAR_plots.R)
# out: figures/response_figure_2_covariate_balance.pdf
#
# Each |SMD| is scaled by its own permutation null 95th percentile, so 1 is a
# significance boundary rather than a fixed convention.
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})
set.seed(3)

repo_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript"
fig_dir  <- file.path(repo_dir, "figures")
out_dir  <- "/icbb/projects/igunduz/finalize_echo_050824/confounder_adjusted/"

bal_csv      <- file.path(fig_dir, "covariate_balance_by_comparison.csv")
dar_sum_csv  <- file.path(fig_dir, "confounder_adjusted_DAR_summary.csv")
dar_reg_csv  <- file.path(out_dir, "confounder_adjusted_DAR_regions.csv")
lola_sum_csv <- file.path(fig_dir, "confounder_adjusted_lola_summary.csv")
lola_trm_csv <- file.path(fig_dir, "confounder_adjusted_lola_terms.csv")
out_pdf      <- file.path(fig_dir, "response_figure_2_covariate_balance.pdf")

stopifnot(file.exists(bal_csv), file.exists(dar_sum_csv), file.exists(dar_reg_csv),
          file.exists(lola_sum_csv), file.exists(lola_trm_csv))

focus_set  <- "TSS_FRIP"
focus_cell <- "Mono_CD14"
n_label    <- 3
n_tf_label <- 8
tf_colls   <- c("TF_motif_clusters", "TF_motifs", "codex", "encode_tfbs")

cohort_cols <- c("COVID-19" = "#3b6ea5", "HIV" = "#c1741f",
                 "Influenza" = "#4f8a5b", "OP" = "#7a5c9e")
status_cols <- c("DAR in both" = "#4A2377", "DAR unadjusted only" = "#3C5488",
                 "DAR adjusted only" = "#E64B35")

theme_resp <- theme_classic(base_size = 11) +
  theme(panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
        panel.grid.minor = element_line(colour = "grey92", linewidth = 0.2),
        plot.title = element_text(face = "bold", hjust = 0, size = 13),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

## ---- Panel A: |SMD| relative to the null 95th percentile -------------------
bal <- read.csv(bal_csv, stringsAsFactors = FALSE) %>%
  dplyr::filter(!p_uninformative, !is.na(SMD)) %>%
  dplyr::mutate(
    ratio  = abs(SMD) / null_q95_abs_SMD,
    label  = paste0(cell, ", ", comparison),
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

top_bal <- bal %>%
  dplyr::arrange(dplyr::desc(ratio)) %>%
  dplyr::filter(ratio > 1 | dplyr::row_number() <= n_label)

jit <- position_jitter(width = 0.24, height = 0, seed = 3)

pA <- ggplot(bal, aes(x = covariate, y = ratio, colour = cohort)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0, ymax = 1,
           fill = "#eaf2ea", colour = NA) +
  geom_point(position = jit, size = 1.9, alpha = 0.8) +
  geom_hline(yintercept = 1, colour = "#c0392b", linewidth = 0.7) +
  geom_text_repel(data = top_bal, aes(label = label), position = jit,
                  size = 2.7, colour = "grey20", segment.colour = "grey60",
                  segment.size = 0.3, min.segment.length = 0,
                  box.padding = 0.45, max.overlaps = Inf,
                  show.legend = FALSE) +
  annotate("text", x = 0.6, y = 1.02, label = "chance ceiling",
           colour = "#c0392b", fontface = "italic", size = 3.3, hjust = 0, vjust = 0) +
  scale_colour_manual(values = cohort_cols, name = NULL) +
  coord_cartesian(ylim = c(0, 1.45)) +
  labs(x = NULL,
       y = "Observed difference between groups\n(relative to its own chance ceiling)",
       title = "A") +
  theme_resp +
  theme(legend.position = "bottom")

if (n_above > 0) {
  pA <- pA + annotate("text", x = 1.1, y = 1.38, size = 3.1, colour = "#444444",
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
  theme_resp

## ---- Panel C: fold change vs fold change, COVID-19 CD14+ monocytes ---------
reg <- read.csv(dar_reg_csv, stringsAsFactors = FALSE)
if (!"adj_set" %in% names(reg)) {
  stop("confounder_adjusted_DAR_regions.csv has no adj_set column; re-run ",
       "07_4_confounder_adjusted_DAR_plots.R to write it.")
}

reg <- reg %>%
  dplyr::filter(cell == focus_cell, adj_set == focus_set, grepl("^C19", comparison)) %>%
  dplyr::mutate(status = factor(status, levels = names(status_cols))) %>%
  dplyr::arrange(status)

dar_r <- read.csv(dar_sum_csv, stringsAsFactors = FALSE) %>%
  dplyr::filter(cell == focus_cell, adj_set == focus_set,
                comparison %in% unique(reg$comparison)) %>%
  dplyr::mutate(lab = sprintf("r = %.3f\n%d of %d recovered",
                              lfc_pearson_all, shared, DARs_unadjusted))

message(sprintf("panel C: %d regions in %d %s comparisons",
                nrow(reg), dplyr::n_distinct(reg$comparison), focus_cell))

pC <- ggplot(reg, aes(x = log2FC_unadj, y = log2FC_adj, colour = status)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 0.7, alpha = 0.6) +
  geom_text(data = dar_r, inherit.aes = FALSE, aes(x = -Inf, y = Inf, label = lab),
            hjust = -0.08, vjust = 1.2, size = 3, colour = "grey20") +
  facet_wrap(~comparison, scales = "free") +
  scale_colour_manual(values = status_cols, name = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(x = "log2 fold change, unadjusted model",
       y = "log2 fold change, QC-adjusted model", title = "C") +
  theme_resp +
  theme(legend.position = "bottom")

## ---- Panel D: TF enrichment, adjusted vs unadjusted ------------------------
lola <- read.csv(lola_trm_csv, stringsAsFactors = FALSE) %>%
  dplyr::filter(collection %in% tf_colls,
                is.finite(log2OR_unadj), is.finite(log2OR_adj)) %>%
  dplyr::mutate(panel = paste0(cell, ", ", comparison),
                motif = sub("^MA[0-9.]+_", "", description))

top_tf <- lola %>%
  dplyr::group_by(panel, motif) %>%
  dplyr::slice_max(abs(log2OR_unadj), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::slice_max(abs(log2OR_unadj) + abs(log2OR_adj), n = n_tf_label)

panel_cols <- stats::setNames(
  c("#3b6ea5", "#7fa6cd", "#c1741f", "#e0a15c", "#4f8a5b", "#7a5c9e")[
    seq_along(sort(unique(lola$panel)))],
  sort(unique(lola$panel)))

lola_r <- read.csv(lola_sum_csv, stringsAsFactors = FALSE) %>%
  dplyr::filter(adj_set == focus_set,
                paste0(cell, ", ", comparison) %in% unique(lola$panel)) %>%
  dplyr::mutate(txt = sprintf("%s, %s: r = %.2f", cell, comparison, log2OR_pearson))

message(sprintf("panel D: %d TF region sets in %d comparisons",
                nrow(lola), dplyr::n_distinct(lola$panel)))

pD <- ggplot(lola, aes(x = log2OR_unadj, y = log2OR_adj, colour = panel)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
  geom_point(size = 1.3, alpha = 0.75) +
  geom_text_repel(data = top_tf, aes(label = motif), size = 2.7,
                  colour = "grey20", segment.colour = "grey60",
                  segment.size = 0.3, min.segment.length = 0,
                  box.padding = 0.4, max.overlaps = Inf,
                  show.legend = FALSE) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2, size = 3,
           colour = "grey20", label = paste(lola_r$txt, collapse = "\n")) +
  scale_colour_manual(values = panel_cols, name = NULL) +
  scale_x_continuous(expand = expansion(mult = 0.08)) +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 2.5, alpha = 1))) +
  labs(x = "log2 odds ratio, unadjusted DARs",
       y = "log2 odds ratio, QC-adjusted DARs", title = "D") +
  theme_resp +
  theme(legend.position = "bottom")

## ---- Compose ---------------------------------------------------------------
row1 <- (pA | pB) + patchwork::plot_layout(widths = c(1.75, 1))
row2 <- (pC | pD) + patchwork::plot_layout(widths = c(1.5, 1))
fig  <- row1 / row2

ggsave(out_pdf, fig, width = 13.5, height = 10.4, useDingbats = FALSE)
message("Wrote ", out_pdf)
#####################################################################
