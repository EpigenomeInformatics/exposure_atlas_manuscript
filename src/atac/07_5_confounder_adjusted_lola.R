#!/usr/bin/env Rscript

#####################################################################
# 07_5_confounder_adjusted_lola.R
# created on 2026-08-12 by Irem B. Gunduz
# LOLA enrichment on the covariate-adjusted DARs vs the unadjusted ones
#
# DAR counts surviving adjustment is not the claim; the enrichment is. Targets:
#   COVID-19  Mono_CD14  C19_sev vs C19_ctrl, C19_mod vs C19_ctrl
#   HIV       T_mem_CD8  HIV_ctrl vs HIV_acu, HIV_ctrl vs HIV_chr
#
# LOLA is recomputed here rather than read from ChrAccR: the adjusted runs in
# 07_3 had no lolaDbPaths set, and both sides need the same region cutoff and
# the same universe or a change in enrichment is confounded with a change in
# how the sets were built.
#
# Nothing is refitted; this reads the diffTabs 07_3 wrote.
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(GenomicRanges)
  library(LOLA)
})
set.seed(12)

repo_dir <- "/scratch/icbb/igunduz/irem_github/exposure_atlas_manuscript"
fig_dir  <- file.path(repo_dir, "figures")
out_dir  <- "/icbb/projects/igunduz/finalize_echo_050824/confounder_adjusted/"
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# unadjusted ChrAccR runs, as in 07_3 / 07_4
covid_dir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/"
other_dir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_2023-10-02/"

# which adjustment design to compare against the unadjusted fit
focus_set <- "TSS_FRIP"

# DAR definition, identical on both sides
l2fc_cut <- 0.5
padj_cut <- 0.05

# collections the manuscript discusses
collections_of_interest <- c("TF_motif_clusters", "TF_motifs", "codex", "encode_tfbs")

## ---- LOLA database ----------------------------------------------------------
lola_candidates <- c(
  "/icbb/projects/share/annotations/lolaDB/hg38/",
  "/icbb/projects/igunduz/annotation/lolaDB/hg38/"
)
lola_path <- lola_candidates[dir.exists(lola_candidates)][1]
if (is.na(lola_path)) {
  stop("No LOLA database found. Looked in:\n  ",
    paste(lola_candidates, collapse = "\n  "))
}
message("Loading LOLA database: ", lola_path)
lolaDb <- LOLA::loadRegionDB(lola_path)

## ---- Targets ----------------------------------------------------------------
targets <- list(
  list(dir = covid_dir, cell = "Mono_CD14", grp1 = "C19_sev", grp2 = "C19_ctrl"),
  list(dir = covid_dir, cell = "Mono_CD14", grp1 = "C19_mod", grp2 = "C19_ctrl"),
  list(dir = other_dir, cell = "T_mem_CD8", grp1 = "HIV_ctrl", grp2 = "HIV_acu"),
  list(dir = other_dir, cell = "T_mem_CD8", grp1 = "HIV_ctrl", grp2 = "HIV_chr")
)

## ---- diffTab lookup ---------------------------------------------------------
list_diff_tabs <- function(dir) {
  f <- list.files(dir, pattern = "diffTab.*archrPeaks.*\\.tsv$",
    recursive = TRUE, full.names = TRUE)
  if (length(f) == 0) {
    f <- list.files(dir, pattern = "diffTab.*\\.tsv$",
      recursive = TRUE, full.names = TRUE)
  }
  sort(f)
}

adjusted_dir_for <- function(cell, grp1, grp2, set_name) {
  tag <- paste0(grp1, "_vs_", grp2)
  cands <- c(
    file.path(out_dir, paste0(cell, "__", tag, "__adj-", set_name)),
    # legacy "__adjusted" dirs are the TSS + FRIP design
    if (identical(set_name, "TSS_FRIP")) {
      file.path(out_dir, paste0(cell, "__", tag, "__adjusted"))
    }
  )
  cands <- cands[dir.exists(cands)]
  if (length(cands) == 0) NA_character_ else cands[1]
}

unadjusted_tab_for <- function(anaDir, cell, grp1, grp2) {
  ddir <- file.path(anaDir, cell, "reports", "differential_data")
  if (!dir.exists(ddir)) return(NA_character_)
  files <- list_diff_tabs(ddir)
  if (length(files) == 0) return(NA_character_)
  hit <- files[grepl(grp1, basename(files), fixed = TRUE) &
    grepl(grp2, basename(files), fixed = TRUE)]
  if (length(hit) >= 1) return(hit[1])
  ct_file <- file.path(ddir, "comparisonTable.rds")
  if (!file.exists(ct_file)) return(NA_character_)
  ct <- readRDS(ct_file)
  i <- which(paste0(ct$grp1, " vs ", ct$grp2) == paste(grp1, "vs", grp2))
  if (length(i) != 1 || length(files) < i) return(NA_character_)
  files[i]
}

read_regions <- function(f) {
  dm <- read.delim(f)
  gr <- GenomicRanges::GRanges(
    seqnames = dm$chrom,
    ranges = IRanges::IRanges(start = dm$chromStart + 1, end = dm$chromEnd)
  )
  S4Vectors::mcols(gr)$log2FoldChange <- dm$log2FoldChange
  S4Vectors::mcols(gr)$padj <- dm$padj
  gr
}

# gain = more accessible in grp1, loss = less accessible in grp1
split_dars <- function(gr) {
  sig <- !is.na(gr$padj) & gr$padj < padj_cut & abs(gr$log2FoldChange) > l2fc_cut
  list(
    gain = gr[sig & gr$log2FoldChange > 0],
    loss = gr[sig & gr$log2FoldChange < 0]
  )
}

## ---- LOLA for one region set ------------------------------------------------
# Universe = regions tested in that comparison, not the genome. Both models then
# ask the same question: among the regions we could have called, are the ones we
# did call enriched?
run_lola_set <- function(userSets, universe, label) {
  userSets <- userSets[vapply(userSets, length, integer(1)) > 0]
  if (length(userSets) == 0) {
    message("  ", label, ": no DARs, skipping LOLA")
    return(NULL)
  }
  message("  ", label, ": ",
    paste(sprintf("%s=%d", names(userSets), vapply(userSets, length, integer(1))),
      collapse = ", "), " regions")
  res <- LOLA::runLOLA(
    userSets = GenomicRanges::GRangesList(userSets),
    userUniverse = universe,
    regionDB = lolaDb,
    cores = 1
  )
  res <- as.data.frame(res)
  # older LOLA versions have no qValue column
  if (!"qValue" %in% names(res)) {
    res$qValue <- stats::p.adjust(10^(-res$pValueLog), method = "BH")
  }
  res$log2OR <- log2(pmax(res$oddsRatio, .Machine$double.eps))
  res$model <- label
  res
}

## ---- Compare adjusted vs unadjusted for one target --------------------------
term_key <- function(d) paste(d$collection, d$filename, d$userSet, sep = " | ")

process_target <- function(tg) {
  comp_label <- paste0(tg$cell, " | ", tg$grp1, " vs ", tg$grp2)
  message("=== ", comp_label)

  adj_dir <- adjusted_dir_for(tg$cell, tg$grp1, tg$grp2, focus_set)
  if (is.na(adj_dir)) {
    message("  no adjusted run found; run 07_3 first")
    return(NULL)
  }
  adj_files <- list_diff_tabs(adj_dir)
  unadj_file <- unadjusted_tab_for(tg$dir, tg$cell, tg$grp1, tg$grp2)
  if (length(adj_files) == 0 || is.na(unadj_file)) {
    message("  missing diffTab on one side; skipping")
    return(NULL)
  }

  adj_gr <- read_regions(adj_files[1])
  unadj_gr <- read_regions(unadj_file)

  # Universe = regions tested by both models, matched on exact coordinates.
  # GenomicRanges::intersect() merges overlapping ranges and would return
  # boundaries that match no tested peak.
  region_id <- function(gr) paste0(
    as.character(GenomicRanges::seqnames(gr)), ":",
    GenomicRanges::start(gr), "-", GenomicRanges::end(gr))
  id_adj <- region_id(adj_gr)
  id_unadj <- region_id(unadj_gr)
  shared_ids <- intersect(id_adj, id_unadj)
  if (length(shared_ids) == 0) {
    message("  the two models share no tested region; skipping")
    return(NULL)
  }
  adj_keep <- adj_gr[id_adj %in% shared_ids]
  unadj_keep <- unadj_gr[id_unadj %in% shared_ids]
  universe <- GenomicRanges::granges(adj_keep)
  message("  universe: ", length(universe), " regions tested by both models (",
    length(adj_gr), " adjusted / ", length(unadj_gr), " unadjusted)")

  adj_sets_gr <- split_dars(adj_keep)
  unadj_sets_gr <- split_dars(unadj_keep)

  res_adj <- run_lola_set(adj_sets_gr, universe, "QC-adjusted")
  res_unadj <- run_lola_set(unadj_sets_gr, universe, "Unadjusted")
  if (is.null(res_adj) || is.null(res_unadj)) {
    message("  one side has no DARs; nothing to compare")
    return(NULL)
  }

  res_adj$key <- term_key(res_adj)
  res_unadj$key <- term_key(res_unadj)
  both <- dplyr::inner_join(
    res_unadj %>% dplyr::select(key, collection, filename, description, userSet,
      log2OR_unadj = log2OR, qValue_unadj = qValue, support_unadj = support),
    res_adj %>% dplyr::select(key,
      log2OR_adj = log2OR, qValue_adj = qValue, support_adj = support),
    by = "key"
  ) %>%
    dplyr::mutate(
      cell = tg$cell, comparison = paste(tg$grp1, "vs", tg$grp2),
      sig_unadj = qValue_unadj < 0.05,
      sig_adj   = qValue_adj < 0.05,
      agreement = dplyr::case_when(
        sig_unadj & sig_adj ~ "significant in both",
        sig_unadj ~ "unadjusted only",
        sig_adj ~ "adjusted only",
        TRUE ~ "neither"
      )
    )

  n_both <- sum(both$agreement == "significant in both")
  n_un <- sum(both$agreement == "unadjusted only")
  n_ad <- sum(both$agreement == "adjusted only")
  # cor() returns NA on any non-finite value; LOLA gives OR = 0 for some terms
  fin <- is.finite(both$log2OR_unadj) & is.finite(both$log2OR_adj)
  r <- if (sum(fin) >= 3) {
    suppressWarnings(stats::cor(both$log2OR_unadj[fin], both$log2OR_adj[fin]))
  } else NA_real_
  if (sum(fin) < nrow(both)) {
    message("  ", nrow(both) - sum(fin),
      " term(s) had a non-finite log2 odds ratio and were excluded from the correlation")
  }
  message(sprintf("  enrichment terms: %d significant in both, %d unadjusted only, %d adjusted only; log2OR r = %.3f",
    n_both, n_un, n_ad, r))

  list(
    comparison = comp_label,
    table = both,
    summary = data.frame(
      cell = tg$cell, comparison = paste(tg$grp1, "vs", tg$grp2),
      adj_set = focus_set,
      universe_regions = length(universe),
      DARs_unadjusted = sum(vapply(unadj_sets_gr, length, integer(1))),
      DARs_adjusted = sum(vapply(adj_sets_gr, length, integer(1))),
      terms_tested = nrow(both),
      terms_sig_both = n_both,
      terms_sig_unadjusted_only = n_un,
      terms_sig_adjusted_only = n_ad,
      terms_recovered_pct = ifelse((n_both + n_un) > 0,
        round(100 * n_both / (n_both + n_un), 1), NA_real_),
      # ORs are not comparable across sets of very different size: a bigger DAR
      # set includes weaker regions and dilutes the enrichment
      DAR_size_ratio = round(sum(vapply(adj_sets_gr, length, integer(1))) /
        pmax(1, sum(vapply(unadj_sets_gr, length, integer(1)))), 2),
      log2OR_pearson = round(r, 3),
      stringsAsFactors = FALSE
    )
  )
}

results <- lapply(targets, function(tg) {
  tryCatch(process_target(tg), error = function(e) {
    message("  FAILED: ", conditionMessage(e))
    NULL
  })
})
results <- Filter(Negate(is.null), results)
if (length(results) == 0) {
  stop("No target could be compared. Check that 07_3 has produced adjusted runs.")
}

lola_summary <- dplyr::bind_rows(lapply(results, `[[`, "summary"))
print(lola_summary)
write.csv(lola_summary,
  file.path(fig_dir, "confounder_adjusted_lola_summary.csv"), row.names = FALSE)

lola_terms <- dplyr::bind_rows(lapply(results, `[[`, "table"))
write.csv(
  lola_terms %>% dplyr::filter(sig_unadj | sig_adj) %>% dplyr::select(-key),
  file.path(fig_dir, "confounder_adjusted_lola_terms.csv"), row.names = FALSE
)

## ---- Figure 1: enrichment odds ratios, adjusted vs unadjusted ---------------
plot_terms <- lola_terms %>%
  dplyr::filter(collection %in% collections_of_interest) %>%
  dplyr::mutate(
    panel = paste0(cell, "\n", comparison, " (", userSet, ")"),
    agreement = factor(agreement,
      levels = c("neither", "unadjusted only", "adjusted only", "significant in both"))
  ) %>%
  dplyr::arrange(agreement)

if (nrow(plot_terms) > 0) {
  cor_lab <- plot_terms %>%
    dplyr::group_by(panel) %>%
    dplyr::summarise(
      r = {
        f <- is.finite(log2OR_unadj) & is.finite(log2OR_adj)
        if (sum(f) >= 3) suppressWarnings(stats::cor(log2OR_unadj[f], log2OR_adj[f])) else NA_real_
      },
      n_both = sum(agreement == "significant in both"), .groups = "drop") %>%
    dplyr::mutate(lab = paste0("r = ", sprintf("%.3f", r), "\n", n_both, " sig. in both"))

  p_or <- ggplot(plot_terms, aes(x = log2OR_unadj, y = log2OR_adj, colour = agreement)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60") +
    geom_point(size = 1, alpha = 0.7) +
    geom_text(data = cor_lab, inherit.aes = FALSE,
      aes(x = -Inf, y = Inf, label = lab), hjust = -0.1, vjust = 1.2, size = 3) +
    facet_wrap(~panel, scales = "free") +
    scale_colour_manual(values = c(
      "neither" = "grey85", "unadjusted only" = "#3C5488",
      "adjusted only" = "#E64B35", "significant in both" = "#4A2377")) +
    guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1))) +
    labs(
      title = "LOLA enrichment is unchanged by technical-covariate adjustment",
      subtitle = paste0(
        "One point per region set in ",
        paste(collections_of_interest, collapse = ", "),
        ". Dashed line, y = x. Significance at LOLA q < 0.05.\\n",
        "Both models use the same DAR definition (|log2FC| > ", l2fc_cut,
        ", padj < ", padj_cut, ") and the same universe of tested regions."),
      x = "log2 odds ratio, unadjusted DARs",
      y = "log2 odds ratio, QC-adjusted DARs", colour = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(legend.position = "top", plot.subtitle = element_text(size = 8),
      strip.background = element_blank())

  n_panels <- dplyr::n_distinct(plot_terms$panel)
  ggsave(file.path(fig_dir, "confounder_adjusted_lola_scatter.pdf"), p_or,
    width = max(8, 3.4 * min(3, n_panels)),
    height = 3.4 * ceiling(n_panels / 3) + 2, limitsize = FALSE)
}

## ---- Figure 2: top terms side by side ---------------------------------------
# The scatter says whether enrichment changed; this names the terms.
top_terms <- lola_terms %>%
  dplyr::filter(collection %in% collections_of_interest, sig_unadj | sig_adj) %>%
  dplyr::mutate(panel = paste0(cell, " | ", comparison, " (", userSet, ")")) %>%
  dplyr::group_by(panel) %>%
  dplyr::slice_max(order_by = pmax(log2OR_unadj, log2OR_adj), n = 15,
    with_ties = FALSE) %>%
  dplyr::ungroup()

if (nrow(top_terms) > 0) {
  top_long <- top_terms %>%
    dplyr::select(panel, description, filename, log2OR_unadj, log2OR_adj) %>%
    tidyr::pivot_longer(c(log2OR_unadj, log2OR_adj),
      names_to = "model", values_to = "log2OR") %>%
    dplyr::mutate(
      model = factor(ifelse(model == "log2OR_unadj", "Unadjusted", "QC-adjusted"),
        levels = c("Unadjusted", "QC-adjusted")),
      term = ifelse(is.na(description) | description == "",
        basename(as.character(filename)), as.character(description))
    )

  p_top <- ggplot(top_long, aes(x = reorder(term, log2OR), y = log2OR, fill = model)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    coord_flip() +
    facet_wrap(~panel, scales = "free") +
    scale_fill_manual(values = c("Unadjusted" = "#3C5488", "QC-adjusted" = "#E64B35")) +
    labs(title = "Top LOLA terms before and after covariate adjustment",
      subtitle = "Top 15 region sets per comparison, ranked by the larger of the two odds ratios",
      x = NULL, y = "log2 odds ratio", fill = NULL) +
    theme_classic(base_size = 10) +
    theme(legend.position = "top", plot.subtitle = element_text(size = 8),
      strip.background = element_blank(), axis.text.y = element_text(size = 7))

  n_panels <- dplyr::n_distinct(top_long$panel)
  ggsave(file.path(fig_dir, "confounder_adjusted_lola_top_terms.pdf"), p_top,
    width = max(9, 5 * min(2, n_panels)),
    height = 4.5 * ceiling(n_panels / 2) + 1.5, limitsize = FALSE)
}

message("Done. LOLA comparison figures and tables in ", fig_dir)

#####################################################################
