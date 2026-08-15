#!/usr/bin/env Rscript

#####################################################################
# 07_4_confounder_adjusted_DAR_plots.R
# created on 2026-08-12 by Irem B. Gunduz
# Plot-only companion to 07_3_confounder_adjusted_DARs.R
#
# Why this exists
#  - 07_3 builds its summary table inside a tryCatch and only draws the
#    figure when that table is non-empty. If a comparison errors AFTER the
#    adjusted DESeq2 run has already written its diffTab (or if `res` comes
#    back NULL for any other reason), the adjusted results sit on disk but
#    no PDF is produced.
#  - This script does the reading and plotting half only: it discovers the
#    "<cell>__<grp1>_vs_<grp2>__adjusted" directories 07_3 has already
#    written, pairs each with the corresponding unadjusted diffTab from the
#    original ChrAccR run, and writes the summary table and figures.
#  - Nothing is re-fitted. No ArchR project, no ChrAccR object, no DESeq2.
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

repo_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript"
source(file.path(repo_dir, "utils/helpers.R")) # cutL0.5fc2Padj05

# must match 07_3
out_dir <- "/icbb/projects/igunduz/finalize_echo_050824/confounder_adjusted/"
qc_adj_cols <- c("mean_TSS", "mean_FRIP")

# Figures and the small summary tables belong in the repo figures/ directory with
# every other panel; only the bulky per-region table stays in the scratch dir.
fig_dir <- file.path(repo_dir, "figures")
if (!dir.exists(fig_dir)) dir.create(fig_dir, recursive = TRUE)

# The sweep in 07_3 covers every cell type, which is far too many facets for one
# readable scatter. Draw the panels with the most unadjusted DARs and say so.
max_scatter_panels <- 24

# 07_3 fits each comparison under several adjustment designs. The overlap and
# scatter figures show ONE design so they stay readable; the design-comparison
# figure below covers the rest.
focus_set <- "TSS_FRIP"

# original (unadjusted) ChrAccR analysis directories, as in 07_3
covid_dir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/"
other_dir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_2023-10-02/"

# the C19 comparisons live in the covid run, everything else in the general run
source_dir_for <- function(grp1) if (grepl("^C19", grp1)) covid_dir else other_dir

stopifnot(dir.exists(out_dir))

## ---- 1. Discover the adjusted runs 07_3 has already written -----------------
# 07_3 names them paste0(cell, "__", tag, "__adjusted") where tag is the
# comparison with every non-alphanumeric character replaced by "_", e.g.
# "Mono_CD14__C19_sev_vs_C19_ctrl__adjusted".
adj_dirs <- list.dirs(out_dir, recursive = FALSE, full.names = TRUE)
adj_dirs <- adj_dirs[grepl("__adjusted$|__adj-", basename(adj_dirs))]
if (length(adj_dirs) == 0) {
  stop("No '*__adjusted' directories found in ", out_dir,
    " -- run 07_3_confounder_adjusted_DARs.R first.")
}
message("Found ", length(adj_dirs), " adjusted run(s):\n  ",
  paste(basename(adj_dirs), collapse = "\n  "))

# 07_3 writes "<cell>__<grp1>_vs_<grp2>__adj-<set>" once the adjustment-set grid
# is in play. The older "<cell>__<grp1>_vs_<grp2>__adjusted" directories are the
# TSS + FRIP design, so they are read as that set rather than ignored.
parse_adj_dir <- function(d) {
  parts <- strsplit(basename(d), "__", fixed = TRUE)[[1]]
  if (length(parts) < 3) return(NULL)
  cell <- parts[1]
  tag <- parts[2]
  suffix <- parts[3]
  adj_set <- if (identical(suffix, "adjusted")) "TSS_FRIP" else sub("^adj-", "", suffix)
  grps <- strsplit(tag, "_vs_", fixed = TRUE)[[1]]
  if (length(grps) != 2) return(NULL)
  list(dir = d, cell = cell, grp1 = grps[1], grp2 = grps[2],
    adj_set = adj_set, comparison = paste(grps[1], "vs", grps[2]))
}

## ---- 2. Locate diffTab files ------------------------------------------------
# Prefer the archrPeaks region set (the one the manuscript reports); fall back to
# any diffTab if that naming is not used in a given run.
list_diff_tabs <- function(dir) {
  f <- list.files(dir, pattern = "diffTab.*archrPeaks.*\\.tsv$",
    recursive = TRUE, full.names = TRUE)
  if (length(f) == 0) {
    f <- list.files(dir, pattern = "diffTab.*\\.tsv$",
      recursive = TRUE, full.names = TRUE)
  }
  sort(f)
}

# The unadjusted run holds every comparison in one directory, so pick the file
# for this comparison by name; if the filenames do not carry the group labels,
# fall back to the comparisonTable row index (the indexing 07_3 used).
pick_unadjusted <- function(anaDir, cell, grp1, grp2) {
  ddir <- file.path(anaDir, cell, "reports", "differential_data")
  if (!dir.exists(ddir)) {
    warning("no differential_data directory at ", ddir)
    return(NA_character_)
  }
  files <- list_diff_tabs(ddir)
  if (length(files) == 0) {
    warning("no diffTab files under ", ddir)
    return(NA_character_)
  }
  bn <- basename(files)
  hit <- files[grepl(grp1, bn, fixed = TRUE) & grepl(grp2, bn, fixed = TRUE)]
  if (length(hit) >= 1) {
    if (length(hit) > 1) {
      message("  several unadjusted tables match ", grp1, " vs ", grp2,
        "; using ", basename(hit[1]))
    }
    return(hit[1])
  }
  ct_file <- file.path(ddir, "comparisonTable.rds")
  if (!file.exists(ct_file)) {
    warning("cannot identify the unadjusted table for ", grp1, " vs ", grp2,
      " (no name match and no comparisonTable.rds)")
    return(NA_character_)
  }
  ct <- readRDS(ct_file)
  i <- which(paste0(ct$grp1, " vs ", ct$grp2) == paste(grp1, "vs", grp2))
  if (length(i) != 1 || length(files) < i) {
    warning("comparisonTable does not resolve ", grp1, " vs ", grp2)
    return(NA_character_)
  }
  message("  resolved the unadjusted table by comparisonTable index (", i, ")")
  files[i]
}

key_of <- function(cell, comparison, adj_set = focus_set) {
  paste(cell, comparison, adj_set, sep = " | ")
}

read_diff <- function(f) {
  dm <- read.delim(f)
  isDiff <- cutL0.5fc2Padj05(dm[, c("log2FoldChange", "padj")])
  isDiff[is.na(isDiff)] <- FALSE
  data.frame(
    id = paste0(dm$chrom, ":", dm$chromStart + 1, "-", dm$chromEnd),
    log2FC = dm$log2FoldChange, padj = dm$padj, isDiff = isDiff,
    stringsAsFactors = FALSE
  )
}

## ---- 3. Summarise each comparison -------------------------------------------
# Pass 1 computes the summary for EVERY comparison and deliberately discards each
# merged table afterwards. Holding all of them at once would be tens of millions
# of rows once the sweep covers every cell type; the scatter panels re-read the
# handful of tables they need in pass 2.
build_merged <- function(j) {
  adj_files <- list_diff_tabs(j$dir)
  if (length(adj_files) == 0) {
    message("  skipped: no diffTab under ", j$dir)
    return(NULL)
  }
  if (length(adj_files) > 1) {
    message("  ", length(adj_files), " adjusted tables present; using ",
      basename(adj_files[1]))
  }
  adj <- read_diff(adj_files[1])

  unadj_file <- pick_unadjusted(source_dir_for(j$grp1), j$cell, j$grp1, j$grp2)
  if (is.na(unadj_file)) {
    message("  skipped: no matching unadjusted table")
    return(NULL)
  }
  unadj <- read_diff(unadj_file)

  m <- merge(unadj, adj, by = "id", suffixes = c("_unadj", "_adj"))
  if (nrow(m) == 0) {
    message("  skipped: adjusted and unadjusted tables share no region ids")
    return(NULL)
  }
  m$cell <- j$cell
  m$comparison <- j$comparison
  m
}

# Regions called a DAR by either model are kept from pass 1 (a small fraction of
# the merge), so the per-region table covers every comparison rather than only
# the ones drawn in the scatter.
dar_list <- list()

# Per-comparison counts broken down by overlap category and direction, used for
# the stacked bars. Computed in pass 1 so nothing has to be re-read for them.
cat_list <- list()

summarise_one <- function(j) {
  message("=== ", j$cell, " | ", j$comparison, "  {", j$adj_set, "}")
  m <- build_merged(j)
  if (is.null(m)) return(NULL)
  m$adj_set <- j$adj_set

  keep <- m$isDiff_unadj | m$isDiff_adj
  if (any(keep)) {
    dar_list[[key_of(j$cell, j$comparison, j$adj_set)]] <<- m[keep, , drop = FALSE]
  }

  a <- m$id[m$isDiff_unadj]
  b <- m$id[m$isDiff_adj]
  shared <- intersect(a, b)
  conc <- if (length(shared)) {
    mean(sign(m$log2FC_unadj[match(shared, m$id)]) ==
      sign(m$log2FC_adj[match(shared, m$id)])) * 100
  } else NA_real_

  message("  ", length(a), " unadjusted / ", length(b), " adjusted DARs, ",
    length(shared), " shared")

  # Classify every DAR by which model called it and in which direction. Shared
  # and unadjusted-only regions take their direction from the unadjusted fit,
  # adjusted-only regions from the adjusted fit; the sign concordance above says
  # how often the two agree where both called a region.
  cls <- dplyr::case_when(
    m$isDiff_unadj & m$isDiff_adj ~ "Shared",
    m$isDiff_unadj ~ "Unadjusted only",
    m$isDiff_adj ~ "Adjusted only",
    TRUE ~ NA_character_
  )
  sel <- !is.na(cls)
  l2 <- ifelse(m$isDiff_unadj, m$log2FC_unadj, m$log2FC_adj)
  dirn <- ifelse(l2 > 0, "Hyper-accessible", "Hypo-accessible")
  cnt_df <- as.data.frame(
    table(
      overlap = factor(cls[sel],
        levels = c("Shared", "Unadjusted only", "Adjusted only")),
      direction = factor(dirn[sel],
        levels = c("Hyper-accessible", "Hypo-accessible"))
    ),
    stringsAsFactors = FALSE
  )
  names(cnt_df)[names(cnt_df) == "Freq"] <- "n"
  cnt_df$cell <- j$cell
  cnt_df$comparison <- j$comparison
  cnt_df$adj_set <- j$adj_set
  cat_list[[key_of(j$cell, j$comparison, j$adj_set)]] <<- cnt_df

  get_n <- function(o, d) {
    v <- cnt_df$n[cnt_df$overlap == o & cnt_df$direction == d]
    if (length(v) == 1) v else 0L
  }

  data.frame(
    cell = j$cell, comparison = j$comparison, adj_set = j$adj_set,
    DARs_unadjusted = length(a), DARs_adjusted = length(b),
    shared = length(shared),
    recovered_pct = ifelse(length(a) > 0, round(100 * length(shared) / length(a), 1), NA),
    sign_concordance_pct = round(conc, 1),
    lfc_pearson_all = round(cor(m$log2FC_unadj, m$log2FC_adj, use = "complete.obs"), 3),
    n_regions_tested = nrow(m),
    # direction breakdown, so the stacked figure can be read off the table
    shared_hyper = get_n("Shared", "Hyper-accessible"),
    shared_hypo = get_n("Shared", "Hypo-accessible"),
    unadj_only_hyper = get_n("Unadjusted only", "Hyper-accessible"),
    unadj_only_hypo = get_n("Unadjusted only", "Hypo-accessible"),
    adj_only_hyper = get_n("Adjusted only", "Hyper-accessible"),
    adj_only_hypo = get_n("Adjusted only", "Hypo-accessible"),
    stringsAsFactors = FALSE
  )
}

jobs <- Filter(Negate(is.null), lapply(adj_dirs, parse_adj_dir))
res <- do.call(rbind, lapply(jobs, function(j) {
  out <- tryCatch(summarise_one(j),
    error = function(e) {
      message("  FAILED: ", conditionMessage(e))
      NULL
    }
  )
  gc()
  out
}))

if (is.null(res) || nrow(res) == 0) {
  stop("None of the adjusted runs could be paired with an unadjusted table; ",
    "nothing to plot. See the messages above for which step failed.")
}

## ---- 3b. Split the results: one design for the detail figures ---------------
# res_all keeps every comparison x adjustment design and is what gets written to
# the summary table; res is narrowed to a single design so the overlap and
# scatter figures stay readable.
res_all <- res
if (!focus_set %in% res_all$adj_set) {
  message("focus_set '", focus_set, "' not present; using the most common design")
  focus_set <- names(sort(table(res_all$adj_set), decreasing = TRUE))[1]
}
res <- res_all[res_all$adj_set == focus_set, , drop = FALSE]
if (nrow(res) == 0) {
  stop("No comparison used the '", focus_set, "' design; nothing to plot.")
}
message("Detail figures use the '", focus_set, "' design (", nrow(res),
  " comparisons); ", nrow(res_all), " comparison x design fits in total")
cat_list <- cat_list[vapply(cat_list, function(d) d$adj_set[1] == focus_set, logical(1))]

print(res_all)
write.csv(res_all, file.path(fig_dir, "confounder_adjusted_DAR_summary.csv"),
  row.names = FALSE)

## ---- 4. Figure 1: overlap and direction of the DAR calls --------------------
# A side-by-side count of "unadjusted DARs" against "adjusted DARs" is hard to
# read, because two similar-height bars can still be two different region sets.
# This shows the partition instead: every DAR is one of shared / unadjusted-only
# / adjusted-only, stacked, and split by direction. A result that is robust to
# the adjustment looks like a tall "Shared" block with thin slivers either side,
# balanced across the two direction panels.
cat_all <- dplyr::bind_rows(cat_list) %>%
  dplyr::mutate(
    label = paste0(cell, " | ", comparison),
    overlap = factor(overlap,
      levels = c("Unadjusted only", "Shared", "Adjusted only"))
  )

# order the comparisons by total DAR count, largest at the top
bar_order <- cat_all %>%
  dplyr::group_by(label) %>%
  dplyr::summarise(total = sum(n), .groups = "drop") %>%
  dplyr::arrange(total)
cat_all$label <- factor(cat_all$label, levels = bar_order$label)

# per-panel totals, annotated at the end of each bar
tot_lab <- cat_all %>%
  dplyr::group_by(label, direction) %>%
  dplyr::summarise(n = sum(n), .groups = "drop")

overlap_cols <- c(
  "Unadjusted only" = "#3C5488", # called only without the QC covariates
  "Shared"          = "#B8B8B8", # called by both models
  "Adjusted only"   = "#E64B35"  # called only with the QC covariates
)

p_stack <- ggplot(cat_all, aes(x = label, y = n, fill = overlap)) +
  geom_col(width = 0.7, colour = "grey25", linewidth = 0.15) +
  geom_text(data = tot_lab, aes(x = label, y = n, label = n),
    inherit.aes = FALSE, hjust = -0.15, size = 2.4, colour = "grey20") +
  coord_flip() +
  facet_wrap(~direction) +
  scale_fill_manual(values = overlap_cols) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(
    title = "Overlap of DAR calls with and without technical-covariate adjustment",
    subtitle = paste0(
      "Adjusted model adds ", paste(qc_adj_cols, collapse = " + "),
      " to the differential design. DAR = |log2FC| > 0.5 and padj < 0.05.\n",
      "Bar segments partition the DARs of each comparison; the number at the end of each bar is the union. ",
      nrow(res), " comparisons across ", dplyr::n_distinct(res$cell), " cell types."
    ),
    x = NULL, y = "Number of DARs", fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", plot.subtitle = element_text(size = 8),
    strip.background = element_blank(), strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "confounder_adjusted_DARs.pdf"), p_stack,
  width = 11, height = max(5, 0.30 * nrow(res) + 2.2), limitsize = FALSE)

# The same partition scaled to 100%, so a comparison with 60 DARs is as readable
# as one with 6000. Read the two together: this one shows what fraction survives
# the adjustment, the one above shows how many regions that is.
p_stack_pct <- ggplot(cat_all, aes(x = label, y = n, fill = overlap)) +
  geom_col(width = 0.7, colour = "grey25", linewidth = 0.15, position = "fill") +
  coord_flip() +
  facet_wrap(~direction) +
  scale_fill_manual(values = overlap_cols) +
  scale_y_continuous(labels = scales::percent,
    expand = expansion(mult = c(0, 0.02))) +
  labs(
    title = "Overlap of DAR calls, as a share of each comparison",
    subtitle = "Same partition as the count figure, scaled to 100% per comparison and direction.",
    x = NULL, y = "Share of DARs", fill = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", plot.subtitle = element_text(size = 8),
    strip.background = element_blank(), strip.text = element_text(face = "bold"))

ggsave(file.path(fig_dir, "confounder_adjusted_DARs_proportion.pdf"), p_stack_pct,
  width = 11, height = max(5, 0.30 * nrow(res) + 2.2), limitsize = FALSE)

## ---- 5. Figure 2: per-region effect sizes, adjusted vs unadjusted -----------
# The count bar chart alone does not show that the SAME regions move the same
# way. This panel does: every tested region, unadjusted log2FC on x against
# adjusted log2FC on y, with the Pearson correlation annotated per comparison.
#
# Across the full sweep there are far more comparisons than fit in one readable
# figure, so the panels are the max_scatter_panels comparisons with the most
# unadjusted DARs -- the ones the manuscript's claims actually rest on. The
# summary table above covers all of them.
res_scatter <- res[order(-res$DARs_unadjusted), , drop = FALSE]
res_scatter <- utils::head(res_scatter, max_scatter_panels)
message("Drawing scatter panels for ", nrow(res_scatter), " of ", nrow(res),
  " comparisons (most unadjusted DARs first)")

want_keys <- key_of(res_scatter$cell, res_scatter$comparison)
sel_jobs <- Filter(function(j) key_of(j$cell, j$comparison) %in% want_keys, jobs)

merged_df <- dplyr::bind_rows(lapply(sel_jobs, function(j) {
  m <- tryCatch(build_merged(j), error = function(e) NULL)
  gc()
  m
})) %>%
  dplyr::mutate(
    panel = paste0(cell, "\n", comparison),
    status = dplyr::case_when(
      isDiff_unadj & isDiff_adj ~ "DAR in both",
      isDiff_unadj ~ "DAR unadjusted only",
      isDiff_adj ~ "DAR adjusted only",
      TRUE ~ "not significant"
    )
  )
merged_df$status <- factor(merged_df$status,
  levels = c("not significant", "DAR unadjusted only", "DAR adjusted only", "DAR in both"))
merged_df <- merged_df[order(merged_df$status), ] # significant points drawn on top

# A vector PDF with one point per tested region is unusably large, so thin the
# grey background: every significant region is kept, the non-significant ones are
# subsampled to at most 30k per panel. This changes the density of the grey cloud
# only, not any reported number.
set.seed(12)
max_bg <- 30000
sig_pts <- merged_df[merged_df$status != "not significant", , drop = FALSE]
bg_pts <- merged_df[merged_df$status == "not significant", , drop = FALSE] %>%
  dplyr::group_by(panel) %>%
  dplyr::slice_sample(n = max_bg) %>% # slice_sample keeps all rows if n > group size
  dplyr::ungroup()
plot_pts <- dplyr::bind_rows(bg_pts, sig_pts) # background first, DARs on top

cor_lab <- res %>%
  dplyr::mutate(panel = paste0(cell, "\n", comparison),
    lab = paste0("r = ", sprintf("%.3f", lfc_pearson_all),
      "\n", shared, "/", DARs_unadjusted, " recovered"))

p_scatter <- ggplot(plot_pts, aes(x = log2FC_unadj, y = log2FC_adj, colour = status)) +
  geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = "dashed") +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_text(data = cor_lab, inherit.aes = FALSE,
    aes(x = -Inf, y = Inf, label = lab), hjust = -0.1, vjust = 1.2, size = 3) +
  facet_wrap(~panel, scales = "free", ncol = min(4, max(1, nrow(res_scatter)))) +
  scale_colour_manual(values = c(
    "not significant" = "grey80", "DAR unadjusted only" = "#3C5488",
    "DAR adjusted only" = "#E64B35", "DAR in both" = "#4A2377"
  )) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(
    title = "Per-region effect sizes are unchanged by technical-covariate adjustment",
    subtitle = paste0("Dashed line, y = x. Each point is one tested region. ",
      nrow(res_scatter), " of ", nrow(res),
      " comparisons shown, ranked by unadjusted DAR count."),
    x = "log2 fold change, unadjusted model",
    y = "log2 fold change, QC-adjusted model", colour = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", plot.subtitle = element_text(size = 9))

n_col <- min(4, max(1, nrow(res_scatter)))
n_row <- ceiling(nrow(res_scatter) / n_col)
ggsave(file.path(fig_dir, "confounder_adjusted_DAR_l2fc_scatter.pdf"), p_scatter,
  width = 2.6 * n_col + 1.5, height = 2.6 * n_row + 1.5, limitsize = FALSE)

# Per-region table for every region called a DAR by either model, across ALL
# comparisons, so each number in the two figures can be traced back. Regions
# significant in neither model are left out; the full merge would be tens of
# millions of rows. This one stays in the scratch dir, not the repo.
dar_regions <- dplyr::bind_rows(dar_list)
if (nrow(dar_regions) > 0) {
  dar_regions <- dar_regions %>%
    dplyr::mutate(status = dplyr::case_when(
      isDiff_unadj & isDiff_adj ~ "DAR in both",
      isDiff_unadj ~ "DAR unadjusted only",
      TRUE ~ "DAR adjusted only"
    )) %>%
    dplyr::select(cell, comparison, id, log2FC_unadj, padj_unadj, isDiff_unadj,
      log2FC_adj, padj_adj, isDiff_adj, status)
}
write.csv(dar_regions,
  file.path(out_dir, "confounder_adjusted_DAR_regions.csv"), row.names = FALSE)
message("Wrote ", nrow(dar_regions), " region rows (DAR in either model)")

## ---- 6. Figure 3: does the choice of adjustment covariates matter? ----------
# For the comparisons fitted under more than one design, plot the DAR count per
# design against the unadjusted count (dashed line). Flat across designs means
# the result does not depend on which covariates were adjusted for -- which is
# the claim the reviewer is asking us to support.
multi <- res_all %>%
  dplyr::group_by(cell, comparison) %>%
  dplyr::filter(dplyr::n_distinct(adj_set) > 1) %>%
  dplyr::ungroup()

if (nrow(multi) > 0) {
  multi <- multi %>%
    dplyr::mutate(
      panel = paste0(cell, " | ", comparison),
      adj_set = factor(adj_set, levels = unique(adj_set[order(nchar(adj_set), adj_set)]))
    )
  unadj_ref <- multi %>%
    dplyr::group_by(panel) %>%
    dplyr::summarise(DARs_unadjusted = dplyr::first(DARs_unadjusted), .groups = "drop")

  p_design <- ggplot(multi, aes(x = adj_set, y = DARs_adjusted)) +
    geom_hline(data = unadj_ref, aes(yintercept = DARs_unadjusted),
      linetype = "dashed", colour = "grey50") +
    geom_col(fill = "#3C5488", width = 0.7) +
    geom_text(aes(label = DARs_adjusted), vjust = -0.35, size = 2.6) +
    facet_wrap(~panel, scales = "free_y") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "DAR count by adjustment design",
      subtitle = paste0(
        "Dashed line = unadjusted DAR count. Designs that could not be fitted ",
        "(constant, collinear or nested covariates, or too few residual\n",
        "degrees of freedom) are absent by design and are listed with their ",
        "reason in confounder_adjusted_DAR_fit_log.csv, written by 07_3."),
      x = "Adjustment covariates", y = "Number of DARs"
    ) +
    theme_classic(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
      plot.subtitle = element_text(size = 8),
      strip.background = element_blank(), strip.text = element_text(face = "bold"))

  n_panels <- dplyr::n_distinct(multi$panel)
  ggsave(file.path(fig_dir, "confounder_adjusted_DARs_by_design.pdf"), p_design,
    width = max(8, 3.2 * min(3, n_panels)),
    height = 3.2 * ceiling(n_panels / 3) + 2, limitsize = FALSE)
  message("Design-comparison figure covers ", n_panels, " comparison(s)")
} else {
  message("Only one adjustment design present; skipping the design-comparison ",
    "figure. Set combo_scope in 07_3 to fit more designs.")
}

message("Done. Figures and summary tables in ", fig_dir,
  "; per-region table in ", out_dir)

#####################################################################
