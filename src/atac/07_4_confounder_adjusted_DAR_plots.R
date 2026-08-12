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
adj_dirs <- adj_dirs[grepl("__adjusted$", basename(adj_dirs))]
if (length(adj_dirs) == 0) {
  stop("No '*__adjusted' directories found in ", out_dir,
    " -- run 07_3_confounder_adjusted_DARs.R first.")
}
message("Found ", length(adj_dirs), " adjusted run(s):\n  ",
  paste(basename(adj_dirs), collapse = "\n  "))

parse_adj_dir <- function(d) {
  parts <- strsplit(basename(d), "__", fixed = TRUE)[[1]]
  if (length(parts) < 3) return(NULL)
  cell <- parts[1]
  tag <- parts[2]
  grps <- strsplit(tag, "_vs_", fixed = TRUE)[[1]]
  if (length(grps) != 2) return(NULL)
  list(dir = d, cell = cell, grp1 = grps[1], grp2 = grps[2],
    comparison = paste(grps[1], "vs", grps[2]))
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
merged_list <- list()

summarise_one <- function(j) {
  message("=== ", j$cell, " | ", j$comparison)

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
  merged_list[[paste(j$cell, j$comparison)]] <<- m

  a <- m$id[m$isDiff_unadj]
  b <- m$id[m$isDiff_adj]
  shared <- intersect(a, b)
  conc <- if (length(shared)) {
    mean(sign(m$log2FC_unadj[match(shared, m$id)]) ==
      sign(m$log2FC_adj[match(shared, m$id)])) * 100
  } else NA_real_

  message("  ", length(a), " unadjusted / ", length(b), " adjusted DARs, ",
    length(shared), " shared")

  data.frame(
    cell = j$cell, comparison = j$comparison,
    DARs_unadjusted = length(a), DARs_adjusted = length(b),
    shared = length(shared),
    recovered_pct = ifelse(length(a) > 0, round(100 * length(shared) / length(a), 1), NA),
    sign_concordance_pct = round(conc, 1),
    lfc_pearson_all = round(cor(m$log2FC_unadj, m$log2FC_adj, use = "complete.obs"), 3),
    n_regions_tested = nrow(m),
    stringsAsFactors = FALSE
  )
}

jobs <- Filter(Negate(is.null), lapply(adj_dirs, parse_adj_dir))
res <- do.call(rbind, lapply(jobs, function(j) {
  tryCatch(summarise_one(j),
    error = function(e) {
      message("  FAILED: ", conditionMessage(e))
      NULL
    }
  )
}))

if (is.null(res) || nrow(res) == 0) {
  stop("None of the adjusted runs could be paired with an unadjusted table; ",
    "nothing to plot. See the messages above for which step failed.")
}

print(res)
write.csv(res, file.path(out_dir, "confounder_adjusted_DAR_summary.csv"),
  row.names = FALSE)

## ---- 4. Figure 1: DAR counts, unadjusted vs adjusted ------------------------
plot_df <- res %>%
  dplyr::select(comparison, cell, DARs_unadjusted, DARs_adjusted) %>%
  tidyr::pivot_longer(c(DARs_unadjusted, DARs_adjusted),
    names_to = "model", values_to = "n_DARs"
  ) %>%
  dplyr::mutate(
    model = factor(ifelse(model == "DARs_unadjusted", "Unadjusted", "QC-adjusted"),
      levels = c("Unadjusted", "QC-adjusted")),
    label = paste0(cell, "\n", comparison)
  )

p_counts <- ggplot(plot_df, aes(x = label, y = n_DARs, fill = model)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65, colour = "grey20") +
  geom_text(aes(label = n_DARs),
    position = position_dodge(width = 0.7), vjust = -0.3, size = 3
  ) +
  scale_fill_manual(values = c("Unadjusted" = "#3C5488", "QC-adjusted" = "#E64B35")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Differential accessibility with and without technical-covariate adjustment",
    subtitle = paste0("Adjusted model adds ", paste(qc_adj_cols, collapse = " + "),
      " to the differential design; DAR = |log2FC| > 0.5 and padj < 0.05"),
    x = NULL, y = "Number of DARs", fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", plot.subtitle = element_text(size = 9))

ggsave(file.path(out_dir, "confounder_adjusted_DARs.pdf"), p_counts,
  width = 7.5, height = 5)

## ---- 5. Figure 2: per-region effect sizes, adjusted vs unadjusted -----------
# The count bar chart alone does not show that the SAME regions move the same
# way. This panel does: every tested region, unadjusted log2FC on x against
# adjusted log2FC on y, with the Pearson correlation annotated per comparison.
merged_df <- dplyr::bind_rows(merged_list) %>%
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
  facet_wrap(~panel, scales = "free") +
  scale_colour_manual(values = c(
    "not significant" = "grey80", "DAR unadjusted only" = "#3C5488",
    "DAR adjusted only" = "#E64B35", "DAR in both" = "#4A2377"
  )) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(
    title = "Per-region effect sizes are unchanged by technical-covariate adjustment",
    subtitle = "Dashed line, y = x. Each point is one tested region.",
    x = "log2 fold change, unadjusted model",
    y = "log2 fold change, QC-adjusted model", colour = NULL
  ) +
  theme_classic(base_size = 11) +
  theme(legend.position = "top", plot.subtitle = element_text(size = 9))

ggsave(file.path(out_dir, "confounder_adjusted_DAR_l2fc_scatter.pdf"), p_scatter,
  width = 9, height = 7)

# Per-region table for every region called a DAR by either model, so each number
# in the two figures can be traced back. Regions significant in neither model are
# left out; the full merge is tens of millions of rows across the comparisons.
dar_regions <- merged_df[merged_df$isDiff_unadj | merged_df$isDiff_adj,
  c("cell", "comparison", "id", "log2FC_unadj", "padj_unadj", "isDiff_unadj",
    "log2FC_adj", "padj_adj", "isDiff_adj", "status")]
write.csv(dar_regions,
  file.path(out_dir, "confounder_adjusted_DAR_regions.csv"), row.names = FALSE)
message("Wrote ", nrow(dar_regions), " region rows (DAR in either model)")

message("Done. Figures and tables in ", out_dir)

#####################################################################
