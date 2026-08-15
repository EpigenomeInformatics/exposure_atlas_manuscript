#!/usr/bin/env Rscript

#####################################################################
# 10_3_tex_leave_one_donor_out.R
# created on 2026-08-15 by Irem B. Gunduz
# Leave-one-donor-out check on the Tex vs other CD8 signature
#
# The Tex contrast in 10_tcells.R is run on cells, so with 4 donors the
# signature could be carried by one person. Recompute it 4 times, dropping one
# donor each time, and compare log2FCs, DAR recovery and motif enrichment.
#
# Peaks are not re-called per fold: subsetting cells keeps the project peak set,
# so all folds are scored on the same features and the log2FCs line up.
#####################################################################

suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
})
set.seed(12)
addArchRThreads(threads = 8)

repo_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript"
fig_dir <- paste0(repo_dir, "/figures/")
outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/"
hiv_rds <- file.path(outputDir, "hiv_tcell_project.rds")

if (!file.exists(hiv_rds)) {
  stop("HIV T-cell project not found at ", hiv_rds,
    "\nRun src/atac/10_tcells.R first to build and cache it.")
}

# same cutoffs as the plotMarkers volcano in 10_tcells.R
fdr_cut <- 0.1
l2fc_cut <- 1

motif_fdr_cut <- 0.1
motif_l2fc_cut <- 0.5
n_top_motifs <- 25 # tracked across folds

# subject assignment and palette as in 10_tcells.R
sample_to_subject <- c(
  "hiv6_fragments.tsv.gz" = "sub1", "hiv12_fragments.tsv.gz" = "sub1", "hiv9_fragments.tsv.gz" = "sub1",
  "hiv8_fragments.tsv.gz" = "sub2", "hiv4_fragments.tsv.gz" = "sub2", "hiv1_fragments.tsv.gz" = "sub2",
  "hiv2_fragments.tsv.gz" = "sub3", "hiv7_fragments.tsv.gz" = "sub3", "hiv3_fragments.tsv.gz" = "sub3",
  "hiv11_fragments.tsv.gz" = "sub4", "hiv10_fragments.tsv.gz" = "sub4", "hiv5_fragments.tsv.gz" = "sub4"
)
colorPalette <- c(
  "sub1" = "#6A1B9A", "sub2" = "#EA80FC",
  "sub3" = "#EF7E23", "sub4" = "#ffbb00"
)
subject_label <- c(sub1 = "Subject 1", sub2 = "Subject 2",
  sub3 = "Subject 3", sub4 = "Subject 4")

## ---- 1. Project and cell annotation -----------------------------------------
message("Loading cached HIV T-cell project: ", hiv_rds)
project <- readRDS(hiv_rds)

project$Status <- ifelse(project$ClustersHIV %in% c("C1", "C2"), "Tex", "Other")
cell_sample <- getCellColData(project, select = "Sample", drop = TRUE)
project$Subject <- unname(sample_to_subject[as.character(cell_sample)])

if (anyNA(project$Subject)) {
  stop(sum(is.na(project$Subject)), " cell(s) could not be assigned to a donor; ",
    "check sample_to_subject against the Sample values in the project.")
}

comp <- table(project$Subject, project$Status)
message("Cells per donor and status:")
print(comp)
donors <- sort(unique(project$Subject))

## ---- 2. Tex vs Other, same call as 10_tcells.R -----------------------------
tex_contrast <- function(proj, label) {
  n_tex <- sum(proj$Status == "Tex")
  n_oth <- sum(proj$Status == "Other")
  message("  ", label, ": ", n_tex, " Tex vs ", n_oth, " Other cells")
  if (n_tex < 50 || n_oth < 50) {
    message("  too few cells in one group, skipping")
    return(NULL)
  }
  getMarkerFeatures(
    ArchRProj = proj,
    useMatrix = "PeakMatrix",
    groupBy = "Status",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    useGroups = c("Tex"),
    bgdGroups = "Other"
  )
}

# match folds feature by feature
peak_ids <- function(se) {
  rd <- SummarizedExperiment::rowData(se)
  paste0(rd$seqnames, ":", rd$start, "-", rd$end)
}

marker_df <- function(se) {
  data.frame(
    peak = peak_ids(se),
    Log2FC = SummarizedExperiment::assays(se)[["Log2FC"]][, 1],
    FDR = SummarizedExperiment::assays(se)[["FDR"]][, 1],
    stringsAsFactors = FALSE
  )
}

motif_df <- function(se, proj, direction = c("up", "down")) {
  direction <- match.arg(direction)
  cut <- if (direction == "up") {
    paste0("FDR <= ", motif_fdr_cut, " & Log2FC >= ", motif_l2fc_cut)
  } else {
    paste0("FDR <= ", motif_fdr_cut, " & Log2FC <= -", motif_l2fc_cut)
  }
  enr <- peakAnnoEnrichment(seMarker = se, ArchRProj = proj,
    peakAnnotation = "Motif", cutOff = cut)
  data.frame(
    TF = rownames(enr),
    mlog10Padj = SummarizedExperiment::assay(enr)[, 1],
    direction = direction,
    stringsAsFactors = FALSE
  )
}

## ---- 3. Full-cohort reference ----------------------------------------------
message("\n=== Full cohort (reference)")
se_full <- tex_contrast(project, "all four donors")
if (is.null(se_full)) stop("The full-cohort contrast failed; nothing to compare against.")
full <- marker_df(se_full)
full$isDAR <- full$FDR <= fdr_cut & abs(full$Log2FC) >= l2fc_cut
message("  ", sum(full$isDAR), " DARs at FDR <= ", fdr_cut,
  " & |log2FC| >= ", l2fc_cut,
  " (", sum(full$isDAR & full$Log2FC > 0), " up, ",
  sum(full$isDAR & full$Log2FC < 0), " down)")

motifs_full <- dplyr::bind_rows(
  motif_df(se_full, project, "up"),
  motif_df(se_full, project, "down")
)

## ---- 4. Leave-one-donor-out -------------------------------------------------
fold_markers <- list()
fold_motifs <- list()

for (d in donors) {
  message("\n=== Holding out ", subject_label[[d]])
  keep <- project$cellNames[project$Subject != d]
  proj_d <- ArchR::subsetCells(project, cellNames = keep)
  se_d <- tryCatch(tex_contrast(proj_d, paste("without", subject_label[[d]])),
    error = function(e) {
      message("  FAILED: ", conditionMessage(e))
      NULL
    })
  if (is.null(se_d)) next

  md <- marker_df(se_d)
  md$held_out <- d
  fold_markers[[d]] <- md

  mo <- tryCatch(
    dplyr::bind_rows(motif_df(se_d, proj_d, "up"), motif_df(se_d, proj_d, "down")),
    error = function(e) {
      message("  motif enrichment failed: ", conditionMessage(e))
      NULL
    })
  if (!is.null(mo)) {
    mo$held_out <- d
    fold_motifs[[d]] <- mo
  }
  rm(proj_d, se_d); gc()
}

if (length(fold_markers) == 0) stop("No leave-one-out fold completed.")

## ---- 5. Agreement with the full analysis ------------------------------------
folds <- dplyr::bind_rows(fold_markers) %>%
  dplyr::inner_join(
    full %>% dplyr::select(peak, Log2FC_full = Log2FC, FDR_full = FDR, isDAR_full = isDAR),
    by = "peak"
  ) %>%
  dplyr::mutate(isDAR_fold = FDR <= fdr_cut & abs(Log2FC) >= l2fc_cut)

# Computed per fold in a loop rather than group_by + summarise: the summarise
# version subsets inside the call (Log2FC[isDAR_full]) which dplyr 1.1 treats as
# a multi-row return, and the columns silently go missing.
fold_stats <- function(d) {
  x <- folds[folds$held_out == d, , drop = FALSE]
  dar <- x$isDAR_full
  ok <- sum(dar) > 2
  data.frame(
    held_out = unname(subject_label[[d]]),
    peaks_tested = nrow(x),
    DARs_full = sum(dar),
    DARs_fold = sum(x$isDAR_fold),
    recovered = sum(dar & x$isDAR_fold),
    recovered_pct = if (sum(dar) > 0) {
      round(100 * sum(dar & x$isDAR_fold) / sum(dar), 1)
    } else NA_real_,
    sign_concordance_pct = if (sum(dar) > 0) {
      round(100 * mean(sign(x$Log2FC[dar]) == sign(x$Log2FC_full[dar])), 1)
    } else NA_real_,
    r_all_peaks = round(stats::cor(x$Log2FC, x$Log2FC_full, use = "complete.obs"), 3),
    r_full_DARs = if (ok) {
      round(stats::cor(x$Log2FC[dar], x$Log2FC_full[dar], use = "complete.obs"), 3)
    } else NA_real_,
    stringsAsFactors = FALSE
  )
}
summary_tbl <- do.call(rbind, lapply(sort(unique(folds$held_out)), fold_stats))

message("\nAgreement with the full-cohort contrast:")
print(as.data.frame(summary_tbl))
write.csv(summary_tbl, paste0(fig_dir, "tex_leave_one_donor_out_summary.csv"),
  row.names = FALSE)

# recovered in every fold
core <- folds %>%
  dplyr::filter(isDAR_full) %>%
  dplyr::group_by(peak) %>%
  dplyr::summarise(n_folds_recovered = sum(isDAR_fold), .groups = "drop")
n_all <- sum(core$n_folds_recovered == length(fold_markers))
message(n_all, " of ", nrow(core), " full-analysis DARs (",
  round(100 * n_all / max(1, nrow(core)), 1),
  "%) are recovered in every leave-one-donor-out fold")
write.csv(core, paste0(fig_dir, "tex_leave_one_donor_out_dar_recovery.csv"),
  row.names = FALSE)

## ---- 6. Panel A: fold vs full log2 fold change ------------------------------
# thin the non-DAR background for file size; all DARs kept, no number affected
set.seed(12)
bg <- folds %>% dplyr::filter(!isDAR_full) %>%
  dplyr::group_by(held_out) %>% dplyr::slice_sample(n = 20000) %>% dplyr::ungroup()
plot_pts <- dplyr::bind_rows(bg, folds %>% dplyr::filter(isDAR_full)) %>%
  dplyr::mutate(
    panel = unname(subject_label[held_out]),
    class = ifelse(isDAR_full, "DAR in full analysis", "other tested peak")
  )
lab_a <- summary_tbl
lab_a$panel <- lab_a$held_out
lab_a$lab <- paste0("r = ", sprintf("%.3f", lab_a$r_full_DARs), " (DARs)\n",
  lab_a$recovered_pct, "% recovered")

p_a <- ggplot(plot_pts, aes(x = Log2FC_full, y = Log2FC, colour = class)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey55") +
  geom_point(size = 0.5, alpha = 0.5) +
  geom_text(data = lab_a, inherit.aes = FALSE,
    aes(x = -Inf, y = Inf, label = lab), hjust = -0.1, vjust = 1.2, size = 3) +
  facet_wrap(~panel, nrow = 1) +
  scale_colour_manual(values = c("DAR in full analysis" = "#B2182B",
    "other tested peak" = "grey80")) +
  guides(colour = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  labs(
    title = "Tex signature with each donor held out",
    subtitle = paste0("log2 fold change (Tex vs other CD8+ clusters) recomputed without the named donor, against the full-cohort value.\n",
      "Dashed line, y = x. Peaks are not re-called per fold, so every point is the same genomic region in both analyses."),
    x = "log2FC, all four donors", y = "log2FC, donor held out", colour = NULL
  ) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", plot.subtitle = element_text(size = 7),
    strip.background = element_blank(), strip.text = element_text(face = "bold"))

# response-letter material, not the supplementary figure
ggsave(p_a, file = paste0(fig_dir, "tex_leave_one_donor_out_l2fc.pdf"),
  width = 10, height = 3.2)

## ---- 7. Panel B: motif enrichment across folds -----------------------------
# Figure 2G is a motif claim, so this is the panel that has to hold. Track the
# top motifs from the full analysis through each fold.
if (length(fold_motifs) > 0) {
  top_up <- motifs_full %>%
    dplyr::filter(direction == "up") %>%
    dplyr::arrange(dplyr::desc(mlog10Padj)) %>%
    utils::head(n_top_motifs) %>%
    dplyr::pull(TF)

  mo_all <- dplyr::bind_rows(fold_motifs) %>%
    dplyr::filter(direction == "up", TF %in% top_up) %>%
    dplyr::mutate(panel = unname(subject_label[held_out]))
  mo_full <- motifs_full %>%
    dplyr::filter(direction == "up", TF %in% top_up) %>%
    dplyr::mutate(panel = "All four donors")

  mo_plot <- dplyr::bind_rows(
    mo_all %>% dplyr::select(TF, mlog10Padj, panel),
    mo_full %>% dplyr::select(TF, mlog10Padj, panel)
  )
  mo_plot$TF <- factor(mo_plot$TF, levels = rev(top_up))
  mo_plot$panel <- factor(mo_plot$panel,
    levels = c("All four donors", unname(subject_label[donors])))
  # FOX family is the one named in the text
  mo_plot$family <- ifelse(grepl("^FOX", mo_plot$TF), "FOX family", "other")

  p_b <- ggplot(mo_plot, aes(x = mlog10Padj, y = TF, colour = family)) +
    geom_point(size = 2, alpha = 0.85) +
    facet_wrap(~panel, nrow = 1) +
    scale_colour_manual(values = c("FOX family" = "#B2182B", "other" = "grey45")) +
    labs(
      title = "Motif enrichment in Tex hyper-accessible regions, with each donor held out",
      subtitle = paste0("Top ", n_top_motifs,
        " motifs from the full-cohort analysis. A motif that keeps its enrichment in every fold is not carried by one donor."),
      x = "-log10(adjusted p), motif enrichment", y = NULL, colour = NULL
    ) +
    theme_classic(base_size = 9) +
    theme(legend.position = "top", plot.subtitle = element_text(size = 7),
      strip.background = element_blank(), strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 6))

  # full 25-motif version, superseded for the paper by the compact panel below
  ggsave(p_b, file = paste0(fig_dir, "tex_leave_one_donor_out_motifs.pdf"),
    width = 11, height = 4.5)

  write.csv(
    dplyr::bind_rows(fold_motifs) %>%
      dplyr::mutate(held_out = unname(subject_label[held_out])) %>%
      dplyr::bind_rows(motifs_full %>% dplyr::mutate(held_out = "All four donors")),
    paste0(fig_dir, "tex_leave_one_donor_out_motifs.csv"), row.names = FALSE
  )

  # FOXP specifically, since that is what Figure 2G claims
  fox_all <- dplyr::bind_rows(fold_motifs) %>%
    dplyr::filter(direction == "up", grepl("^FOXP", TF))
  if (nrow(fox_all) > 0) {
    fox_sig <- fox_all %>%
      dplyr::group_by(TF) %>%
      dplyr::summarise(folds_enriched = sum(mlog10Padj > -log10(0.05)),
        min_mlog10Padj = round(min(mlog10Padj), 2), .groups = "drop")
    message("\nFOXP motif enrichment across the ", length(fold_motifs), " folds:")
    print(as.data.frame(fox_sig))
  }
}

## ---- 8. Supplementary panel -------------------------------------------------
# One scatter, footprint comparable to the donor-composition plot. Effect sizes
# with each donor held out against the full cohort, coloured by which donor was
# dropped. Points are the DARs that carry a FOX-family motif -- the regions the
# Figure 2G enrichment is built from -- so the panel shows in one place that
# both the enrichment and the regions behind it survive dropping any donor.
# Everything else is grey background.
fox_peaks <- tryCatch({
  mm <- ArchR::getMatches(project, "Motif")
  fox_cols <- grep("^FOX", colnames(mm))
  if (length(fox_cols) == 0) stop("no FOX motifs in the annotation")
  hit <- Matrix::rowSums(SummarizedExperiment::assay(mm)[, fox_cols, drop = FALSE]) > 0
  gr <- SummarizedExperiment::rowRanges(mm)
  data.frame(
    peak = paste0(as.character(GenomicRanges::seqnames(gr)), ":",
      GenomicRanges::start(gr), "-", GenomicRanges::end(gr)),
    has_FOX = as.logical(hit), stringsAsFactors = FALSE
  )
}, error = function(e) {
  message("Could not read motif matches (", conditionMessage(e),
    "); the panel will show all DARs instead.")
  NULL
})

panel_dat <- folds %>% dplyr::filter(isDAR_full)
if (!is.null(fox_peaks)) {
  panel_dat <- dplyr::left_join(panel_dat, fox_peaks, by = "peak")
  panel_dat$has_FOX[is.na(panel_dat$has_FOX)] <- FALSE
} else {
  panel_dat$has_FOX <- TRUE
}
panel_dat$Subject <- unname(subject_label[panel_dat$held_out])

n_fox <- dplyr::n_distinct(panel_dat$peak[panel_dat$has_FOX])
message("FOX-motif DARs: ", n_fox, " of ", dplyr::n_distinct(panel_dat$peak))

# per-donor correlation on the FOX-motif DARs, for the legend labels
fox_r <- panel_dat %>%
  dplyr::filter(has_FOX) %>%
  dplyr::group_by(Subject) %>%
  dplyr::summarise(r = stats::cor(Log2FC, Log2FC_full, use = "complete.obs"),
    .groups = "drop") %>%
  dplyr::mutate(lab = paste0(Subject, "  (r = ", sprintf("%.2f", r), ")"))

lab_map <- stats::setNames(fox_r$lab, fox_r$Subject)
panel_dat$SubjectLab <- unname(lab_map[panel_dat$Subject])
pal <- stats::setNames(unname(colorPalette[donors]),
  unname(lab_map[unname(subject_label[donors])]))

p_panel <- ggplot() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
    colour = "grey50", linewidth = 0.4) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(data = dplyr::filter(panel_dat, !has_FOX),
    aes(x = Log2FC_full, y = Log2FC), colour = "grey85", size = 0.4, alpha = 0.4) +
  geom_point(data = dplyr::filter(panel_dat, has_FOX),
    aes(x = Log2FC_full, y = Log2FC, colour = SubjectLab),
    size = 0.7, alpha = 0.65) +
  scale_colour_manual(values = pal, name = "donor held out") +
  guides(colour = guide_legend(override.aes = list(size = 2.5, alpha = 1), ncol = 1)) +
  labs(
    x = "log2FC, all four donors",
    y = "log2FC, one donor held out"
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "right",
    legend.key.height = unit(0.38, "cm"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7)
  )

ggsave(p_panel, file = paste0(fig_dir, "tex_leave_one_donor_out_panel.pdf"),
  width = 5.4, height = 3.4)
message("Supplementary panel: tex_leave_one_donor_out_panel.pdf")

# numbers for the legend, so they are quoted rather than recalled
message("\nFor the legend:")
message("  DARs (full cohort): ", sum(full$isDAR), "; carrying a FOX motif: ", n_fox)
message("  FOX-motif DAR log2FC r per fold: ",
  paste0(sprintf("%.2f", fox_r$r), collapse = ", "))
message("  recovered per fold: ",
  paste0(summary_tbl$recovered_pct, "%", collapse = ", "))
message("  recovered in all folds: ", n_all, " (",
  round(100 * n_all / max(1, nrow(core)), 1), "%)")
message("  sign concordance: ",
  paste0(range(summary_tbl$sign_concordance_pct), collapse = "-"), "%")
message("  FOXP motifs enriched in all ", length(fold_motifs), " folds (see CSV)")

message("\nDone. Figures and tables in ", fig_dir)

#####################################################################
