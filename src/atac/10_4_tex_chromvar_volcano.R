#!/usr/bin/env Rscript

#####################################################################
# 10_4_tex_chromvar_volcano.R
# created on 2026-08-15 by Irem B. Gunduz
# chromVAR z-score volcano, Tex vs other CD8+ clusters
#
# CTLA4, HAVCR2 and TIGIT all rise on ACTIVATED CD8 cells, so marker activity
# does not separate exhaustion from activation. Every motif is tested, so the
# exhaustion regulators have to come out of the data rather than be asserted.
#
# Run on both motif sets: cisbp (named TFs, for TOX / NR4A / TCF7) and the
# Altius archetypes (same style as the zdiff volcano elsewhere in the paper).
#
# !! The Wilcoxon runs on CELLS, which are not independent within a donor. With
# ~7000 cells the p-value largely reports the cell count, so a motif is only
# called differential if it ALSO moves the same way in all four donors. Both
# criteria are in the figure caption.
#####################################################################

suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})
set.seed(12)
addArchRThreads(threads = 8)

repo_dir <- "/scratch/icbb/igunduz/irem_github/exposure_atlas_manuscript"
fig_dir <- paste0(repo_dir, "/figures/")
outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/"
hiv_rds <- file.path(outputDir, "hiv_tcell_project.rds")

if (!file.exists(hiv_rds)) {
  stop("HIV T-cell project not found at ", hiv_rds,
    "\nRun src/atac/10_tcells.R first to build and cache it.")
}

tex_clusters <- c("C1", "C2")

# matrices to run; any not present in the project is skipped
matrices <- c("MotifMatrix", "altiusMatrix")

zdiff_cut <- 0.25   # |difference in mean chromVAR z| to call differential
padj_cut <- 0.01    # threshold used to CALL motifs differential

# The red line marks padj_cut and therefore sits low in the panel. With ~7000
# cells the weakest motif we call is at -log10(padj) = 6.6 (cisbp) / 11.6
# (Altius), so padj does not bind at any conventional threshold; |zdiff| and
# donor consistency do the selecting. Keep the line reporting the real cutoff.
padj_line <- padj_cut
padj_cap <- 50      # -log10 cap; cell-level p underflows to 0
n_label_each <- 10  # top motifs labelled per direction
x_lim <- c(-2, 2)   # clipped, not filtered; anything outside is reported below

# always labelled if present, whatever their rank; motif names differ between
# the two sets (cisbp keeps gene symbols, Altius collapses families into
# lowercase archetypes), so one list per matrix
tf_always <- list(
  MotifMatrix = c("TOX", "TOX2", "NR4A1", "NR4A2", "NR4A3", "TCF7", "TCF7L2",
    "LEF1", "TBX21", "EOMES", "FOXP1", "FOXP2", "FOXP3", "FOXP4", "FOXO1",
    "FOXK1", "RUNX1", "BATF", "PRDM1"),
  altiusMatrix = c("fox_4", "fox_6", "tcf_lef", "lef1", "runx_1", "runx_2",
    "batf", "ap1_1", "ap1_2", "nfat_1", "nfat_2", "tbx_1", "nr_11", "prdm1")
)

sample_to_subject <- c(
  "hiv6_fragments.tsv.gz" = "sub1", "hiv12_fragments.tsv.gz" = "sub1", "hiv9_fragments.tsv.gz" = "sub1",
  "hiv8_fragments.tsv.gz" = "sub2", "hiv4_fragments.tsv.gz" = "sub2", "hiv1_fragments.tsv.gz" = "sub2",
  "hiv2_fragments.tsv.gz" = "sub3", "hiv7_fragments.tsv.gz" = "sub3", "hiv3_fragments.tsv.gz" = "sub3",
  "hiv11_fragments.tsv.gz" = "sub4", "hiv10_fragments.tsv.gz" = "sub4", "hiv5_fragments.tsv.gz" = "sub4"
)

## ---- 1. Project and cell annotation -----------------------------------------
message("Loading cached HIV T-cell project: ", hiv_rds)
project <- readRDS(hiv_rds)

avail <- getAvailableMatrices(project)
message("Matrices in the project: ", paste(avail, collapse = ", "))
use_matrices <- intersect(matrices, avail)
if (length(use_matrices) == 0) {
  stop("None of ", paste(matrices, collapse = ", "), " are in the project.")
}

## ---- 2. Pull the chromVAR z-scores from one matrix ---------------------------
# The z-scores and raw deviations live either in two ASSAYS ("z" / "deviations")
# or in two sets of ROWS marked by rowData$seqnames, depending on ArchR version.
# Try the assay route first and report what was found on failure.
get_z <- function(mat_name) {
  mm <- getMatrixFromProject(project, useMatrix = mat_name)
  rd <- SummarizedExperiment::rowData(mm)
  an <- SummarizedExperiment::assayNames(mm)
  message("  assays: ", paste(an, collapse = ", "),
    " | rowData: ", paste(colnames(rd), collapse = ", "))

  nm <- if ("name" %in% colnames(rd)) as.character(rd$name) else rownames(mm)

  if (!is.null(an) && "z" %in% an) {
    z <- SummarizedExperiment::assay(mm, "z")
    rownames(z) <- sub("^[^:]*:", "", nm)
  } else {
    kind <- if ("seqnames" %in% colnames(rd)) as.character(rd$seqnames) else NULL
    if ((is.null(kind) || !any(kind == "z")) && !is.null(rownames(mm)) &&
        any(grepl("^z:", rownames(mm)))) {
      kind <- ifelse(grepl("^z:", rownames(mm)), "z", "deviations")
      nm <- sub("^[^:]*:", "", rownames(mm))
    }
    if (is.null(kind) || !any(kind == "z")) {
      stop("Could not locate z-scores in ", mat_name,
        ".\n  assays: ", paste(an, collapse = ", "),
        "\n  rowData: ", paste(colnames(rd), collapse = ", "),
        "\n  first names: ", paste(utils::head(nm, 3), collapse = ", "))
    }
    z <- SummarizedExperiment::assay(mm)[kind == "z", , drop = FALSE]
    rownames(z) <- nm[kind == "z"]
  }
  list(z = z, cells = colnames(mm))
}

## ---- 3. Test every motif ----------------------------------------------------
test_matrix <- function(zmat, status, subj, donors) {
  is_tex <- status == "Tex"
  res <- do.call(rbind, lapply(seq_len(nrow(zmat)), function(i) {
    z <- zmat[i, ]
    d <- mean(z[is_tex], na.rm = TRUE) - mean(z[!is_tex], na.rm = TRUE)
    p <- tryCatch(stats::wilcox.test(z[is_tex], z[!is_tex])$p.value,
      error = function(e) NA_real_)
    same <- vapply(donors, function(dd) {
      k <- subj == dd
      isTRUE(sign(mean(z[k & is_tex], na.rm = TRUE) -
        mean(z[k & !is_tex], na.rm = TRUE)) == sign(d))
    }, logical(1))
    c(zdiff = d, p = p, donors_consistent = sum(same))
  }))
  data.frame(
    motif = rownames(zmat),
    zdiff = res[, "zdiff"],
    p = res[, "p"],
    donors_consistent = as.integer(res[, "donors_consistent"]),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      padj = stats::p.adjust(p, method = "BH"),
      # cell-level p underflows to 0; floor before -log10 or the strongest
      # motifs become Inf and drop off the plot
      mlog10padj = pmin(-log10(pmax(padj, .Machine$double.xmin)), padj_cap),
      robust = donors_consistent == length(donors),
      group = dplyr::case_when(
        robust & padj < padj_cut & zdiff > zdiff_cut ~ "Tex",
        robust & padj < padj_cut & zdiff < -zdiff_cut ~ "Other",
        TRUE ~ "NO"
      )
    ) %>%
    dplyr::arrange(dplyr::desc(abs(zdiff)))
}

## ---- 4. Volcano, same style as the zdiff volcanoes elsewhere ----------------
volcano_plot <- function(v, mat_name) {
  keep_lab <- tf_always[[mat_name]]
  if (is.null(keep_lab)) keep_lab <- character(0)
  lab <- dplyr::bind_rows(
    v %>% dplyr::filter(group != "NO") %>%
      dplyr::group_by(group) %>%
      dplyr::slice_max(order_by = abs(zdiff), n = n_label_each, with_ties = FALSE) %>%
      dplyr::ungroup(),
    v %>% dplyr::filter(motif %in% keep_lab)
  ) %>% dplyr::distinct()

  ggplot(v, aes(x = zdiff, y = mlog10padj, colour = group)) +
    geom_vline(xintercept = c(-zdiff_cut, zdiff_cut),
      colour = "grey40", linewidth = 0.3) +
    geom_hline(yintercept = -log10(padj_line),
      colour = "#B22222", linewidth = 0.3, linetype = "dashed") +
    annotate("text", x = x_lim[2], y = -log10(padj_line),
      label = paste0("padj = ", padj_line), colour = "#B22222",
      size = 3.2, hjust = 1, vjust = -0.6) +
    geom_point(size = 1.4, alpha = 0.85) +
    ggrepel::geom_text_repel(data = lab, aes(label = motif), size = 3.6,
      max.overlaps = Inf, min.segment.length = Inf, segment.size = 0,
      # no leader lines, so keep the repulsion weak or a label drifts away
      # from the point it belongs to
      box.padding = 0.4, point.padding = 0.2, force = 0.6,
      seed = 12, show.legend = FALSE) +
    scale_colour_manual(
      values = c("Tex" = "#B22222", "NO" = "grey75", "Other" = "#3C6EB4"),
      breaks = c("Tex", "NO", "Other"),
      labels = c("Tex", "NO", "other clusters"),
      name = NULL) +
    coord_cartesian(xlim = x_lim) +
    scale_x_continuous(breaks = seq(x_lim[1], x_lim[2], by = 0.5)) +
    labs(x = "chromVAR z-score difference (Tex - other CD8+ clusters)",
      y = expression(-log[10] * "(padj)"),
      title = mat_name) +
    theme_classic(base_size = 11) +
    theme(plot.title = element_text(size = 10, face = "bold"),
      legend.position = "right")
}

## ---- 5. Run ------------------------------------------------------------------
for (mat_name in use_matrices) {
  message("\n=== ", mat_name)
  got <- get_z(mat_name)
  zmat <- got$z

  ccd <- getCellColData(project)
  idx <- match(got$cells, rownames(ccd))
  clust <- as.character(ccd$ClustersHIV)[idx]
  samp <- as.character(ccd$Sample)[idx]
  subj <- unname(sample_to_subject[samp])
  ok <- !is.na(clust) & !is.na(subj)

  zmat <- zmat[, ok, drop = FALSE]
  status <- ifelse(clust[ok] %in% tex_clusters, "Tex", "Other")
  subj <- subj[ok]
  donors <- sort(unique(subj))
  message("  ", nrow(zmat), " motifs | ", sum(status == "Tex"), " Tex vs ",
    sum(status == "Other"), " other cells | ", length(donors), " donors")

  v <- test_matrix(zmat, status, subj, donors)
  tag <- sub("Matrix$", "", mat_name)

  write.csv(v %>% dplyr::mutate(dplyr::across(where(is.numeric), ~ signif(.x, 4))),
    paste0(fig_dir, "tex_chromvar_volcano_", tag, ".csv"), row.names = FALSE)

  ggsave(volcano_plot(v, mat_name),
    file = paste0(fig_dir, "tex_chromvar_volcano_", tag, ".pdf"),
    width = 6.4, height = 4.6)

  off <- v %>% dplyr::filter(zdiff < x_lim[1] | zdiff > x_lim[2])
  if (nrow(off) > 0) {
    message("  outside the x range and therefore not drawn: ",
      paste0(off$motif, " (", round(off$zdiff, 2), ")", collapse = ", "))
  }

  message("  differential: ", sum(v$group == "Tex"), " up in Tex, ",
    sum(v$group == "Other"), " up in other clusters, of ", nrow(v),
    " (|zdiff| > ", zdiff_cut, ", padj < ", padj_cut,
    ", consistent in all ", length(donors), " donors)")
  message("  top by |zdiff| among the differential ones:")
  print(as.data.frame(v %>% dplyr::filter(group != "NO") %>%
    dplyr::transmute(motif, zdiff = round(zdiff, 3), padj = signif(padj, 3),
      donors_consistent, group) %>% utils::head(12)))

  hits <- v %>% dplyr::filter(motif %in% tf_always[[mat_name]])
  if (nrow(hits) > 0) {
    message("  exhaustion vs activation discriminators:")
    print(as.data.frame(hits %>%
      dplyr::transmute(motif, zdiff = round(zdiff, 3), padj = signif(padj, 3),
        donors_consistent, group) %>% dplyr::arrange(dplyr::desc(zdiff))))
  }
  rm(zmat, v); gc()
}

message("\nDone. Volcanoes and tables in ", fig_dir)

#####################################################################
