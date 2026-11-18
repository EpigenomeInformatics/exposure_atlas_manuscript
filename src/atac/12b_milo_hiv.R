#!/usr/bin/env Rscript

#####################################################################
# 12b_milo_hiv.R
# miloR neighbourhood differential-abundance on the HIV CD8+ T-cell subset
#####################################################################

set.seed(1)
suppressPackageStartupMessages({
  library(ArchR)
  library(miloR)
  library(SingleCellExperiment)
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
})

fig_dir   <- "/scratch/icbb/igunduz/irem_github/exposure_atlas_manuscript/figures/"
outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/"
hiv_rds   <- file.path(outputDir, "hiv_tcell_project.rds")
stopifnot(file.exists(hiv_rds))
project <- readRDS(hiv_rds)

# --- sample -> subject / stage maps (same as 10_tcells.R) ---------------------
sample_to_subject <- c(
  "hiv6_fragments.tsv.gz"  = "sub1", "hiv12_fragments.tsv.gz" = "sub1", "hiv9_fragments.tsv.gz"  = "sub1",
  "hiv8_fragments.tsv.gz"  = "sub2", "hiv4_fragments.tsv.gz"  = "sub2", "hiv1_fragments.tsv.gz"  = "sub2",
  "hiv2_fragments.tsv.gz"  = "sub3", "hiv7_fragments.tsv.gz"  = "sub3", "hiv3_fragments.tsv.gz"  = "sub3",
  "hiv11_fragments.tsv.gz" = "sub4", "hiv10_fragments.tsv.gz" = "sub4", "hiv5_fragments.tsv.gz"  = "sub4"
)
sample_files        <- names(sample_to_subject)
time_point          <- c("pre", "acute", "chronic", "pre", "acute", "chronic",
                         "pre", "acute", "chronic", "pre", "acute", "chronic")
sample_to_timepoint <- setNames(time_point, sample_files)

# --- Build a SingleCellExperiment carrying the HIV LSI embedding ---------------
lsi   <- getReducedDims(project, reducedDims = "IterativeLSI_hiv")
cells <- rownames(lsi)
samp  <- sapply(strsplit(cells, "#"), `[`, 1)

# provide per-cell depth (nFrags) as a 1-row counts assay so Milo's depth filter is defined
nfrags <- getCellColData(project, select = "nFrags", drop = TRUE)[match(cells, project$cellNames)]
sce <- SingleCellExperiment(
  assays = list(counts = matrix(as.numeric(nfrags), nrow = 1,
                                dimnames = list("depth", cells)))
)
reducedDim(sce, "LSI") <- lsi
sce$Sample  <- samp
sce$Subject <- sample_to_subject[samp]
sce$Stage   <- factor(sample_to_timepoint[samp], levels = c("pre", "acute", "chronic"))
sce$Cluster <- as.character(project$ClustersHIV)[match(cells, project$cellNames)]

d_use <- min(30, ncol(lsi))   # never request more LSI dims than exist

# --- Build the Milo object and neighbourhoods ---------------------------------
milo <- Milo(sce)
milo <- buildGraph(milo, k = 30, d = d_use, reduced.dim = "LSI")
milo <- makeNhoods(milo, prop = 0.1, k = 30, d = d_use, refined = TRUE, reduced_dims = "LSI")

pdf(file = paste0(fig_dir, "milo_hiv_nhood_sizes.pdf"), width = 5, height = 4)
print(plotNhoodSizeHist(milo))
dev.off()

milo <- calcNhoodDistance(milo, d = d_use, reduced.dim = "LSI")

#####################################################################
# (1) REAL replicates: the 12 samples, design ~ Subject + Stage
#####################################################################
milo <- countCells(milo, meta.data = as.data.frame(colData(milo)), samples = "Sample")

design_real <- distinct(as.data.frame(colData(milo))[, c("Sample", "Subject", "Stage")])
rownames(design_real) <- design_real$Sample

res_real <- tryCatch(
  testNhoods(milo, design = ~ Subject + Stage, design.df = design_real, reduced.dim = "LSI"),
  error = function(e) { message("Real-replicate testNhoods failed: ", conditionMessage(e)); NULL }
)
if (!is.null(res_real)) {
  n_sig_real <- sum(res_real$SpatialFDR < 0.1, na.rm = TRUE)
  message(sprintf("REAL 12-sample model (~Subject + Stage): %d / %d neighbourhoods at SpatialFDR < 0.1",
                  n_sig_real, nrow(res_real)))
  # majority CD8 cluster per neighbourhood (adds Cluster, Cluster_fraction)
  res_real <- annotateNhoods(milo, res_real, coldata_col = "Cluster")
  write.csv(res_real, paste0(fig_dir, "milo_hiv_real_results.csv"), row.names = FALSE)
}

#####################################################################
# (2) PSEUDO-replicates: split each sample's cells into n_rep fake replicates
#####################################################################
run_pseudo <- function(seed, n_rep = 3, permute = FALSE) {
  set.seed(seed)
  m2 <- milo

  stage_use <- as.character(m2$Stage)
  if (permute) {
    # shuffle the 3 stages within each subject's 3 samples -> no real stage effect
    perm_map <- sample_to_timepoint
    for (s in unique(sample_to_subject)) {
      smp <- names(sample_to_subject)[sample_to_subject == s]
      perm_map[smp] <- sample(sample_to_timepoint[smp])
    }
    stage_use <- perm_map[m2$Sample]
  }
  m2$StageTest    <- factor(stage_use, levels = c("pre", "acute", "chronic"))
  m2$PseudoSample <- paste0(m2$Sample, "_r", sample(seq_len(n_rep), ncol(m2), replace = TRUE))
  m2 <- countCells(m2, meta.data = as.data.frame(colData(m2)), samples = "PseudoSample")
  dd <- distinct(as.data.frame(colData(m2))[, c("PseudoSample", "Subject", "StageTest")])
  rownames(dd) <- dd$PseudoSample
  r <- tryCatch(
    testNhoods(m2, design = ~ Subject + StageTest, design.df = dd, reduced.dim = "LSI"),
    error = function(e) {
      message("  seed ", seed, " (permute=", permute, ") failed: ", conditionMessage(e)); NULL })
  if (is.null(r)) return(NA_integer_)
  sum(r$SpatialFDR < 0.1, na.rm = TRUE)
}

seeds      <- 1:5
n_sig_real <- sapply(seeds, run_pseudo, permute = FALSE)  # real stage labels
n_sig_null <- sapply(seeds, run_pseudo, permute = TRUE)   # permuted (null) labels

pseudo_summary <- data.frame(seed = seeds,
                             n_sig_real_labels   = n_sig_real,
                             n_sig_permuted_null = n_sig_null)
write.csv(pseudo_summary, paste0(fig_dir, "milo_hiv_pseudo_seed_summary.csv"), row.names = FALSE)
print(pseudo_summary)

message(sprintf(
  "PSEUDO-replicate DA: real-label median %d vs permuted-null median %d significant neighbourhoods (SpatialFDR < 0.1).",
  as.integer(median(n_sig_real, na.rm = TRUE)), as.integer(median(n_sig_null, na.rm = TRUE))))
message("A large permuted-null count means the significance is a pseudoreplication ",
        "artefact of the fake replicates, not a real stage effect.")

#####################################################################
# Interpretation: neighbourhood DA is not trustworthy at 4 subjects; keep the subject-level cluster test
#####################################################################

#####################################################################
# (3) Response figure: A = nhood logFC by CD8 cluster, B = real vs permuted per seed
# out: figures/response_figure_3_milo.pdf / .png
#####################################################################
suppressPackageStartupMessages(library(patchwork))

real_csv   <- file.path(fig_dir, "milo_hiv_real_results.csv")
pseudo_csv <- file.path(fig_dir, "milo_hiv_pseudo_seed_summary.csv")
out_pdf    <- file.path(fig_dir, "response_figure_3_milo.pdf")
out_png    <- file.path(fig_dir, "response_figure_3_milo.png")
stopifnot(file.exists(real_csv), file.exists(pseudo_csv))

fdr_cut <- 0.1
mix_cut <- 0.7  # nhoods whose majority cluster is < 70% of cells are labelled "Mixed"

# copied from 07_7_covariate_balance_response_figure.R
theme_resp <- theme_classic(base_size = 11) +
  theme(panel.grid.major = element_line(colour = "grey85", linewidth = 0.3),
        panel.grid.minor = element_line(colour = "grey92", linewidth = 0.2),
        plot.title = element_text(face = "bold", hjust = 0, size = 13),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"))

## ---- Panel A: neighbourhood logFC by CD8 cluster --------------------------
nh <- read.csv(real_csv, stringsAsFactors = FALSE)
stopifnot(all(c("logFC", "SpatialFDR", "Cluster", "Cluster_fraction") %in% names(nh)))
nh <- nh %>%
  dplyr::filter(!is.na(SpatialFDR)) %>%
  dplyr::mutate(Cluster = ifelse(Cluster_fraction < mix_cut, "Mixed", Cluster),
                sig     = SpatialFDR < fdr_cut)

n_sig_nh <- sum(nh$sig)
message(sprintf("panel A: %d / %d neighbourhoods at SpatialFDR < %.1f (%d labelled Mixed)",
                n_sig_nh, nrow(nh), fdr_cut, sum(nh$Cluster == "Mixed")))
print(table(Cluster = nh$Cluster, significant = nh$sig))

# testNhoods() tests the last design coefficient, i.e. Stage chronic vs pre (adjusted for Subject)
jit <- position_jitter(width = 0.24, height = 0, seed = 3)

# sort so significant nhoods are drawn last (on top); red below 0.1, grey above
nh <- dplyr::arrange(nh, dplyr::desc(SpatialFDR))

pA <- ggplot(nh, aes(x = Cluster, y = logFC, colour = SpatialFDR)) +
  geom_hline(yintercept = 0, colour = "grey60") +
  geom_point(aes(size = sig), position = jit, alpha = 0.85, show.legend = c(size = FALSE)) +
  scale_size_manual(values = c(`FALSE` = 1.4, `TRUE` = 2.8)) +
  scale_colour_gradientn(colours = c("#7b1d14", "#e74c3c", "grey70", "grey88"),
                         values = c(0, fdr_cut - 1e-3, fdr_cut, 1), limits = c(0, 1),
                         breaks = c(0, fdr_cut, 0.5, 1), name = "SpatialFDR") +
  guides(colour = guide_colourbar(barheight = grid::unit(8, "lines"))) +
  labs(x = "CD8+ T-cell cluster (majority label)",
       y = "Neighbourhood log fold change\n(chronic vs pre, adjusted for subject)",
       title = "A") +
  theme_resp +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# dashed line at the |logFC| where neighbourhoods start passing SpatialFDR = 0.1
if (n_sig_nh > 0) {
  lfc_cut <- min(abs(nh$logFC[nh$sig]))
  pA <- pA +
    geom_hline(yintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", colour = "#c0392b") +
    annotate("text", x = 0.6, y = lfc_cut, label = "SpatialFDR = 0.1", colour = "#c0392b",
             fontface = "italic", size = 3.1, hjust = 0, vjust = -0.4)
}

## ---- Panel B: significant nhoods per seed, real vs permuted labels --------
ps <- read.csv(pseudo_csv, stringsAsFactors = FALSE)
stopifnot(all(c("seed", "n_sig_real_labels", "n_sig_permuted_null") %in% names(ps)))
ps_long <- dplyr::bind_rows(
  data.frame(seed = ps$seed, labels = "Real stage labels",     n_sig = ps$n_sig_real_labels),
  data.frame(seed = ps$seed, labels = "Permuted stage labels", n_sig = ps$n_sig_permuted_null)
) %>%
  dplyr::mutate(labels = factor(labels, levels = c("Real stage labels", "Permuted stage labels")))

med_real <- median(ps$n_sig_real_labels, na.rm = TRUE)
med_perm <- median(ps$n_sig_permuted_null, na.rm = TRUE)
message(sprintf("panel B: median significant nhoods, real %.0f vs permuted %.0f (%d seeds)",
                med_real, med_perm, nrow(ps)))

pB <- ggplot(ps_long, aes(x = factor(seed), y = n_sig, fill = labels)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_hline(yintercept = med_real, linetype = "dashed", colour = "#3b6ea5") +
  geom_hline(yintercept = med_perm, linetype = "dashed", colour = "#c1741f") +
  scale_fill_manual(values = c("Real stage labels" = "#3b6ea5",
                               "Permuted stage labels" = "#c1741f"), name = NULL) +
  labs(x = "Pseudo-replicate seed",
       y = sprintf("Neighbourhoods at SpatialFDR < %.1f", fdr_cut),
       title = "B") +
  theme_resp +
  theme(legend.position = "bottom")

## ---- Compose ---------------------------------------------------------------
fig <- (pA | pB) + patchwork::plot_layout(widths = c(1.6, 1))

fig_w <- 12
fig_h <- 5

ggsave(out_pdf, fig, width = fig_w, height = fig_h, useDingbats = FALSE)
ggsave(out_png, fig, width = fig_w, height = fig_h, dpi = 300, bg = "white")
message("Wrote ", out_pdf, "\n       ", out_png)
#####################################################################
