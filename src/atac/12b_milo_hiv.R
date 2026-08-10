#!/usr/bin/env Rscript

#####################################################################
# 12b_milo_hiv.R
# miloR neighbourhood differential-abundance on the HIV CD8+ T-cell subset.
#
# PURPOSE (reviewer R3.2): a reviewer asked whether neighbourhood-based
# differential-abundance testing (Milo) would give better resolution than
# predefined clusters. Our HIV design is 4 subjects x 3 stages = 12 samples,
# which is at or below Milo's minimum for a mixed model. This script shows that
# Milo does not give a usable result here. We run it two ways:
#   (1) with the real 12 samples as replicates (design ~ Subject + Stage), and
#   (2) with pseudo-replicates (each sample's cells split into fake replicates,
#       as suggested), which is the only way to give the model "replication".
# We then show that (1) yields essentially no significant neighbourhoods and
# (2) the pseudo-replicate result is unstable: the set of "significant"
# neighbourhoods changes with the random split (seed), i.e. the significance is
# an artefact of the fake replicates, not a real signal.
#
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

fig_dir   <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/figures/"
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

# Milo's buildGraph filters LSI dims that correlate with sequencing depth, using
# colSums of the counts assay. A dummy all-zero assay makes depth constant, which
# yields NA correlations and a "subscript out of bounds" error in findKNN. We give
# it a real per-cell depth (nFrags) as a 1-row counts assay so the depth filter is
# well defined (it will drop LSI dim 1, as intended for LSI embeddings).
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
  write.csv(res_real, paste0(fig_dir, "milo_hiv_real_results.csv"), row.names = FALSE)
}

#####################################################################
# (2) PSEUDO-replicates: split each sample's cells into n_rep fake replicates,
#     repeated over several seeds to show the DA calls are unstable.
#####################################################################
run_pseudo <- function(seed, n_rep = 3) {
  set.seed(seed)
  m2 <- milo
  m2$PseudoSample <- paste0(m2$Sample, "_r", sample(seq_len(n_rep), ncol(m2), replace = TRUE))
  m2 <- countCells(m2, meta.data = as.data.frame(colData(m2)), samples = "PseudoSample")
  dd <- distinct(as.data.frame(colData(m2))[, c("PseudoSample", "Subject", "Stage")])
  rownames(dd) <- dd$PseudoSample
  r <- tryCatch(
    testNhoods(m2, design = ~ Subject + Stage, design.df = dd, reduced.dim = "LSI"),
    error = function(e) { message("  seed ", seed, " failed: ", conditionMessage(e)); NULL })
  if (is.null(r)) return(NULL)
  sig <- which(r$SpatialFDR < 0.1)
  list(n_sig = length(sig), sig_nhoods = r$Nhood[sig])
}

seeds       <- 1:5
pseudo_runs <- lapply(seeds, run_pseudo)
n_sig       <- sapply(pseudo_runs, function(x) if (is.null(x)) NA_integer_ else x$n_sig)

# Jaccard overlap of the significant-neighbourhood sets between seeds (stability)
sig_sets <- lapply(pseudo_runs, function(x) if (is.null(x)) integer(0) else x$sig_nhoods)
jac <- function(a, b) if (length(union(a, b)) == 0) NA_real_ else length(intersect(a, b)) / length(union(a, b))
jac_mat <- outer(seq_along(sig_sets), seq_along(sig_sets),
                 Vectorize(function(i, j) jac(sig_sets[[i]], sig_sets[[j]])))

write.csv(data.frame(seed = seeds, n_significant_nhoods = n_sig),
          paste0(fig_dir, "milo_hiv_pseudo_seed_summary.csv"), row.names = FALSE)
write.csv(round(jac_mat, 3),
          paste0(fig_dir, "milo_hiv_pseudo_seed_jaccard.csv"), row.names = FALSE)

message("Pseudo-replicate significant-neighbourhood counts by seed: ",
        paste(n_sig, collapse = ", "))
message("Mean pairwise Jaccard overlap of the significant sets across seeds: ",
        round(mean(jac_mat[upper.tri(jac_mat)], na.rm = TRUE), 3),
        "  (low overlap = the DA calls are an artefact of the random split)")

#####################################################################
# Interpretation for the response letter:
#  - REAL model: report n_sig_real (expected ~0), i.e. no neighbourhood reaches
#    significance once donor structure is respected at n = 4 subjects.
#  - PSEUDO model: apparent hits appear, but n varies by seed and the sets barely
#    overlap (low mean Jaccard), so the significance is manufactured by the fake
#    replicates. This is why we keep the subject-level cluster test rather than a
#    neighbourhood-based DA for this cohort.
#####################################################################
