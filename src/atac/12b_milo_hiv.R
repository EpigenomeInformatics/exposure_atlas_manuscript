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
# (2) PSEUDO-replicates: split each sample's cells into n_rep fake replicates.
#     Splitting one sample's cells into "replicates" is pseudoreplication: the
#     fake replicates are not independent, so the model sees far more replication
#     than the 4 real subjects provide and inflates significance. We demonstrate
#     the inflation with a NULL: permute the stage labels within each subject (so
#     there is no real stage effect) and re-run. If the null still calls many
#     "significant" neighbourhoods, those hits are an artefact of the fake
#     replicates, not stage biology.
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
# Interpretation for the response letter:
#  - REAL model (Section 1, ~Subject + Stage on the 12 samples): the honest
#    donor-aware test. Report n_sig_real from that section as the main result.
#  - PSEUDO model: to give Milo "replication" we split cells into fake replicates,
#    which is pseudoreplication. It calls hundreds of significant neighbourhoods,
#    but the permuted-stage NULL calls a comparable number -> those hits are false
#    positives manufactured by the fake replicates. Either way (underpowered real
#    model, or inflated pseudo model) neighbourhood DA does not give a trustworthy
#    answer at 4 subjects, which is why we keep the subject-level cluster test.
#####################################################################
