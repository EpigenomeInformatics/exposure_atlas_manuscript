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
