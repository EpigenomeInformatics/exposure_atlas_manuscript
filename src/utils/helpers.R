# Compute the z-score, helper for clusterATACzscores
computeZScore <- function(counts) {
  counts <- DelayedArray::DelayedArray(counts)
  counts <- (counts - DelayedMatrixStats::rowMeans2(counts)) / DelayedMatrixStats::rowSds(counts)
  counts[base::is.nan(counts)] <- 0
  return(counts)
}

# helper function for atac
cutL0.5fc2Padj05 <- function(dm) {
  abs(dm[, "log2FoldChange"]) > 0.5 & dm[, "padj"] < 0.05
}
