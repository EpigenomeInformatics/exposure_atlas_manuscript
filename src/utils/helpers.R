# Compute the z-score, helper for clusterATACzscores
computeZScore <- function(counts) {
  counts <- (counts - matrixStats::rowMeans2(counts)) / matrixStats::rowSds(counts)
  counts[base::is.nan(counts)] <- 0
  return(counts)
}

# helper function for atac
cutL0.5fc2Padj05 <- function(dm, padj = 0.05) {
  abs(dm[, "log2FoldChange"]) > 0.5 & dm[, "padj"] < padj
}

consensus_kmeans <- function(counts, k, km_repeats) {
  partition_list <- lapply(seq_len(km_repeats), function(i) {
    as.cl_hard_partition(kmeans(counts, k, iter.max = 50))
  })
  partition_list <- cl_ensemble(list = partition_list)
  partition_consensus <- cl_consensus(partition_list)
  as.vector(cl_class_ids(partition_consensus))
}
