# load libraries
library(tidyr)
library(dplyr)
library(pheatmap)

jaccard_helper <- function(mat, i, j) {
  AiB <- mat[i, j]
  AuB <- sum(mat[i, ]) + sum(mat[, j]) - AiB
  AiB / AuB
}


Jaccard <- function(data) {
  data <- t(data)
  df <- data # to assign names to later
  jaccards <- matrix(data = NA, nrow = NROW(data), ncol = ncol(data))
  for (r in 1:NROW(data)) {
    for (c in 1:ncol(data)) {
      jaccards[r, c] <- jaccard_helper(data, r, c)
    }
  }

  colnames(jaccards) <- colnames(data)
  row.names(jaccards) <- row.names(data)
  pheatmap <- pheatmap::pheatmap(jaccards,
    show_rownames = T,
    cluster_rows = F,
    display_numbers = as.matrix(df)
  )
  return(pheatmap)
}

