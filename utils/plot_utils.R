## function to prepare differentially methylated regions for scatter plot
## Written by: Irem B. Gündüz & Fabian Müller
prepareDMRforPlot <- function(cell, comp, path, global_peaklist, region, changeMethod = "meandiff") {
  # organise global_peakset and set as an annotation
  colnames(global_peaklist) <- c("Chromosome", "Start", "End")
  rnb.set.annotation(type = region, regions = global_peaklist, assembly = "hg38")

  # create analysis dir
  analysis.dir <- paste0(path, cell, "/", comp)
  # load the diffmeth object
  diffMeth <- load.rnb.diffmeth(paste0(analysis.dir, "/reports/differential_methylation_data/differential_rnbDiffMeth/"))
  # load the rnbset object
  rnb <- load.rnb.set(paste0(path, cell, "/", comp, "/reports/data_import_data/rnb.set_preprocessed/"))
  # load the annotation
  aa <- annotation(rnb, type = region, add.names = FALSE, include.regions = FALSE)

  # get the rank cuts
  rank.cuts.auto <- lapply(1:length(get.comparisons(diffMeth)), FUN = function(i) {
    cc <- names(get.comparisons(diffMeth))[i]
    ccc <- get.comparisons(diffMeth)[cc]
    dmt <- get.table(diffMeth, ccc, region, return.data.frame = TRUE)
    res <- RnBeads:::auto.select.rank.cut(dmt$diffMeth.p.adj.fdr, dmt$combinedRank, alpha = 0.1)
    return(as.integer(res))
  })

  # get the differential table
  for (i in 1:length(get.comparisons(diffMeth))) {
    cc <- names(get.comparisons(diffMeth))[i]
    ccc <- get.comparisons(diffMeth)[cc]
    dmt <- get.table(diffMeth, ccc, region, return.data.frame = TRUE)
    ccn <- ifelse(RnBeads:::is.valid.fname(cc), cc, paste("cmp", i, sep = ""))
    grp.names <- get.comparison.grouplabels(diffMeth)[ccc, ]
    auto.rank.cut <- rank.cuts.auto[[i]]
    ChrAccR:::cleanMem()
  }

  # add annotation to data.frame
  dmt <- cbind(aa, dmt)
  pval <- ifelse(region == "sites", "diffmeth.p.adj.fdr", "comb.p.adj.fdr")
  # add significant based on adjusted p-value
  dmt$isDiff <- dmt[, pval] < 0.05


  if (changeMethod == "l2fc") {
    if (region == "sites") {
      # select the necessary rows
      dmt <- data.table::as.data.table(dmt) %>%
        dplyr::mutate(log2FC = log2((mean.g2 + 1) / (mean.g1 + 1))) %>%
        dplyr::select(Chromosome, Start, End, Strand, log2FC, isDiff)
    } else {
      # select the necessary rows
      dmt <- data.table::as.data.table(dmt) %>%
        dplyr::mutate(log2FC = log2((mean.mean.g2 + 1) / (mean.mean.g1 + 1))) %>%
        dplyr::select(Chromosome, Start, End, Strand, log2FC, isDiff)
    }
  } else {
    if (region == "sites") {
      # select the necessary rows
      dmt <- data.table::as.data.table(dmt) %>%
        dplyr::mutate(log2FC = mean.g2 - mean.g1) %>%
        dplyr::select(Chromosome, Start, End, Strand, log2FC, isDiff)
    } else {
      # select the necessary rows
      dmt <- data.table::as.data.table(dmt) %>%
        dplyr::mutate(log2FC = mean.mean.g2 - mean.mean.g1) %>%
        dplyr::select(Chromosome, Start, End, Strand, log2FC, isDiff)
    }
  }
  return(na.omit(dmt))
}


# helper function for atac
cutL0.5fc2Padj05 <- function(dm) {
  abs(dm[, "log2FoldChange"]) > 0.5 & dm[, "padj"] < 0.05
}


### Helper function to compute overlap region betwen methylation and
### Written by: Irem B. Gündüz
prepareScatterMethAtac <- function(dt1, dt2, threshold = NULL) {
  if (!is.null(threshold)) {
    dt1 <- dt1 %>%
      dplyr::filter(abs(log2FC) < threshold)
    dt2 <- dt2 %>%
      dplyr::filter(abs(log2FC) < threshold)
  }
  # convert them to Granges
  dt1_gr <- dt1 %>%
    dplyr::rename(log2FC_1 = log2FC, isDiff_1 = isDiff) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T) %>%
    GenomicRanges::resize(500, "center") # group 1
  dt2_gr <- dt2 %>%
    dplyr::rename(log2FC_2 = log2FC, isDiff_2 = isDiff) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = T) %>%
    GenomicRanges::resize(500, "center") # group 1

  # merge the datasets
  comb <- mergeByOverlaps(dt1_gr, dt2_gr, minoverlap = 1, type = "within")
  ChrAccR:::cleanMem() # clean the memory
  return(as.data.frame(comb))
}


rnbeadsDensityScatter <- function(cell, comp, path, region) {
  if (!file.exists(paste0("/icbb/projects/igunduz/DARPA/Plots/scatter_", cell, "_", comp, ".pdf"))) {
    # get the analysis dir
    analysis.dir <- paste0(path, cell, "/", comp)

    # load the diffmeth object
    diffMeth <- load.rnb.diffmeth(paste0(analysis.dir, "/reports/differential_methylation_data/differential_rnbDiffMeth/"))

    # the rank cuts
    rank.cuts.auto <- 0

    for (i in 1:length(get.comparisons(diffMeth))) {
      cc <- names(get.comparisons(diffMeth))[i]
      ccc <- get.comparisons(diffMeth)[cc]
      dmt <- get.table(diffMeth, ccc, region, return.data.frame = TRUE)
      ccn <- ifelse(RnBeads:::is.valid.fname(cc), cc, paste("cmp", i, sep = ""))
      grp.names <- get.comparison.grouplabels(diffMeth)[ccc, ]
      auto.rank.cut <- rank.cuts.auto[[i]]
      df2p <- dmt # data frame to plot
      ChrAccR:::cleanMem()
    }
    # scatterplot based on adjusted p-value significance
    if (is.element("comb.p.adj.fdr", colnames(df2p))) {
      df2p$isDMP <- df2p[, "comb.p.adj.fdr"] < 0.05
      sparse.points <- 0.001
      dens.subsample <- TRUE
      pp <- create.densityScatter(df2p[, c("mean.mean.g2", "mean.mean.g1")],
        is.special = df2p$isDMP,
        dens.subsample = dens.subsample, sparse.points = sparse.points, add.text.cor = TRUE
      ) +
        labs(x = paste("Mean-Methylation (", grp.names[2], ")", sep = " "), y = paste("Mean-Methylation (", grp.names[1], ")", sep = "  ")) +
        coord_fixed() +
        theme_classic() +
        theme(legend.position = "none")
    }
    ChrAccR:::cleanMem()
  }
  return(pp)
}
