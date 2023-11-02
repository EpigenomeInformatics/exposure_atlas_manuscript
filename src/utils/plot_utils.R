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
      pp <- create.densityScatter(df2p[, c("mean.mean.g2","mean.mean.g1")],
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



lolaVolcanoPlotC19 <- function(cell, lolaDb, outputDir, pValCut = 2, region = "archrPeaks", database = "TF_motif_clusters",
                               signifCol = "qValue", cnames = c("loss", "gain"), n = 3, colorpanel = c()) {
  if (endsWith(outputDir, ".rds")) {
    bed.files <- readRDS(outputDir)
    df <- bed.files$region
    bed <- names(df)
    df <- df[[bed]]$archr_peaks %>%
      dplyr::filter(collection == database) %>%
      dplyr::filter(userSet %in% c("rankCut_500_hyper", "rankCut_500_hypo")) %>%
      dplyr::mutate(condition = ifelse(userSet == "rankCut_500_hyper", "gain", "loss"))
    df <- muRtools:::lolaPrepareDataFrameForPlot(lolaDb, df,
      scoreCol = "log2OR", signifCol = signifCol,
      orderCol = "maxRnk", includedCollections = database,
      pvalCut = pValCut, maxTerms = Inf, perUserSet = FALSE,
      groupByCollection = TRUE, orderDecreasing = NULL
    )
    df$log2OR <- ifelse(df$condition == "loss", -df$log2OR, df$log2OR)
    df$term <- sub("\\s*\\[.*", "", df$term)
    df$name <- df$term
  } else {
    bed.files <- list.files(paste0(outputDir, cell, "/reports/differential_data"), pattern = "lolaRes", full.names = TRUE)
    bed.files <- grep(bed.files, pattern = region, value = TRUE)
    bed.files_l <- grep(pattern = paste0("_cutL2fcPadj05loss.rds"), bed.files, value = TRUE)[n]
    bed.files_g <- grep(pattern = paste0("_cutL2fcPadj05gain.rds"), bed.files, value = TRUE)[n]
    bed.files <- c(bed.files_l, bed.files_g)
    names(bed.files) <- cnames
    data <- lapply(names(bed.files), function(bed) {
      a <- readRDS(bed.files[[bed]]) %>%
        dplyr::filter(collection == database)
      a <- muRtools:::lolaPrepareDataFrameForPlot(lolaDb, a,
        scoreCol = "log2OR", signifCol = signifCol,
        orderCol = "maxRnk", includedCollections = database,
        pvalCut = pValCut, maxTerms = Inf, perUserSet = FALSE,
        groupByCollection = TRUE, orderDecreasing = NULL
      )
      if (!is.null(a)) {
        a <- a %>%
          dplyr::mutate(condition = rep(bed, nrow(a)), log2OR = log2(oddsRatio))
        if (bed == "loss") {
          a <- a %>%
            dplyr::mutate(log2OR = -log2OR)
        }
      }
      a
    })
    df <- rlist::list.rbind(data)
    df$name <- gsub("\\s*\\[.*\\]", "", df$name)
  }
  if (length(bed.files) == 0) {
    stop("No files found with the specified pattern.")
  }
  if (signifCol == "qValue") {
    signifCol <- "qValueLog"
    df$qValueLog <- -log10(df$qValue)
    df$differential <- ifelse(df$qValueLog > 1.3, "Differential", "Non-differential")
  }


  top_motifs_gain <- df %>%
    dplyr::filter(differential == "Differential", condition == "gain") %>%
    arrange(desc(qValueLog)) %>%
    dplyr::filter(qValueLog > 2) %>%
    top_n(10, qValueLog)

  top_motifs_loss <- df %>%
    dplyr::filter(differential == "Differential", condition == "loss") %>%
    arrange(desc(qValueLog)) %>%
    dplyr::filter(qValueLog > 2) %>%
    top_n(10, qValueLog)

  top_motifs <- bind_rows(top_motifs_gain, top_motifs_loss)

  oddsRatioCol <- "log2OR"
  if (!is.element(oddsRatioCol, colnames(df))) {
    logger.error("Invalid LOLA result. Could not find a valid column containing odds ratios.")
  }
  pp <- plotVolcano(df, top_motifs, oddsRatioCol, signifCol)

  return(list(plot = pp, df = df))
}

plotVolcano <- function(df, top_motifs, oddsRatioCol, signifCol) {
  library(ggplot2)
  library(ggrepel)

  pp <- ggplot(df) +
    aes_string(oddsRatioCol, signifCol) +
    geom_point()
  pp <- pp +
    theme_classic() +
    geom_point(data = df %>% filter(differential == "Non-differential"), color = "black") +
    ggrepel::geom_text_repel(aes(label = name),
      data = top_motifs,
      size = 20,
      color = ifelse(top_motifs$condition == "gain", "red", "blue"),
      min.segment.length = 0,
      seed = 42,
      box.padding = 1,
      max.overlaps = Inf,
      segment.size = 0.5
    ) +
    ylab("-log10(Q-Value)") +
    xlab("log2(Odds-Ratio)") +
    theme(
      axis.text = element_text(size = 14),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_blank()
    ) +
    guides(alpha = "bottom", size = guide_legend(order = 1)) +
    coord_cartesian(ylim = c(0, max(df$qValueLog)), xlim = c(-max(df$log2OR), max(df$log2OR)))

  return(pp)
}
