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


custom.densityScatter <- function(df2p, is.special = NULL,
                                  dens.subsample = FALSE, sparse.points = 0.01,
                                  dens.n = 100, add.text.cor = TRUE,
                                  color.by.direction = NULL, color.map = NULL) {
    if (!is.null(is.special)) {
        is.special[is.na(is.special)] <- FALSE
        df2p$is.special <- is.special
    }

    df2p <- na.omit(df2p)
    if (nrow(df2p) < 1) return(ggplot() + ggtitle("No data"))

    # Density estimation
    stable.bandwidth.fun <- function(x, eps = 1e-04) {
        r <- quantile(x, c(0.25, 0.75))
        h <- (r[2] - r[1]) / 1.34
        if (h == 0) h <- eps
        4 * 1.06 * min(sqrt(var(x)), h) * length(x)^(-1/5)
    }
    stable.h <- c(
        stable.bandwidth.fun(df2p[, 1]),
        stable.bandwidth.fun(df2p[, 2])
    )

    # Base plot with density
    pp <- ggplot(df2p, aes(x = !!sym(colnames(df2p)[1]), y = !!sym(colnames(df2p)[2]))) +
        stat_density2d(geom = "tile",
                       aes(alpha = after_stat(density)^0.25),
                       contour = FALSE, n = dens.n, h = stable.h, fill = "lightblue") +
        scale_alpha(range = c(0, 1))

    # Add sparse points (non-DMPs, low density)
    if (sparse.points > 0) {
        dens.ranks <- muRtools:::densRanks(df2p[, 1], df2p[, 2])
        thres <- if (sparse.points <= 1) ceiling(nrow(df2p) * sparse.points) else sparse.points
        df2p.loose <- df2p[dens.ranks <= thres, ]
        pp <- pp + geom_point(data = df2p.loose,
                              aes(x = !!sym(colnames(df2p)[1]), y = !!sym(colnames(df2p)[2])),
                              color = "1F78B4", size = 0.4)
    }

    # Add differential points with custom coloring
    if (!is.null(is.special) && any(is.special)) {
        df2p.special <- df2p[df2p$is.special, ]

        # If user supplied coloring variable
        if (!is.null(color.by.direction) && !is.null(color.map)) {
            df2p.special$direction <- color.by.direction[is.special]
            pp <- pp + geom_point(
                data = df2p.special,
                aes(x = !!sym(colnames(df2p)[1]), y = !!sym(colnames(df2p)[2]), color = direction),
                size = 1
            ) +
            scale_color_manual(values = color.map, name = "DMP direction")
        } else {
            # fallback: single color
            pp <- pp + geom_point(
                data = df2p.special,
                aes(x = !!sym(colnames(df2p)[1]), y = !!sym(colnames(df2p)[2])),
                color = "red", size = 1
            )
        }
    }

    # Add correlation text
    if (add.text.cor) {
        cc <- cor(df2p[, 1], df2p[, 2], use = "pairwise.complete.obs")
        txt.cor <- paste0("rho==", round(cc, 4))
        pp <- pp + annotate("text", x = max(df2p[, 1], na.rm = TRUE),
                            y = min(df2p[, 2], na.rm = TRUE),
                            label = txt.cor, parse = TRUE,
                            hjust = 1, vjust = 1, size = 4)
    }

    return(pp + coord_fixed() + theme_classic())
}

rnbeadsDensityScatter_sub <- function(diffMeth, region, plot_path, color_mapping) {
    plot_list <- list()

    for (i in seq_along(get.comparisons(diffMeth))) {
        cc <- names(get.comparisons(diffMeth))[i]
        ccc <- get.comparisons(diffMeth)[[cc]]
        dmt <- get.table(diffMeth, ccc, region, return.data.frame = TRUE)
        ccn <- ifelse(RnBeads:::is.valid.fname(cc), cc, paste0("cmp", i))
        grp.names <- get.comparison.grouplabels(diffMeth)[ccc, ]
        ChrAccR:::cleanMem()

        if ("comb.p.adj.fdr" %in% colnames(dmt)) {
            dmt$isDMP <- dmt$comb.p.adj.fdr <= 0.05

            # Define DMP direction
            dmt$direction <- NA
            dmt$direction[dmt$isDMP & dmt$mean.mean.g1 > dmt$mean.mean.g2] <- paste0(grp.names[1], " higher")
            dmt$direction[dmt$isDMP & dmt$mean.mean.g1 < dmt$mean.mean.g2] <- paste0(grp.names[2], " higher")

            # Build color map from global color_mapping
            group_labels <- c(
                paste0(grp.names[1], " higher"),
                paste0(grp.names[2], " higher")
            )

            group_colors <- setNames(
                c(color_mapping[grp.names[1]], color_mapping[grp.names[2]]),
                group_labels
            )

            # Create plot
            pp <- custom.densityScatter(
                dmt[, c("mean.mean.g2", "mean.mean.g1")],
                is.special = dmt$isDMP,
                color.by.direction = dmt$direction,
                color.map = group_colors
            ) 

            plot_list[[ccn]] <- pp
            pdf_file <- file.path(plot_path, paste0(ccn, "_density_scatter.pdf"))
            ggsave(pdf_file, plot = pp, width = 8, height = 6)
        }
    }

    return(plot_list)
}
