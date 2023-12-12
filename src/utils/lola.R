lolaVolcanoPlotC19 <- function(cell, lolaDb, outputDir, pValCut = 2, region = "archrPeaks", database = "TF_motif_clusters",
                               signifCol = "qValue", cnames = c("loss", "gain"), n = 3, colorpanel = c()) {
  if (endsWith(outputDir, ".rds")) {
    bed.files <- readRDS(outputDir)
    df <- bed.files$region
    bed <- names(df)
    df <- df[[bed]]$archr_peaks %>%
      dplyr::filter(collection == database) # %>%
    # dplyr::filter(userSet %in% c("rankCut_500_hyper", "rankCut_500_hypo")) %>%
    # dplyr::mutate(condition = ifelse(userSet == "rankCut_500_hyper", "gain", "loss"))
    df <- muRtools:::lolaPrepareDataFrameForPlot(lolaDb, df,
      scoreCol = "log2OR", signifCol = signifCol,
      orderCol = "maxRnk", includedCollections = database,
      pvalCut = pValCut, maxTerms = Inf, perUserSet = FALSE,
      groupByCollection = TRUE, orderDecreasing = NULL
    )
    df$condition <- ifelse(grepl("hypo", df$userSet, ignore.case = TRUE), "loss", "gain")
    df$log2OR <- ifelse(df$condition == "loss", -df$log2OR, df$log2OR)
    df$name <- ifelse(database == "TF_motif_clusters", df$description, df$target)
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
    distinct(name, .keep_all = TRUE) %>%
    dplyr::filter(qValueLog > 2) %>%
    top_n(12, qValueLog)

  top_motifs_loss <- df %>%
    dplyr::filter(differential == "Differential", condition == "loss") %>%
    arrange(desc(qValueLog)) %>%
    distinct(name, .keep_all = TRUE) %>%
    dplyr::filter(qValueLog > 2) %>%
    top_n(12, qValueLog)

  top_motifs <- bind_rows(top_motifs_gain, top_motifs_loss)

  oddsRatioCol <- "log2OR"
  if (!is.element(oddsRatioCol, colnames(df))) {
    logger.error("Invalid LOLA result. Could not find a valid column containing odds ratios.")
  }
  pp <- plotVolcano(df, top_motifs, oddsRatioCol, signifCol)

  return(list(plot = pp, df = df))
}

plotVolcano <- function(df, top_motifs, oddsRatioCol, signifCol) {
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
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 14),
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_blank()
    ) +
    guides(alpha = "none", size = guide_legend(order = 1)) +
    coord_cartesian(ylim = c(0, max(df$qValueLog)), xlim = c(-max(df$log2OR), max(df$log2OR)))

  return(pp)
}
