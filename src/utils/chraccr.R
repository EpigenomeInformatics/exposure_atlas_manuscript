plotMAwithChrAccR <- function(cell, outputDir, i, region) {
  bed.files <- list.files(paste0(outputDir, cell, "/reports/differential_data"), pattern = "diffTab", full.names = T)
  bed.files <- grep(bed.files, pattern = region, fixed = T, value = T)[i]
  dm <- read.delim(bed.files)
  df2p.ma <- dm[, c("log2BaseMean", "log2FoldChange")]
  isDiff <- cutL0.5fc2Padj05(dm)
  isDiff[is.na(isDiff)] <- FALSE
  t <- as.character(as.data.frame(table(isDiff))$Freq[2])
  t <- ifelse(is.na(t), 0, t)
  pp <- muRtools::create.densityScatter(df2p.ma, is.special = isDiff, sparse.points = 0.001) +
    theme_classic() +
    guides(alpha = "none") +
    annotate("text", x = max(df2p.ma) - 1, y = min(df2p.ma), label = paste0("# of diff. peaks = ", t), size = 10) + # TODO change this
    theme(
      axis.title.x = element_text(size = 15),
      axis.title.y = element_text(size = 15)
    )
  return(pp)
}

### Scatter plot of log2FoldChange per two group condition
### Written by: Irem B. Gündüz
plotScatterL2FC <- function(datatable, y_lab, x_lab, comb, group1, group2, textsize = 10, point_size = 3, label = FALSE,
                            meth_label = "differential in METH", atac_label = "differential in ATAC", max_overlaps = 500) {
  p1 <- ggplot(datatable, aes(
    x = .data[[group1]], y = .data[[group2]],
    color = interaction(isDiff_1, isDiff_2, sep = "-", lex.order = TRUE)
  )) +
    geom_point(size = point_size) +
    theme_classic() +
    ylab(y_lab) +
    xlab(x_lab) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    theme(legend.position = "bottom", axis.title = element_text(size = 16), legend.text = element_text(size = 14)) +
    scale_color_manual(paste0("Is differential? ", comb),
      values = c("#a6cee3", "#006400", "#1a1ce3", "#e31a1c"),
      labels = c("Not Differential", meth_label, atac_label, "differential in both")
    ) +
    scale_alpha_continuous(range = c(0.1, 1))
  if (label) {
    p1 <- p1 +
      geom_text_repel(
        data = datatable[datatable$isDiff_1 | datatable$isDiff_2, ],
        aes(x = .data[[group1]], y = .data[[group2]], label = name),
        color = "black", size = textsize, box.padding = 0.5,
        segment.color = "black", segment.size = 0.1,
        max.overlaps = max_overlaps # Increase max.overlaps as needed
      )
  }
  ChrAccR:::cleanMem()
  return(p1)
}


### Function to prepare differentially accessible regions for scatter plots
### Written by : Irem B. Gündüz  & Fabian Müller
prepareDARforPlot <- function(cell, outputDir, i, region = "archrPeaks") {
  # get the path of bed.files
  bed.files <- list.files(paste0(outputDir, cell, "/reports/differential_data"), pattern = "diffTab", full.names = T)
  bed.files <- grep(bed.files, pattern = region, fixed = T, value = T)[i]

  dm <- read.delim(bed.files)
  isDiff <- cutL0.5fc2Padj05(dm[, c("log2FoldChange", "padj")]) # ,padj=0.05)
  isDiff[is.na(isDiff)] <- FALSE # fill missing as FALSE
  dm$isDiff <- isDiff # Add differentials
  dm <- data.table::as.data.table(dm) %>%
    dplyr::rename(log2FC = log2FoldChange, Chromosome = chrom, Start = chromStart, End = chromEnd, Strand = strand) %>%
    dplyr::select(Chromosome, Start, End, Strand, log2FC, isDiff)
  ChrAccR:::cleanMem()
  return(na.omit(dm))
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
