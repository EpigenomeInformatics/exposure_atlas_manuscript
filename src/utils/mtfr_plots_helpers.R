computeL2FCdevs <- function(deviations, computezscore = FALSE, differential_deviation_test = TRUE,
                            alternative = c("two.sided"), group = NULL, parametric = TRUE, chromvar_obj = NULL, grp1name = NULL, grp2name = NULL) {
  if (is.null(row.names(deviations))) {
    stop("Please provide motifs as rownames in deviations data.frame!")
  }
  if (is.null(grp1name) & is.null(group)) {
    stop("Please provide at least one group names!")
  }
  if (!is.null(grp1name)) {
    # Identify the columns corresponding to each group
    sidx.grp1 <- grep(grp1name, colnames(deviations))
    sidx.grp2 <- setdiff(1:ncol(deviations), sidx.grp1)
    grp2name <- "group2"
  }
  if (!is.logical(computezscore)) {
    stop("Please provide a logical value for computezscore!")
  }
  if (!is.null(group)) {
    unique_groups <- unique(group)
    grp1name <- unique_groups[grep("trl$", unique_groups)]
    grp2name <- unique_groups[which(unique_groups != grp1name)]
    sidx.grp1 <- which(group == grp1name)
    sidx.grp2 <- which(group == grp2name)
  }
  # run differential deviation test
  diff <- methylTFR:::differential_deviation_test(deviations, groups = group, parametric = parametric, alternative = alternative)
  colnames(diff) <- c("motifs", "p_value", "p_value_adjusted")

  # Subset the matrix based on the group columns
  # deviations <- as.data.frame(deviations)
  m1 <- deviations[, sidx.grp1, drop = FALSE]
  m2 <- deviations[, sidx.grp2, drop = FALSE]

  # Calculate the row means for each group
  cn_mean1 <- paste0("meanDeviationgrp", "1_", grp1name)
  cn_mean2 <- paste0("meanDeviationgrp", "2_", grp2name)
  diff[, cn_mean1] <- rowMeans(m1, na.rm = FALSE)
  diff[, cn_mean2] <- rowMeans(m2, na.rm = FALSE)

  # Calculate the log2FC column
  diff[, "log2FC"] <- log2((diff[, cn_mean1] + 1) / (diff[, cn_mean2] + 1))

  if (computezscore) {
    # Calculate the z-scores
    # motifs <- rownames(deviations)
    if (!is.null(chromvar_obj)) {
      deviations <- chromVAR::deviationScores(chromvar_obj)
    } else {
      deviations <- computeZScore(deviations)
    }
    deviations <- as.data.frame(deviations)
    rownames(deviations) <- NULL
    diff[, cn_mean1] <- rowMeans(deviations[, sidx.grp1, drop = FALSE], na.rm = TRUE)
    diff[, cn_mean2] <- rowMeans(deviations[, sidx.grp2, drop = FALSE], na.rm = TRUE)
    diff[, "zDiff"] <- diff[, cn_mean1] - diff[, cn_mean2]
  }
  diff[, "negLog10padj"] <- -log10(diff[, "p_value_adjusted"])
  return(diff)
}

deviationHeatmap <- function(
    mat, ann_df, ann_col = "Condition",
    groups = c("Control", "Severe"), package = "methylTFR",
    fill_col = c(Control = "#4F609C", Severe = "#C03830"), colors = "cptcity.arendal_temperature",
    clustering_method = "within", cluster_rows = TRUE, row_order = NULL) {
  # Cluster rows and columns based on the selected method
  if (clustering_method == "within") {
    dend_cols <- cluster_within_group(mat, ann_df[[ann_col]])
  } else if (clustering_method == "between") {
    dend_cols <- cluster_between_groups(mat, ann_df[[ann_col]])
  }

  # Create a HeatmapAnnotation object for column annotations
  # column_ha <- HeatmapAnnotation(cluster = anno_block(gp = gpar(fill = fill_col),labels = groups ))
  # TODO delete the one above 4F609C

  column_ha <- HeatmapAnnotation("Condition" = ann_df[[ann_col]], col = list("Condition" = fill_col))

  c <- if (package == "methylTFR") {
    muRtools::colpal.cont(length(ann_df[[ann_col]]), colors)
  } else {
    rev(muRtools::colpal.cont(length(ann_df[[ann_col]]), colors))
  }
  if (!is.null(row_order)) {
    mat <- mat[row_order, ]
  }
  # Create a heatmap object
  heatmap_obj <- Heatmap(mat,
    cluster_rows = cluster_rows,
    cluster_columns = dend_cols,
    # column_split = 3,
    top_annotation = column_ha,
    # name = paste0("TF Deviation Z-scores (", package, ")"),
    col = c,
    heatmap_legend_param = list(
      title = paste0("TF Deviation Z-scores (", package, ")") # ,
      # at = c(-max(round(mat)), 0, max(round(mat)))
    ),
    show_column_names = FALSE # ,
    # heatmap_legend_param = list(title = "Z-score", at = seq(-3, 3, by = 1)),
    # width = unit(nrow(mat) / 2, "cm")
    # height = unit(40, "cm")
  )
  row_order <- row_order(heatmap_obj)
  row_order <- rownames(mat[row_order, ])
  return(list(hm = heatmap_obj, row_order = row_order))
}

deviationW2annotHeatmap <- function(
    mat, ann_df,
    cluster_col = "CellType",
    fill_col_cell_type = NULL, # Colors for Cell Type
    fill_col_condition = NULL, # Colors for Condition
    colors = "cptcity.arendal_temperature",
    clustering_method = "within",
    cluster_rows = TRUE,
    package = "methylTFR",
    row_order = NULL,
    show_column_names = FALSE) {
  # Cluster rows and columns based on the selected method
  if (!is.null(cluster_col) & is.character(cluster_col) & length(cluster_col) == 1) {
    if (clustering_method == "within") {
      dend_cols <- cluster_within_group(mat, ann_df[[cluster_col]])
    } else if (clustering_method == "between") {
      dend_cols <- cluster_between_groups(mat, ann_df[[cluster_col]])
    }
  }

  if (!is.null(fill_col_cell_type) & !is.null(fill_col_condition)) {
    logger.info("Using more than one annotation for coloring the heatmap!")
    column_ha <- HeatmapAnnotation(
      df = ann_df, col = list(Cell = fill_col_cell_type, Condition = fill_col_condition)
    )
  } else {
    logger.info("Using only one annotation for coloring the heatmap!")
    column_ha <- HeatmapAnnotation(
      df = ann_df, col = list(Cell = fill_col_cell_type)
    )
  }
  c <- if (package == "methylTFR") {
    muRtools::colpal.cont(nrow(ann_df), colors)
  } else {
    colors.cv <- ChrAccR::getConfigElement("colorSchemesCont")
    colors.cv <- colors.cv[[".default.div"]]
    c <- grDevices::colorRampPalette(colors.cv)(nrow(ann_df))

    # rev(muRtools::colpal.cont(nrow(ann_df), colors))
  }

  if (!is.null(row_order)) {
    mat <- mat[row_order, ]
  }

  # Create a heatmap object
  heatmap_obj <- Heatmap(
    mat,
    cluster_rows = cluster_rows,
    cluster_columns = dend_cols,
    top_annotation = column_ha,
    col = c,
    heatmap_legend_param = list(
      title = paste0("TF Deviation Z-scores (", package, ")") # ,
      # at = c(-2,-1,0,1, 2)
    ),
    show_column_names = FALSE
  )

  row_order <- row_order(heatmap_obj)
  row_order <- rownames(mat[row_order, ])

  return(list(hm = heatmap_obj, row_order = row_order))
}


diagDivCellHeatmap <- function(ml, mr, col.l = NULL, col.r = NULL, name.l = "Lower left", name.r = "Upper right", ...) {
  require(ComplexHeatmap)
  if (nrow(ml) != nrow(mr)) stop("Numbers of rows of the two matrices must match")
  if (ncol(ml) != ncol(mr)) stop("Numbers of columns of the two matrices must match")

  if (!is.function(col.l)) col.l <- getColorFun(ml, col.l)
  if (!is.function(col.r)) col.r <- getColorFun(mr, col.r)
  colPalLegend_l <- col.l
  if (is.character(ml) || is.factor(ml)) {
    colPalLegend_l <- environment(col.l)[["mapVec"]]
  }
  colPalLegend_r <- col.r
  if (is.character(mr) || is.factor(mr)) {
    colPalLegend_r <- environment(col.r)[["mapVec"]]
  }

  dummyColor <- circlize::colorRamp2(seq(0, 1, length.out = 2), rep("grey", 2))
  res <- Heatmap(ml,
    col = colPalLegend_l, rect_gp = gpar(type = "none"),
    cell_fun = function(j, i, x, y, w, h, fill) {
      # grid.rect(x=x, y=y, width=width, height=height, gp=gpar(col="grey", fill="grey"))
      grid.polygon(
        unit.c(x - 0.5 * w, x - 0.5 * w, x + 0.5 * w),
        unit.c(y - 0.5 * h, y + 0.5 * h, y - 0.5 * h),
        gp = gpar(col = NA, fill = col.l(ml[i, j]))
      )
      grid.polygon(
        unit.c(x + 0.5 * w, x + 0.5 * w, x - 0.5 * w),
        unit.c(y + 0.5 * h, y - 0.5 * h, y + 0.5 * h),
        gp = gpar(col = NA, fill = col.r(mr[i, j]))
      )
    },
    name = name.l,
    ...
  )
  # dummy heatmap for color legend
  dummyM <- matrix(rep(NA, length.out = nrow(mr)), nrow = nrow(ml), ncol = 1)
  rownames(dummyM) <- rownames(ml)
  dummyHm <- Heatmap(dummyM, col = colPalLegend_r, width = unit(0, "mm"), name = name.r)

  return(res + dummyHm)
}

computeZScore <- function(mat) {
  mat <- (mat - matrixStats::rowMeans2(mat)) / matrixStats::rowSds(mat)
  mat[base::is.nan(mat)] <- 0
  return(mat)
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
