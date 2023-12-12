set.seed(42)
logger::log_info("Loading libraries...")

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(methylTFR)
  library(muLogR)
  library(gplots)
  library(muRtools)
  library(ComplexHeatmap)
  library(factoextra)
  library(ggfortify)
})
  get_groupname <- function(x) {
    return(unlist(strsplit(x, split = "_"))[1])
  }
plot_dir <- "/icbb/projects/igunduz/exposure_atlas_manuscript/Figures"
if (!dir.exists(plot_dir)) dir.create(plot_dir)
motifsetList <- c("altius", "jaspar2020")
logger.info("Starting analysis for Monocyte in chromVAR and methylTFR")
conditions <- c("HIV_acu", "HIV_chr", "HIV_ctrl", "OP_high", "OP_low", "OP_med")
cells <- c("B-cell", "Monocyte", "NK-cell", "Tc-Mem", "Tc-Naive", "Th-Mem", "Th-Naive")

cell_type_colors <- c(
  "Other-cell" = "#CCCCCC",
  "B-cell" = "#AE017E",
  "Monocyte" = "#CC4C02",
  "NK-cell" = "#A65628",
  "Th-Mem" = "#41B6C4",
  "Tc-Mem" = "#4292C6",
  "Tc-Naive" = "#888FB5",
  "Th-Naive" = "#C7E9B4"
)

for (motifset in motifsetList) {
  mtfr_dir <- paste0("/icbb/projects/igunduz/DARPA_analysis/mtfr_results_241123/shared/updated/", motifset)
  result <- list.files(mtfr_dir, pattern = "deviations.RDS", full.names = TRUE)

  # Loop through combinations and read the data
  mtfr_devs <- list()
  for (x in seq_along(result)) {
    # path <- paste0(mtfr_dir, "/", result[x], "deviations.RDS")
    path <- result[x]
    cur_dev <- readRDS(path)
    cur_dev <- deviations(cur_dev) # deviationZScores(cur_dev)
    mtfr_devs[[x]] <- as.matrix(cur_dev)
  }
  # Combine the lists into data frames
  mtfr_devs <- do.call(base::cbind, mtfr_devs)

  # Get the group names
  groups <- unlist(lapply(FUN = get_groupname, X = colnames(mtfr_devs)))
  tdf <- as.data.frame(t(mtfr_devs))
  tdf$cell_type <- groups

  # skip cell_type column
  k <- ifelse(motifset == "altius", 287, 633)
  pca <- prcomp(tdf[, -k]) # 287 for jaspar2020
  fn_pca_ind <- file.path(plot_dir, paste0("pca_", motifset, "_rawscores.pdf"))
  pdf(fn_pca_ind)
  autoplot(pca,
    data = tdf,
    colour = "cell_type",
    main = "PCA",
    size = 5
  ) +
    theme_classic() +
    scale_color_manual(values = cell_type_colors) +
    theme(legend.position = "bottom")
  dev.off()
}
