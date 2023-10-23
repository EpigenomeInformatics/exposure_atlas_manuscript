suppressPackageStartupMessages({
  library(dplyr)
  library(artemis)
  library(muLogR)
  library(SummarizedExperiment)
  library(Matrix)
  library(DelayedArray)
  library(HDF5Array)
  library(ggplot2)
  library(Seurat)
  library(irlba)
})
set.seed(43)
source("/icbb/projects/igunduz/exposure_atlas_manuscript/src/utils/lsi_utils.R")
k <- "sub11kpc30"

lsidir <- "/icbb/projects/igunduz/DARPA_analysis/artemis_031023/itLSI_res"
if (!dir.exists(lsidir)) {
  dir.create(lsidir)
}
umapdir <- "/icbb/projects/igunduz/DARPA_analysis/artemis_031023/umap"
if (!dir.exists(umapdir)) {
  dir.create(umapdir)
}

logger.start("Loading data")
# create a summarised experiment object
methSe <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/methSe_filtered.rds")
logger.completed()

logger.info(paste0("Number of regions:", nrow(methSe)))
# get the iterative LSI
res <- itLsi(methSe, doImpute = "pcimpute", it0mostVarRegs = 50000L, it1pcs = 1:30, it1mostVarRegs = 50000L, it2pcs = 1:30)

# save the result
logger.start("Saving the result")
saveRDS(res, paste0("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/itLSI_res/itLSI_res_", k, ".rds"))
logger.completed()

logger.start("Plotting UMAP...")
data <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/artemis_obj.rds")
filtered_sampleAnnot <- data@sampleAnnot[data@sampleAnnot$Cell_UID %in% row.names(res$umapCoord), ]

# plot the UMAP
df <- data.frame(
  res$umapCoord,
  cluster = res$clustAss,
  cellId = rownames(res$umapCoord),
  sampleSource = filtered_sampleAnnot,
  stringsAsFactors = FALSE
)

# Create a scatter plot
umap_plot <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = sampleSource.cell_type)) +
  geom_point(size = 0.25) +
  theme_classic() +
  labs(title = "UMAP Visualization by Cell Type")

# Display the plot
png(paste0("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/umap/itLSI_res_", k, "_umap_pca_imp_plot.png"))
umap_plot
dev.off()
logger.completed()
