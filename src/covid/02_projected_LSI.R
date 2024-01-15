set.seed(12) # set seed
outputDir <- "/scratch/icbb/covid_project_191223/"
mpal_dir <- "/icbb/projects/share/datasets/MPAL_2019/scATAC-Healthy-Hematopoiesis-191120.rds"
plot_dir <- "/icbb/projects/igunduz/projectLSI/Plot/"

# load libraries
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(ggplot2)
  library(rescueR)
})
harmony <- "cd34_final"

addArchRThreads(threads = 30) # set the cores
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)
mpal <- readRDS(mpal_dir) 

idxSample <- BiocGenerics::which(project$ClusterCellTypes %in% c("DC", "Mono_CD14", "Mono_CD16"))
cellsSample <- project$cellNames[idxSample]
project <- project[cellsSample, ]

# Subset the MPAL data to the cell types of interest
subset_mpal <- mpal[, mpal$Group %in% c("CD34_D7T1","CD34_D8T1","CD34_D9T1")]
desired_bio_classes <- c("01_HSC","05_CMP.LMPP","06_CLP.1","07_GMP","08_GMP.Neut","09_pDC","10_cDC",
 "11_CD14.Mono.1", "12_CD14.Mono.2","15_CLP.2")
subset_mpal <- subset_mpal[, subset_mpal$BioClassification %in% desired_bio_classes]

# Project the MPAL data onto the ArchR project
projected <- projectBulkATAC(
  ArchRProj = project,
  seATAC = subset_mpal,
  reducedDims = "IterativeLSI_covid",
  embedding = "UMAP_projection",
  n = 1,
  verbose = TRUE
)

# Plot the projected data
saveRDS(projected, file = paste0("/icbb/projects/igunduz/archr_project_011023/projected_", harmony, "_n", n, "_atacv2.rds"))

colData(mpal)$BioClassification <- sub(".*_", "", colData(mpal)$BioClassification)
metadata <- as.data.frame(colData(mpal))
metadata <- dplyr::select(metadata, BioClassification)
metadata <- dplyr::filter(metadata, BioClassification %in% c("cDC", "CD14.Mono.1", "CD14.Mono.2", "CMP.LMPP", "HSC", "GMP.Neut", "GMP")) # != "Unk")
metadata$Type <- rownames(metadata)

projected[[2]]$CellTypes <- project$ClusterCellTypes
projected[[2]] <- as.data.frame(projected[[2]])
projected[[2]] <- dplyr::filter(projected[[2]], CellTypes %in% c("Mono_CD14", "Mono_CD16", "DC")) # != "Unk")
projected[[2]] <- dplyr::select(projected[[2]], UMAP1, UMAP2, Type)
projected[[1]] <- merge(metadata, as.data.frame(projected[[1]]), by = "Type")
projected[[1]] <- dplyr::select(projected[[1]], BioClassification, UMAP1, UMAP2)
projected[[1]] <- dplyr::rename(projected[[1]], Type = BioClassification)

# Plot the projected data
plotProj <- rbind(projected[[2]], projected[[1]])
pal <- paletteDiscrete(unique(as.vector(plotProj[, 3])))
pal["scATAC"] <- "lightgrey"
p <- ggPoint(plotProj[, 1], plotProj[, 2], as.vector(plotProj[, 3]), rastr = TRUE, pal = pal)

# Save the plot
ggsave(paste0(plot_dir, "projection_", harmony, ".pdf"), plot = p, width = 10, height = 10)
ggsave(paste0(plot_dir, "projection_", harmony, ".jpeg"), plot = p, width = 10, height = 10)

# Compute the nearest neighbor in the ArchR data for each cell in the MPAL data
IterativeLSI_covid <- getReducedDims(project, "IterativeLSI_covid")
nn_index <- FNN::get.knnx(IterativeLSI_covid, projected[[3]], k = 5)$nn.index

# Assign the sample_exposure_group of the nearest neighbor to each cell in the MPAL data
projected[[1]]$sample_exposure_group <- project$sample_exposure_group[nn_index]

#Compute the proportions of each cell type in each sample_exposure_group
df <- data.frame(Type = projected[[1]]$Type, sample_exposure_group = projected[[1]]$sample_exposure_group)
contingency_table <- table(df)

# Plot the proportions of each cell type in each sample_exposure_group
stacked <- cellTypeProportionPlot(contingency_table,
  scale = TRUE, center = FALSE,
  groupName = "Exposure", 
  colorPalette = paletteDiscrete(unique(as.vector(plotProj[, 3]))),
  theme=theme_classic()
)
ggsave(paste0(plot_dir, "proportion_", harmony, ".pdf"), plot = stacked, width = 10, height = 10)