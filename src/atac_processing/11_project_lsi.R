set.seed(12) # set seed
outputDir <- "/icbb/projects/igunduz/archr_project_011023/"
mpal_dir <- "/icbb/projects/igunduz/scATAC-Healthy-Hematopoiesis-191120.rds"
source("/icbb/projects/igunduz/exposure_atlas_manuscript/src/utils/archr_utils.R")
# load libraries
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(ggplot2)
})
n <- 1
harmony <- "withHARMONY"

addArchRThreads(threads = 30) # set the cores
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)
mpal <- readRDS(mpal_dir) # Summarised experiment containing counts and metadata for the MPAL data
# colData(mpal)$BioClassification <- sub(".*_", "", colData(mpal)$BioClassification)

# Filter unknown cells
# mpal <- mpal[!colData(mpal)$BioClassification== "Unk", ]

logger::log_info(paste0("Setting n to ", n))
# Project the MPAL data onto the ArchR project
projected <- projectBulkATAC_vs2(
  ArchRProj = project,
  seATAC = mpal,
  reducedDims = "IterativeLSI",
  embedding = "UMAPHarmony",
  n = n,
  verbose = TRUE
)

# Plot the projected data
saveRDS(projected, file = paste0("/icbb/projects/igunduz/archr_project_011023/projected_", harmony, "_n", n, "_atacv2.rds"))
projected <- readRDS(file = paste0("/icbb/projects/igunduz/archr_project_011023/projected_", harmony, "_n", n, "_atacv2.rds"))

colData(mpal)$BioClassification <- sub(".*_", "", colData(mpal)$BioClassification)
metadata <- as.data.frame(colData(mpal))
metadata <- dplyr::select(metadata, BioClassification)
metadata <- dplyr::filter(metadata, BioClassification != "Unk")
metadata$Type <- rownames(metadata)

projected[[1]] <- merge(metadata, as.data.frame(projected[[1]]), by = "Type")
projected[[1]] <- dplyr::select(projected[[1]], BioClassification, UMAP1, UMAP2)
projected[[1]] <- dplyr::rename(projected[[1]], Type = BioClassification)
# Plot the projected data
# projected[[1]][,3] <- sub(":.*", "", as.vector(projected[[1]]$Type))
# stringr::str_split(as.vector(head(projected[[1]]$Type)),pattern="\\_",simplify=TRUE)[,2]
plotProj <- rbind(projected[[2]], projected[[1]])
pal <- paletteDiscrete(unique(as.vector(plotProj[, 3])))
pal["scATAC"] <- "lightgrey"
p <- ggPoint(plotProj[, 1], plotProj[, 2], as.vector(plotProj[, 3]), rastr = TRUE, pal = pal)

ggsave(paste0("/icbb/projects/igunduz/Plot-Bulk-Heme-Overlay-Projection_n", n, "_", harmony, "v2.pdf"), plot = p, width = 10, height = 10)
ggsave(paste0("/icbb/projects/igunduz/Plot-Bulk-Heme-Overlay-Projection_n", n, "_", harmony, "v2.jpeg"), plot = p, width = 10, height = 10)
