suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
})
set.seed(12) # set seed

addArchRThreads(threads = 30) # set the cores
addArchRGenome("hg38") # set the reference genome

project <- ArchR::loadArchRProject("/icbb/projects/igunduz/archr_project_011023/", showLogo = FALSE)

# subset the T cells subtypes
idxSample <- BiocGenerics::which(project$ClusterCellTypes %in% c("T_mix", "T_mem_CD8", "T_mem_CD4", "T_naive", "T_mait"))
cellsSample <- project$cellNames[idxSample]
project2 <- project[cellsSample, ]

# save the project to the subdirectory
outputDir <- "/icbb/projects/igunduz/Tcell_subset/"
saveArchRProject(project2, outputDirectory = outputDir, load = TRUE)
project <- project2
project <- project[cellsSample, ]
arrows <- list.files(paste0(outputDir, "ArrowFiles"), full.names = TRUE)
project@sampleColData <- DataFrame(ArrowFiles = arrows)
rownames(project@sampleColData) <- gsub(x = gsub(x = arrows, ".arrow", ""), paste0(outputDir, "ArrowFiles/"), "")

# Get the row names
row_names <- rownames(project@sampleColData)
# Create a SimpleList with the row names
simple_list <- SimpleList(vector("list", length = 92))
names(simple_list) <- row_names

project@sampleMetadata <- simple_list
project@projectMetadata$outputDirectory <- outputDir
# project@projectMetadata$GroupCoverages$ClusterCellTypes$coverageMetadata$File = gsub("DARPA/ArchRProject_5x/ATAC_processed","ATAC_processed_final",project@projectMetadata$GroupCoverages$ClusterCellTypes$coverageMetadata$File)

# save the project to the subdirectory
saveRDS(project, paste0(outputDir, "Save-ArchR-Project.rds"))
