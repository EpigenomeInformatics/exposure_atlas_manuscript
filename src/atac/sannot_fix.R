## Load Libraries
suppressPackageStartupMessages({library(ArchR)
library(dplyr)
library(readr)
library(openxlsx)
})
set.seed(12) # set seed

# Load ArchR project
outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/"
project <- ArchR::loadArchRProject(outputDir, force = T)

# Extract cellColData
cell_metadata <- as.data.frame(getCellColData(project))

# Add CellNames as column (optional)
cell_metadata$CellNames <- rownames(cell_metadata)

# Summarize by sample
sample_summary <- cell_metadata %>%
  group_by(Sample) %>%
  summarise(
    nCells = n(),                                # Number of cells per sample
    mean_FRIP = mean(FRIP, na.rm = TRUE),        # Average FRIP
    mean_TSS = mean(TSSEnrichment, na.rm = TRUE),# Average TSSEnrichment
    mean_nFrags = mean(nFrags, na.rm = TRUE),  # Average number of fragments
    .groups = "drop"
  )

# For example, exposure type or cell types can be inferred if you have that info in cellColData
sample_annotation <- cell_metadata %>%
  group_by(Sample) %>%
  summarise(
    exposure_type  = first(na.omit(sample_exposure_type)),
    exposure_group = first(na.omit(sample_exposure_group)),
    .groups = "drop"
  )

# Join to summary
sample_summary <- sample_summary %>%
  left_join(sample_annotation, by = "Sample")

# Write the summarized sample-level metadata to Excel
write.xlsx(sample_summary, file = "SampleLevel_Metadata_Supplementary.xlsx", rowNames = FALSE)
