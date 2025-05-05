#!/usr/bin/env Rscript

#####################################################################
# 09_gene_exp_vs_activity.R
# Created on 05-05-2025 by Irem Gunduz
# Correlation analysis between gene activity and gene expression
######################################################################

## Load Libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(Seurat)
  library(rescueR)
  library(ArchR)
  library(ggrepel)
  library(muLogR)
  library(SummarizedExperiment)
  library(Seurat)
  library(biomaRt)
})
set.seed(12) # set seed

# Load the data and set the output directory
plot_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/figures/"
sannot_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/sample_annots/"
archr_dir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/"

# Load the ArchR project
project <- ArchR::loadArchRProject(archr_dir, showLogo = FALSE)

# Subset the monocytes from archr project
idxSample <- BiocGenerics::which(project$ClusterCellTypes %in% c("Mono_CD14"))
cellsSample <- project$cellNames[idxSample]
project_mono <- project[cellsSample, ]

# Subset C19_sev and C19_ctrl
project_mono <- project_mono[project_mono$sample_exposure_group %in% c("C19_sev", "C19_ctrl"), ]
cellsSample <- project_mono$cellNames[idxSample]
project_mono <- project_mono[cellsSample, ]

# Get the protein coding genes
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
protein_coding_genes <- getBM(
  attributes = c("external_gene_name"),
  filters = "biotype",
  values = "protein_coding",
  mart = ensembl
)
protein_coding_genes <- protein_coding_genes$external_gene_name
length(protein_coding_genes)
# [1] 19457

# Get marker scores
markersGS <- getMarkerFeatures(
  ArchRProj = project_mono,
  useMatrix = "GeneScoreMatrix",
  useGroups = "C19_sev",
  bgdGroups = "C19_ctrl",
  groupBy = "sample_exposure_group",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# Get the log2FC info
markerList <- getMarkers(markersGS, cutOff = "FDR <= 100") # Not filter here
atac <- as.data.frame(markerList$C19_sev)

# Filter the data to keep only the protein coding genes
atac_protein_coding <- atac[atac$name %in% protein_coding_genes, ]
dim(atac_protein_coding)
# [1] 18248     9

# Save the results
saveRDS(atac_protein_coding, file.path(sannot_dir, "gene_activity_protein_coding_diffs.rds"))

######################################################################
# scRNA-seq data gene expression for the protein coding genes
######################################################################

seu <- readRDS("/icbb/projects/igunduz/blish_awilk_covid_seurat.rds")
meta <- read.csv(url("https://github.com/ajwilk/COVID_scMultiome/raw/main/data/scRNA/patient_metadata.csv"))

# Merge new metadata
current_meta <- seu@meta.data # [, colnames(seu@meta.data) != "current_severity"]
current_meta$current_severity <- as.numeric(current_meta$current_severity)

# Define the severity breaks and labels
severity_breaks <- c(-Inf, 0, 3, 5, Inf)
severity_labels <- c("C19_ctrl", "C19_mild", "C19_mod", "C19_sev")
current_meta$condition <- cut(current_meta$current_severity, breaks = severity_breaks, labels = severity_labels, include.lowest = TRUE)
seu@meta.data <- current_meta

# Subset the seu object to the common ids
sannot <- readr::read_csv("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/sample_annots/ATAC_metadata_covid.csv", show_col_types = FALSE)
seu_meta <- seu@meta.data
common_ids <- seu_meta$record_id[seu_meta$record_id %in% sannot$record_id]
seu_subset <- subset(x = seu, subset = record_id %in% common_ids)
seu <- seu_subset

# Differential expression analysis
seu$celltype.condition <- paste(seu$cell.type, seu$condition, sep = "_")
Idents(seu) <- "celltype.condition"

# Find the DE genes between ctrl and sev
bulk.mono.de <- FindMarkers(
  object = seu,
  ident.1 = "CD14 Mono_C19_sev",
  ident.2 = "CD14 Mono_C19_ctrl",
  test.use = "DESeq2"
)
bulk.mono.de$name <- rownames(bulk.mono.de)

# Save the results
saveRDS(bulk.mono.de, paste0(sannot_dir, "gene_expression_protein_coding_diffs.rds"))

######################################################################
# Correlation analysis
######################################################################

atac_protein_coding <- readRDS(file.path(sannot_dir, "gene_activity_protein_coding_diffs.rds"))
bulk.mono.de <- readRDS(paste0(sannot_dir, "gene_expression_protein_coding_diffs.rds"))

# Find common genes and subset
common_genes <- intersect(bulk.mono.de$name, atac_protein_coding$name)
bulk.mono_de_common <- bulk.mono.de[bulk.mono.de$name %in% common_genes, ]
atac_common <- atac_protein_coding[atac_protein_coding$name %in% common_genes, ]

# Organize the column names
atac_common <- atac_common[, c("name", "Log2FC", "FDR")]
colnames(atac_common) <- c("name", "Log2FC_atac", "FDR_atac")
bulk.mono_de_common <- bulk.mono_de_common[, c("name", "avg_log2FC", "p_val_adj")]
colnames(bulk.mono_de_common) <- c("name", "Log2FC_rna", "FDR_rna")

# Merge
atac_common$isDiff_1 <- ifelse(atac_common$FDR_atac <= 0.05 & abs(atac_common$Log2FC_atac) > 0.5, T, F)
bulk.mono_de_common$isDiff_2 <- ifelse(bulk.mono_de_common$FDR_rna <= 0.05 & abs(bulk.mono_de_common$Log2FC_rna) > 0.5, T, F)
bulk.mono_de_common <- bulk.mono_de_common[abs(bulk.mono_de_common$Log2FC_rna) <= 5, ]
atac_common <- atac_common[abs(atac_common$Log2FC_atac) <= 2, ]
datatable <- merge(atac_common, bulk.mono_de_common, by = "name")
datatable <- na.omit(datatable)
write.csv(datatable, file = paste0(sannot_dir, "protein_coding_diffs.csv"))

# Organize for plotting
datatable$alpha_level <- ifelse(datatable$isDiff_1 | datatable$isDiff_2, 1, 0.2)

# Correlation between only the differentials
cor_diff <- cor(datatable[datatable$isDiff_1 | datatable$isDiff_2, ]$Log2FC_atac, datatable[datatable$isDiff_1 | datatable$isDiff_2, ]$Log2FC_rna)
# 0.2896502 is the correlation between the differentials

# Plot fold-change vs fold-change plot
p1 <- ggplot(datatable, aes(
  x = .data[["Log2FC_atac"]], y = .data[["Log2FC_rna"]],
  color = interaction(isDiff_1, isDiff_2, sep = "-", lex.order = TRUE),
  alpha = .data[["alpha_level"]]
)) +
  geom_point(size = 2) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "bottom", axis.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  scale_color_manual(paste0("Is differential? "),
    values = c("#a6cee3", "#006400", "#1a1ce3", "#e31a1c"),
    labels = c("Not Differential", "differential in RNA", "differential in ATAC", "differential in both")
  ) +
  scale_alpha_continuous(range = c(0.1, 1), guide = "none") +
  scale_x_continuous(limits = c(min(datatable$Log2FC_atac), max(datatable$Log2FC_atac))) +
  scale_y_continuous(limits = c(min(datatable$Log2FC_rna), max(datatable$Log2FC_rna))) +
  annotate("text", x = max(datatable$Log2FC_atac), y = min(datatable$Log2FC_rna), label = paste("Correlation: ", round(cor_diff, 2)), hjust = 1, vjust = 1)

p1 <- p1 +
  geom_text_repel(
    data = datatable[datatable$isDiff_1 | datatable$isDiff_2, ],
    aes(x = .data[["Log2FC_atac"]], y = .data[["Log2FC_rna"]], label = name),
    color = "black", size = 5,
    box.padding = 0.5,
    segment.color = "black", segment.size = 0.1,
    min.segment.length = 0.2,
    max.overlaps = 15
  )
ggsave(p1, filename = paste0(plot_dir, "protein_coding_diffs.pdf"), width = 10, height = 10, dpi = 300)

######################################################################
# Check the genes from Wilk 2021 are present in the datatable
######################################################################
genes <- c(
  "HLA-E", "CCL2", "CD4", "TNF", "IL6", "IL1B", "TXNIP", "CD74", "HLA-DRA", "HLA-DPA1", "HLA-DRB1",
  "HLA-DRB5", "HLA-DQB1", "EIF1", "ZFP36", "S100A8", "S100A9", "XAF1", "IFI6",
  "OAS1", "OAS2", "IFITM3", "LGALS1", "PLBD1", "S100A12", "FCGR3A", "CCL3",
  "CCL4", "HLA-DPB1", "HLA-DMA", "S100A9", "IRF7", "CXCL2", "ACTB", "CLU",
  "PTMA", "ITM3", "DYF", "MX1", "RNY1", "LGALS1", "PPIA"
)
wimmers_genes <- c(
  "IGKC", "ATF3", "IL1B", "DOCK4", "RBM47", "IFI27", "PLSCR1", "STAT1",
  "OAS1", "DDX60", "PARP9", "DDX58", "TNFSF13B", "APOBEC3A", "SIGLEC1",
  "MX2", "OASL", "EIF2AK2", "OAS3", "IFIH1", "IRF7", "ISG15", "HERC5",
  "RSAD2", "IFIT1", "IFIT3", "IFIT2"
)

genes <- union(genes, wimmers_genes)
genes <- genes[genes %in% datatable$name]

# Add a column to indicate if the gene is in the Wilk 2021 list
datatable$highlight <- ifelse(datatable$name %in% genes, "bold", "plain")

# Plot fold-change vs fold-change plot
p1 <- ggplot(datatable, aes(
  x = .data[["Log2FC_atac"]], y = .data[["Log2FC_rna"]],
  color = interaction(isDiff_1, isDiff_2, sep = "-", lex.order = TRUE),
  alpha = .data[["alpha_level"]]
)) +
  geom_point(size = 2) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "bottom", axis.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  scale_color_manual(paste0("Is differential? "),
    values = c("#a6cee3", "#006400", "#1a1ce3", "#e31a1c"),
    labels = c("Not Differential", "differential in RNA", "differential in ATAC", "differential in both")
  ) +
  scale_alpha_continuous(range = c(0.1, 1), guide = "none") +
  scale_x_continuous(limits = c(min(datatable$Log2FC_atac), max(datatable$Log2FC_atac))) +
  scale_y_continuous(limits = c(min(datatable$Log2FC_rna), max(datatable$Log2FC_rna))) +
  annotate("text", x = max(datatable$Log2FC_atac), y = min(datatable$Log2FC_rna), label = paste("Correlation: ", round(cor_diff, 2)), hjust = 1, vjust = 1)

p1 <- p1 +
  geom_text_repel(
    data = datatable[datatable$isDiff_1 | datatable$isDiff_2, ],
    aes(
      x = .data[["Log2FC_atac"]],
      y = .data[["Log2FC_rna"]],
      label = name,
      fontface = highlight
    ),
    color = "black", size = 5,
    box.padding = 0.5,
    segment.color = "black", segment.size = 0.1,
    min.segment.length = 0.2,
    max.overlaps = 15
  )
ggsave(p1, filename = paste0(plot_dir, "protein_coding_diffs_highlighted.pdf"), width = 10, height = 10, dpi = 300)
######################################################################

wimmers_genes <- wimmers_genes[wimmers_genes %in% datatable$name]
datatable$highlight <- ifelse(datatable$name %in% wimmers_genes, "bold", "plain")

# Plot fold-change vs fold-change plot
p1 <- ggplot(datatable, aes(
  x = .data[["Log2FC_atac"]], y = .data[["Log2FC_rna"]],
  color = interaction(isDiff_1, isDiff_2, sep = "-", lex.order = TRUE),
  alpha = .data[["alpha_level"]]
)) +
  geom_point(size = 2) +
  theme_classic() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 0, linetype = "dashed") +
  theme(legend.position = "bottom", axis.title = element_text(size = 16), legend.text = element_text(size = 14)) +
  scale_color_manual(paste0("Is differential? "),
    values = c("#a6cee3", "#006400", "#1a1ce3", "#e31a1c"),
    labels = c("Not Differential", "differential in RNA", "differential in ATAC", "differential in both")
  ) +
  scale_alpha_continuous(range = c(0.1, 1), guide = "none") +
  scale_x_continuous(limits = c(min(datatable$Log2FC_atac), max(datatable$Log2FC_atac))) +
  scale_y_continuous(limits = c(min(datatable$Log2FC_rna), max(datatable$Log2FC_rna))) +
  annotate("text", x = max(datatable$Log2FC_atac), y = min(datatable$Log2FC_rna), label = paste("Correlation: ", round(cor_diff, 2)), hjust = 1, vjust = 1)

p1 <- p1 +
  geom_text_repel(
    data = datatable[datatable$isDiff_1 | datatable$isDiff_2, ],
    aes(
      x = .data[["Log2FC_atac"]],
      y = .data[["Log2FC_rna"]],
      label = name,
      fontface = highlight
    ),
    color = "black", size = 5,
    box.padding = 0.5,
    segment.color = "black", segment.size = 0.1,
    min.segment.length = 0.2,
    max.overlaps = 15
  )

ggsave(p1, filename = paste0(plot_dir, "protein_coding_diffs_wimmers_highlighted.pdf"), width = 10, height = 10, dpi = 300)

######################################################################
