#!/usr/bin/env Rscript

#####################################################################
# 10_GO_analysis.R
# Created by: IBG on 23-05-2025
#####################################################################


## Load Libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(Seurat)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(DOSE)

})
set.seed(12) # set seed

# Load the data and set the output directory
plot_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/figures/"
sannot_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/sample_annots/"
bulk.mono.de <- readRDS(paste0(sannot_dir, "gene_expression_protein_coding_diffs.rds"))


bulk.mono_de_common$isDiff_2 <- ifelse(bulk.mono_de_common$FDR_rna <= 0.05 & abs(bulk.mono_de_common$Log2FC_rna) > 0.5, T, F)

# Filter significant genes
sig_genes <- subset(bulk.mono.de, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
gene_symbols <- sig_genes$name

# Convert to Entrez IDs
gene_df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Run GO enrichment analysis
go_results <- enrichGO(
  gene          = gene_df$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",            # Use "MF" or "CC" for other ontologies
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.01,
  readable      = TRUE
)

pdf(paste0(plot_dir, "GO_enrichment_CD14_mono_severe_vs_control.pdf"), width = 10, height = 10)
mutate(go_results, qscore = -log(p.adjust, base=100)) %>% 
    barplot(x="qscore",        # Number of top GO terms to display
  title = "GO Enrichment: CD14 Mono Severe vs Control"
)
dev.off()


gop <- go_results %>%
  top_n(25, wt = -p.adjust) %>%
  mutate(Description = factor(Description, levels = rev(Description))) %>%
  ggplot(aes(x = Description, y = Count)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(title = "GO Enrichment: CD14 Mono Severe vs Control")

ggsave(paste0(plot_dir, "GO_enrichment_CD14_mono_severe_vs_control_top20.pdf"), plot = gop, width = 10, height = 8)

# KEGG enrichment
kegg_res <- enrichKEGG(
  gene         = gene_df$ENTREZID,
  organism     = 'hsa',          # human
  pvalueCutoff = 0.01
)

pdf(paste0(plot_dir, "KEGG_enrichment_CD14_mono_severe_vs_control.pdf"), width = 10, height = 10)
barplot(kegg_res, showCategory = 50, title = "KEGG Pathways: CD14 Mono Severe vs Control")
dev.off()

# Convert to data frame and sort by adjusted p-value
kegg_df <- as.data.frame(kegg_res)
kegg_df$log_padj <- -log10(kegg_df$p.adjust)
kegg_top25 <- kegg_df[order(kegg_df$p.adjust), ][1:25, ]

# Plot and save
pdf(paste0(plot_dir, "KEGG_top25_sorted_by_pval.pdf"), width = 10, height = 8)
ggplot(kegg_top25, aes(x = reorder(Description, log_padj), y = log_padj)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 25 KEGG Pathways (Sorted by p.adj)",
    x = "Pathway",
    y = "-log10(adj. p-value)"
  ) +
  theme_minimal(base_size = 12)
dev.off()

#########################################
library(DOSE)
library(enrichplot)

sig_genes <- subset(bulk.mono.de, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
gene_symbols <- sig_genes$name

# Convert to Entrez IDs
gene_df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
edo <- enrichDGN(as.vector(gene_df$ENTREZID))

pdf(paste0(plot_dir, "GO_enrichment_CD14_mono_severe_vs_control.pdf"), width = 10, height = 10)
barplot(edo, showCategory=20) 
dev.off()

pdf(paste0(plot_dir, "GO_enrichment_CD14_mono_severe_vs_control.pdf"), width = 10, height = 10)
mutate(edo, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore")
dev.off()


edo2 <- gseDO(geneList)
pdf(paste0(plot_dir, "ORA_enrichment_CD14_mono_severe_vs_control.pdf"), width = 10, height = 10)
dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
dev.off()
pdf(paste0(plot_dir, "GSEA_enrichment_CD14_mono_severe_vs_control.pdf"), width = 10, height = 10)
dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")
dev.off()

edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')

p1 <- heatplot(edox, showCategory=5)
p2 <- heatplot(edox, foldChange=geneList, showCategory=5)
pdf(paste0(plot_dir, "heatmap_enrichment_CD14_mono_severe_vs_control.pdf"), width = 10, height = 10)
cowplot::plot_grid(p1, p2, ncol=1, labels=LETTERS[1:2],
                   label_size=20, align='v', axis='tb',
                   rel_heights=c(1, 1.2))
dev.off()

##################################################################################
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
project_mono <- project_mono[project_mono$sample_exposure_group %in% c("C19_sev", "C19_mild"), ]
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
  bgdGroups = "C19_mild",
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
#write.csv(datatable, file = paste0(sannot_dir, "protein_coding_diffs.csv"))

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
ggsave(p1, filename = paste0(plot_dir, "protein_coding_diffs_mild.pdf"), width = 10, height = 10, dpi = 300)


# Subset the genes that are differential in both ATAC and RNA
datatable_diff <- datatable[datatable$isDiff_1 & datatable$isDiff_2, ]$name
saveRDS(datatable_diff, file = paste0(sannot_dir, "protein_coding_diffs_mild.rds"))

datatable_diff <- readRDS(paste0(sannot_dir, "protein_coding_diffs_mild.rds"))
  library(clusterProfiler)
  library(org.Hs.eg.db)
# Convert to Entrez IDs
gene_df <- bitr(datatable_diff, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Run GO enrichment analysis
go_results <- enrichGO(
  gene          = gene_df$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",            # Use "MF" or "CC" for other ontologies
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.01,
  readable      = TRUE
)


# Sıralama: Count değerine göre azalan sırada
go_results@result <- go_results@result[order(go_results@result$Count, decreasing = TRUE), ]

# PDF çıktısı için görselleştirme
pdf(paste0(plot_dir, "GO_enrichment_CD14_mono_severe_vs_mild.pdf"), width = 12, height = 8)

# Dotplot ile görselleştirme
dotplot(go_results, 
        showCategory = 20, 
        title = "GO Enrichment: CD14 Mono Severe vs Mild") +
  theme_minimal()

# Alternatif: Barplot ile görselleştirme
barplot(go_results, 
        showCategory = 20, 
        title = "GO Enrichment: CD14 Mono Severe vs Mild") +
  theme_minimal()

dev.off()
