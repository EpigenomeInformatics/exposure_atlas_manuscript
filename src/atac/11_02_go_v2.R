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
datatable_pos <- paste0(sannot_dir, "protein_coding_datatable_pos.csv")
datatable_pos <- read.csv(datatable_pos, row.names = 1)
datatable_neg <- paste0(sannot_dir, "protein_coding_datatable_neg.csv")
datatable_neg <- read.csv(datatable_neg, row.names = 1)
pos_genes <- rownames(datatable_pos)
neg_genes <- rownames(datatable_neg)

# Convert to Entrez IDs
pos_genes <- bitr(pos_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
neg_genes <- bitr(neg_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

# Run GO enrichment analysis
pos_go_results <- enrichGO(
  gene          = pos_genes$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",            # Use "MF" or "CC" for other ontologies
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.01,
  readable      = TRUE
)

# Run GO enrichment analysis
neg_go_results <- enrichGO(
  gene          = neg_genes$E∂NTREZID,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",            # Use "MF" or "CC" for other ontologies
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.01,
  readable      = TRUE
)

# Generate the plot for odds ratio (FoldEnrichment) with coloring based on p-value
pdf(paste0(plot_dir, "GO_enrichment_CD14_mon∂o_severe_vs_control_pos_odds_ratio_colored.pdf"), width = 10, height = 10)
pos_go_results %>%
    mutate(OddsRatio = FoldEnrichment) %>%  # Use FoldEnrichment as OddsRatio
    #filter(p.adjust < pvalue_threshold) %>%   # Apply p-value threshold
    arrange(desc(OddsRatio)) %>%            # Sort by OddsRatio in descending order
    head(20) %>%                            # Select top 20 GO terms
    ggplot(aes(x = reorder(Description, OddsRatio), y = OddsRatio, fill = qvalue)) +  # Color by -log10(pvalue)
    geom_bar(stat = "identity") +           # Create bar plot
    coord_flip() +                          # Flip coordinates for better readability
    scale_fill_gradient(low = "blue", high = "red", name = "qvalue") +  # Gradient color scale
    labs(
        title = "GO Enrichment: CD14 Mono Severe vs Control (Odds Ratio)",
        x = "GO Terms",
        y = "Fold Enrichment"
    ) +
    theme_minimal()
dev.off()

# Generate the plot for odds ratio (FoldEnrichment) with coloring based on p-value
pdf(paste0(plot_dir, "GO_enrichment_CD14_mono_severe_vs_control_neg_odds_ratio_colored.pdf"), width = 10, height = 10)
neg_go_results %>%
    mutate(OddsRatio = FoldEnrichment) %>%  # Use FoldEnrichment as OddsRatio
    #filter(p.adjust < pvalue_threshold) %>%   # Apply p-value threshold
    arrange(desc(OddsRatio)) %>%            # Sort by OddsRatio in descending order
    head(20) %>%                            # Select top 20 GO terms
    ggplot(aes(x = reorder(Description, OddsRatio), y = OddsRatio, fill = qvalue)) +  # Color by -log10(pvalue)
    geom_bar(stat = "identity") +           # Create bar plot
    coord_flip() +                          # Flip coordinates for better readability
    scale_fill_gradient(low = "blue", high = "red", name = "qvalue") +  # Gradient color scale
    labs(
        title = "GO Enrichment: CD14 Mono Severe vs Control (Odds Ratio)",
        x = "GO Terms",
        y = "Fold Enrichment"
    ) +
    theme_minimal()
dev.off()
