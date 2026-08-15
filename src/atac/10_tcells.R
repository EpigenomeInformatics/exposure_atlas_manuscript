#!/usr/bin/env Rscript

#####################################################################
# 10_tcells.R
# Created on 28-04-2025 by Irem B. Gunduz, written by Bei Wei
# Last modified (revision) on 2025 by I.B.G.
# This script is used to analyze T cells
#
# REVISION NOTES (reviewer response):
#  (R1/R2/R3) Replaced pooled Fisher's exact test on cell counts with a
#             subject-level differential abundance test (donor as unit).
#  (R2/R3)    Replaced FOXP2/FOXP3-only motif matching with all FOXP-family
#             motifs; report as family, not FOXP2-specific.
#  (R2/R3)    Added CD4/CD8A/CD8B and Treg-marker (FOXP3/IL2RA) gene-activity
#             panels to support CD8 identity and exclude Treg contamination.
#####################################################################

## Load Libraries
set.seed(42)
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(ArchR)
  library(ggplot2)
  library(dtplyr)
  library(tidyr)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(chromVARmotifs)
  library(methylTFRAnnotationHg38)   # provides the `altius` archetype annotation
  library(ChrAccR)                   # for the diverging chromVAR colour scheme
  library(RColorBrewer)
})
# addPeakAnnotationsNew() lives in the project utils (same as 04_1_markers.R)
source("/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/utils/archr_utils.R")

addArchRThreads(threads = 30)
fig_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/figures/"
outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/"
# --- HIV T-cell project cache -------------------------------------------------
# Build the subsetted/clustered/peak-called project once and save it as an RDS.
# On later runs we read the RDS and skip the expensive steps (LSI, clustering,
# UMAP, peak calling, chromVAR deviations), so only the plots/analyses re-run.
# NB: the Arrow files under outputDir must remain in place for the RDS to load.
hiv_rds <- file.path(outputDir, "hiv_tcell_project.rds")
build_project <- !file.exists(hiv_rds)

if (build_project) {

echo_full <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)

# Get the cell data from echo
cell_data <- as.data.frame(echo_full@cellColData) %>%
  dplyr::filter(sample_exposure_type %in% c("HIV")) %>%
  dplyr::select(Sample, ClusterCellTypes) %>%
  dplyr::filter(ClusterCellTypes %in% c("T_mem_CD8"))

# Subset echo project for this cells
cellsSample <- rownames(cell_data)
project <- echo_full[cellsSample, ]
project <- subsetCells(echo_full, cellNames = cellsSample)

#project <- echo_full[cellsSample, ]

# add iterative LSI
project <- addIterativeLSI(
  ArchRProj = project,
  useMatrix = "TileMatrix",
  name = "IterativeLSI_hiv",
  iterations = 2,
  clusterParams = list(
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:30,
  saveIterations = TRUE,
  nPlot = 100000
)

# add clusters
project <- addClusters(
  input = project,
  reducedDims = "IterativeLSI_hiv",
  method = "Seurat",
  name = "ClustersHIV",
  resolution = 0.8,
  maxClusters=6,
   force = TRUE
)

project <- addUMAP(
  ArchRProj = project,
  reducedDims = "IterativeLSI_hiv",
  name = "UMAP_HIV",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "euclidean",
  force=TRUE
)

} else {
  message("Loading cached HIV T-cell project: ", hiv_rds)
  project <- readRDS(hiv_rds)
}

# Create a vector for sample file names corresponding to each sample
sample_files <- c(
  "hiv6_fragments.tsv.gz", "hiv12_fragments.tsv.gz", "hiv9_fragments.tsv.gz",
  "hiv8_fragments.tsv.gz", "hiv4_fragments.tsv.gz", "hiv1_fragments.tsv.gz",
  "hiv2_fragments.tsv.gz", "hiv7_fragments.tsv.gz", "hiv3_fragments.tsv.gz",
  "hiv11_fragments.tsv.gz", "hiv10_fragments.tsv.gz", "hiv5_fragments.tsv.gz"
)

# Create a vector for the subject assignment corresponding to each sample
subject_info <- c(
  "sub1", "sub1", "sub1",
  "sub2", "sub2", "sub2",
  "sub3", "sub3", "sub3",
  "sub4", "sub4", "sub4"
)

# Create a vector for time points
time_point <- c(
  "pre", "acute", "chronic",
  "pre", "acute", "chronic",
  "pre", "acute", "chronic",
  "pre", "acute", "chronic"
)

# Add the time point information to the ArchR project
project <- addSampleColData(
  ArchRProj = project,
  data = time_point,
  name = "TimePoint",
  samples = sample_files,
  force = TRUE
)

# Add the subject information to the ArchR project
project <- addSampleColData(
  ArchRProj = project,
  data = subject_info,
  name = "Subject",
  samples = sample_files,
  force = TRUE # Set to TRUE to overwrite if the column already exists
)

df <- getEmbedding(project, embedding = "UMAP_HIV", returnDF = TRUE)
colnames(df) <- c("UMAP1", "UMAP2")

# Create a named vector where names are sample files and values are subjects
sample_to_subject <- c(
  "hiv6_fragments.tsv.gz" = "sub1", "hiv12_fragments.tsv.gz" = "sub1", "hiv9_fragments.tsv.gz" = "sub1",
  "hiv8_fragments.tsv.gz" = "sub2", "hiv4_fragments.tsv.gz" = "sub2", "hiv1_fragments.tsv.gz" = "sub2",
  "hiv2_fragments.tsv.gz" = "sub3", "hiv7_fragments.tsv.gz" = "sub3", "hiv3_fragments.tsv.gz" = "sub3",
  "hiv11_fragments.tsv.gz" = "sub4", "hiv10_fragments.tsv.gz" = "sub4", "hiv5_fragments.tsv.gz" = "sub4"
)
df$Sample <- sapply(strsplit(rownames(df), "#"), `[`, 1)
df$Subject <- sample_to_subject[df$Sample]

# Define the color palette for the subjects
colorPalette <- c(
  "sub1" = "#6A1B9A", # Purple for Subject 1
  "sub2" = "#EA80FC", # Violet for Subject 2
  "sub3" = "#EF7E23", # orange for Subject 3
  "sub4" = "#ffbb00"
) # yellow for Subject 4

# Create the UMAP plot with the specified color palette
umap <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = Subject)) +
  geom_point() +
  scale_color_manual(values = colorPalette) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme_void() +
  theme(legend.position = "bottom")

ggsave(umap, file = paste0(fig_dir, "umap_tcell_subject.pdf"), width = 7, height = 7)

df$Clusters <- project$ClustersHIV
tex_palette <- c(
  "C1" = "#B22222", # Red/Brown (Tex)
  "C2" = "#2C3E50", # Dark Blue (Tex)
  "C3" = "#27AE60", # Green
  "C4" = "#8E44AD", # Purple
  "C5" = "#E67E22", # Orange
  "C6" = "#F1C40F"  # Yellow
)

# UMAP colored by clusters
umap_clusters <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = Clusters)) +
  geom_point() +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme_void() +
  theme(legend.position = "bottom")+
    scale_color_manual(values = tex_palette)
ggsave(umap_clusters, file = paste0(fig_dir, "umap_tcell_clusters.pdf"), width = 7, height = 7)

sample_to_timepoint <- setNames(time_point, sample_files)
df$Sample <- sapply(strsplit(rownames(df), "#"), `[`, 1)
df$TimePoint <- sample_to_timepoint[df$Sample]

# UMAP colored by HIV status with distinct colors
umap_hiv <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = TimePoint)) +
  geom_point(size = 0.8) +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme_void() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c(
    "pre" = "#1F78B4",   # Distinct Blue (Baseline)
    "acute" = "#E31A1C", # Vivid Red (Peak/Acute)
    "chronic" = "#33A02C" # Distinct Green (Chronic)
  ))

ggsave(umap_hiv, file = paste0(fig_dir, "umap_tcell_timepoint.pdf"), width = 7, height = 7)


#####################################################################
# HIV Stage/exhaustion associations
#####################################################################

# Define Status
df$Status <- ifelse(df$Clusters %in% c("C1","C2"), "Tex", "Other")

## <<< REVISION (R1/R2/R3): subject-level differential abundance ---------------
## Previously, Tex accumulation was tested with a Fisher's exact test on a
## pooled cell-count contingency table (df$Status x df$TimePoint). That treats
## individual cells as independent observations and ignores donor structure
## (pseudoreplication). All three reviewers flagged this. We instead compute
## the per-subject Tex proportion at each stage and test the change across the
## four subjects, so the donor is the unit of analysis.

# Per-subject, per-timepoint Tex proportion
subj_prop <- df %>%
  dplyr::group_by(Subject, TimePoint) %>%
  dplyr::summarise(
    n_tex   = sum(Status == "Tex"),
    n_total = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(tex_prop = n_tex / n_total)

# Wide format: one row per subject, one column per stage
subj_prop_wide <- subj_prop %>%
  dplyr::select(Subject, TimePoint, tex_prop) %>%
  tidyr::pivot_wider(names_from = TimePoint, values_from = tex_prop)

write.csv(subj_prop_wide,
          file = paste0(fig_dir, "tex_proportion_per_subject.csv"),
          row.names = FALSE)

# Paired, subject-level tests across stages (N = 4 subjects).
# Wilcoxon signed-rank is used as the nonparametric paired test; with n = 4 the
# p-values are necessarily coarse, so we report them alongside the per-subject
# values rather than as standalone significance. We also report the paired
# t-test for reference.
run_subject_level_paired <- function(wide, g1, g2) {
  x <- wide[[g1]]
  y <- wide[[g2]]
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]; y <- y[keep]
  wt <- suppressWarnings(wilcox.test(x, y, paired = TRUE))
  tt <- suppressWarnings(t.test(x, y, paired = TRUE))
  data.frame(
    Comparison   = paste(g1, "vs", g2),
    n_subjects   = length(x),
    mean_diff    = round(mean(y - x), 4),
    Wilcoxon_P   = signif(wt$p.value, 3),
    Paired_t_P   = signif(tt$p.value, 3)
  )
}

subject_level_results <- rbind(
  run_subject_level_paired(subj_prop_wide, "pre",   "acute"),
  run_subject_level_paired(subj_prop_wide, "acute", "chronic"),
  run_subject_level_paired(subj_prop_wide, "pre",   "chronic")
)

write.csv(subject_level_results,
          file = paste0(fig_dir, "tex_subject_level_test.csv"),
          row.names = FALSE)

# Per-subject Tex-proportion trajectory (visual support for the subject-level test)
subj_prop_plot <- subj_prop %>%
  dplyr::mutate(Stage = factor(TimePoint, levels = c("pre", "acute", "chronic")))

tex_traj <- ggplot(subj_prop_plot,
                   aes(x = Stage, y = tex_prop, color = Subject, group = Subject)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  labs(x = "Stage", y = "Tex proportion (per subject)") +
  theme_classic() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = colorPalette)

ggsave(tex_traj, file = paste0(fig_dir, "tex_proportion_trajectory.pdf"),
       width = 6, height = 5)

## <<< REVISION (R2/R3): donor composition of each Tex cluster -----------------
## Reviewer (R2, annotation reproducibility across individuals) and co-author
## comments (Fabian 64/65/66) asked for the percentage of cells in each
## exhausted (Tex) cluster contributed by each HIV subject, to show that both
## Tex clusters (C1, C2) are made up of cells from all four subjects rather than
## being driven by a single individual. This is the ClustersHIV x Subject
## cross-tabulation (as percentages) promised in the response letter
## (Response Figure N).
cluster_levels_all <- paste0("C", 1:6)   # C1/C2 are the exhausted (Tex) clusters
donor_comp_df <- getCellColData(project, select = c("ClustersHIV", "Sample")) %>%
  as.data.frame() %>%
  dplyr::mutate(Subject = sample_to_subject[as.character(Sample)])

# Within each cluster: percentage of its cells contributed by each subject
donor_comp <- donor_comp_df %>%
  dplyr::count(ClustersHIV, Subject, name = "n_cells") %>%
  tidyr::complete(ClustersHIV = cluster_levels_all, Subject = names(colorPalette),
                  fill = list(n_cells = 0)) %>%
  dplyr::group_by(ClustersHIV) %>%
  dplyr::mutate(pct = 100 * n_cells / sum(n_cells)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    Subject     = factor(Subject, levels = names(colorPalette)),
    ClustersHIV = factor(ClustersHIV, levels = cluster_levels_all)
  )

write.csv(donor_comp,
          file = paste0(fig_dir, "hiv_cluster_donor_composition.csv"),
          row.names = FALSE)

# Horizontal stacked barplot: one bar per cluster (all six), split by subject
# (each bar sums to 100%). Clusters on the y-axis (C1 at top), proportion on x.
donor_comp_bar <- ggplot(donor_comp,
                         aes(x = pct, y = ClustersHIV, fill = Subject)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = colorPalette) +
  scale_y_discrete(limits = rev(cluster_levels_all)) +   # C1 at top
  labs(x = "Cells per cluster (%)", y = "HIV CD8+ T-cell cluster (C1, C2 = Tex)",
       fill = "Subject") +
  theme_classic(base_size = 10) +
  theme(legend.position = "right")

ggsave(donor_comp_bar,
       file = paste0(fig_dir, "hiv_cluster_donor_composition_bar.pdf"),
       width = 6, height = 4)
## <<< END REVISION -----------------------------------------------------------

## NOTE: the previous pooled Fisher test is retained below but commented out,
## so the change from the original submission is transparent. Do not use for
## the manuscript claim.
# status_time_table <- table(df$Status, df$TimePoint)
# run_pairwise_fisher <- function(tbl) {
#   stages <- colnames(tbl)
#   pairs <- combn(stages, 2, simplify = FALSE)
#   results <- data.frame()
#   for (pair in pairs) {
#     g1 <- pair[1]; g2 <- pair[2]
#     current_matrix <- tbl[, c(g1, g2)]
#     ft <- fisher.test(current_matrix)
#     results <- rbind(results, data.frame(
#       Comparison = paste(g1, "vs", g2),
#       Odds_Ratio = round(ft$estimate, 2),
#       P_Value = ft$p.value))
#   }
#   results$P_Adj <- p.adjust(results$P_Value, method = "bonferroni")
#   results
# }
# pairwise_results <- run_pairwise_fisher(status_time_table)
# write.csv(pairwise_results, file = paste0(fig_dir, "fisher_test_results.csv"), row.names = FALSE)
## <<< END REVISION -----------------------------------------------------------

###########################################################################################

# Subject | Pre | Acute | Chronic
df <- data.frame(
  Subject = c("S1", "S2", "S3", "S4"),
  Pre = c(0, 0, 40, 40),
  Acute = c(1866453, 852580, 415480, 27250),
  Chronic = c(794, 110691, 244551, 12341)
)


# Reshape the data for ggplot and order the 'Stage' factor levels
df_long <- df %>%
  pivot_longer(cols = -Subject, names_to = "Stage", values_to = "Viral_Load") %>%
  mutate(Stage = factor(Stage, levels = c("Pre", "Acute", "Chronic")))


# Plotting
lp <- ggplot(df_long, aes(x = Stage, y = Viral_Load, color = Subject, group = Subject)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  labs(
    title = "Viral Load per Subject Across Stages",
    x = "Stage",
    y = "Viral Load"
  ) +
  theme_classic() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("S1" = "#1f77b4", "S2" = "#ff7f0e", "S3" = "#7f7f7f", "S4" = "#ffbb00"))
ggsave(lp, file = paste0(fig_dir, "viral_load_subject.pdf"), width = 7, height = 7)


## <<< REVISION (R2/R3): CD8 identity + Treg-exclusion markers -----------------
## Reviewers asked (a) how CD4 vs CD8 T cells were distinguished, and (b) that
## the possibility of CD4 Treg / activated-CD4 contamination in the Tex clusters
## be excluded, since CTLA4 is also a Treg marker. We add lineage markers
## (CD4, CD8A, CD8B) and canonical Treg markers (FOXP3, IL2RA) to the exhaustion
## panel so the Tex clusters can be shown to be CD8+ and FOXP3-/IL2RA-low.

# Exhaustion markers (original) + lineage and Treg markers (added)
markerGenes <- c(
  "HAVCR2", "CTLA4", "NCAM1", "ROBO2", "ROBO1", "TIGIT",  # exhaustion (original)
  "CD8A", "CD8B", "CD4",                                   # lineage (added)
  "FOXP3", "IL2RA"                                         # Treg (added)
)
## <<< END REVISION -----------------------------------------------------------

# Get the marker genes
markersGS <- getMarkerFeatures(
  ArchRProj = project,
  useMatrix = "GeneScoreMatrix",
  groupBy = "ClustersHIV",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1.25")

# Create the heatmap
heatmapGS <- markerHeatmap(
  seMarker = markersGS,
  cutOff = "FDR <= 1 & Log2FC >= 0.25",
  labelMarkers = markerGenes,
  transpose = TRUE
)

pdf(file = paste0(fig_dir, "exhaustion_markers_heatmap.pdf"), width = 10, height = 10)
draw(heatmapGS)
dev.off()

project <- addImputeWeights(project)

## <<< REVISION (R2/R3): keep a compact colour scale on each UMAP panel --------
## Previously each gene-activity UMAP dropped its colour guide
## (guides(color = FALSE, fill = FALSE)), so the panels had no min-to-max scale.
## Reviewers/co-authors need to read the gene-activity range, so we keep a small
## per-panel colour bar (each gene keeps its own scale, since e.g. CD8A and FOXP3
## span different ranges). Axis ticks stay off to keep the panels clean.
style_umap_gene <- function(x) {
  x +
    theme_ArchR(baseSize = 6.5) +
    theme(
      plot.margin       = unit(c(0.05, 0.05, 0.05, 0.05), "cm"),
      axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y  = element_blank(), axis.ticks.y = element_blank(),
      legend.position   = "right",
      legend.key.width  = unit(0.20, "cm"),
      legend.key.height = unit(0.45, "cm"),
      legend.text       = element_text(size = 5),
      legend.title      = element_text(size = 5)
    )
}
## <<< END REVISION -----------------------------------------------------------

p <- plotEmbedding(
  ArchRProj = project,
  colorBy = "GeneScoreMatrix",
  name = markerGenes,
  embedding = "UMAP_HIV",
  quantCut = c(0.01, 0.95)
)

p2 <- lapply(p, style_umap_gene)
pdf(file = paste0(fig_dir, "exhaustion_markers_umap.pdf"), width = 12, height = 10)
do.call(cowplot::plot_grid, c(list(ncol = 3), p2))
dev.off()

## <<< REVISION (R2/R3): dedicated lineage/Treg UMAP panel ---------------------
## A separate, clearly labelled panel for the CD8-identity / Treg-exclusion
## markers, so it can be cited directly in the annotation-justification text.
lineage_treg_genes <- c("CD8A", "CD8B", "CD4", "FOXP3", "IL2RA")
p_lin <- plotEmbedding(
  ArchRProj = project,
  colorBy = "GeneScoreMatrix",
  name = lineage_treg_genes,
  embedding = "UMAP_HIV",
  quantCut = c(0.01, 0.95)
)
p_lin2 <- lapply(p_lin, style_umap_gene)
pdf(file = paste0(fig_dir, "lineage_treg_markers_umap.pdf"), width = 11, height = 6)
do.call(cowplot::plot_grid, c(list(ncol = 3), p_lin2))
dev.off()

# Quantitative check: mean gene-activity of lineage/Treg markers in Tex vs Other,
# to state numerically that Tex clusters are CD8+ and FOXP3-/IL2RA-low.
gsm <- getMatrixFromProject(project, useMatrix = "GeneScoreMatrix")
gs_rownames <- rowData(gsm)$name
want <- intersect(lineage_treg_genes, gs_rownames)
gs_sub <- assay(gsm)[match(want, gs_rownames), , drop = FALSE]
rownames(gs_sub) <- want
status_vec <- ifelse(project$ClustersHIV %in% c("C1","C2"), "Tex", "Other")
tex_summary <- sapply(want, function(g) {
  tapply(gs_sub[g, ], status_vec, mean)
})
tex_summary <- t(tex_summary)  # rows = gene, cols = Tex/Other
write.csv(as.data.frame(tex_summary),
          file = paste0(fig_dir, "tex_lineage_treg_meanactivity.csv"))
## <<< END REVISION -----------------------------------------------------------


## <<< REVISION (R2): per-cluster lineage + exhaustion activity distributions ---
## Side-by-side violin/box plots of the lineage markers (CD8A, CD8B, CD4) and the
## two canonical exhaustion markers (CTLA4, HAVCR2) across the six HIV CD8
## subclusters. Shows that all clusters remain CD8A/CD8B-positive and CD4-low,
## and that CTLA4/HAVCR2 gene activity is elevated in the Tex clusters (C1/C2).
## This is the quantitative, per-cluster companion to the marker UMAPs and shows
## the CD8A gradient within the compartment.
comp_genes_wanted <- c("CD8A", "CD8B", "CD4", "CTLA4", "HAVCR2")
# gs_sub above holds only the lineage/Treg panel, so the exhaustion markers are
# taken straight from the full gene-score matrix instead.
comp_genes <- intersect(comp_genes_wanted, gs_rownames)
comp_mat   <- assay(gsm)[match(comp_genes, gs_rownames), , drop = FALSE]
rownames(comp_mat) <- comp_genes

# align gene-activity columns to the project cell order before attaching clusters
cell_ids         <- colnames(gsm)
clusters_by_cell <- as.character(project$ClustersHIV)[match(cell_ids, project$cellNames)]

comp_df <- do.call(rbind, lapply(comp_genes, function(g) {
  data.frame(Gene = g, Cluster = clusters_by_cell, Activity = comp_mat[g, ])
})) %>%
  dplyr::mutate(
    Cluster = factor(Cluster, levels = paste0("C", 1:6)),
    Gene    = factor(Gene, levels = comp_genes)
  )

comp_violin <- ggplot(comp_df, aes(x = Cluster, y = Activity, fill = Cluster)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.2) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", linewidth = 0.2) +
  facet_wrap(~ Gene, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = tex_palette) +
  labs(x = "HIV CD8+ T-cell cluster (C1, C2 = Tex)", y = "Gene activity") +
  theme_classic(base_size = 10) +
  theme(legend.position = "none")

ggsave(comp_violin,
       file = paste0(fig_dir, "hiv_cluster_CD4_CD8_geneactivity_violin.pdf"),
       width = 14, height = 4)

# Per-cluster mean gene activity table (for the legend/text)
comp_means <- comp_df %>%
  dplyr::group_by(Gene, Cluster) %>%
  dplyr::summarise(mean_activity = round(mean(Activity, na.rm = TRUE), 3), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Cluster, values_from = mean_activity)
write.csv(comp_means,
          file = paste0(fig_dir, "hiv_cluster_CD4_CD8_geneactivity_means.csv"),
          row.names = FALSE)
## <<< END REVISION -----------------------------------------------------------


#####################################################################
# Add peak matrix and motif annotations
#####################################################################
project$TexClusters <- ifelse(project$ClustersHIV %in% c("C1","C2"), "Tex", project$ClustersHIV)

if (build_project) {
  project <- addGroupCoverages(ArchRProj = project, groupBy = "TexClusters", threads = 30, force = TRUE)
  pathToMacs2 <- findMacs2()
  project <- addReproduciblePeakSet(
      ArchRProj = project,
      groupBy = "TexClusters", maxPeaks = 300000,
      pathToMacs2 = pathToMacs2, threads = 30
  )
  project <- addPeakMatrix(project, threads = 30, force = TRUE)
  project <- addBgdPeaks(project, force = TRUE)
}

#####################################################################
# Marker testing
#####################################################################
# Rename clusters as Tex and others
project$Status <- ifelse(project$ClustersHIV %in% c("C1","C2"), "Tex", "Other")

markerTest <- getMarkerFeatures(
  ArchRProj = project, 
  useMatrix = "PeakMatrix",
  groupBy = "Status", 
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = c("Tex"),
  bgdGroups = "Other"
)

pv <- plotMarkers(seMarker = markerTest, name = "Tex", cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", plotAs = "Volcano")
ggsave(pv, file = paste0(fig_dir, "marker_volcano_c1.pdf"), width = 6, height = 5)
if (build_project) project <- addMotifAnnotations(ArchRProj = project, motifSet = "cisbp", name = "Motif", force = TRUE)

motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = project,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggsave(ggUp, file = paste0(fig_dir, "motif_enrichment_tex.pdf"), width = 6, height = 5)

motifsDo <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = project,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
  )

df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))

ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

ggsave(ggDo, file = paste0(fig_dir, "motif_depletion_tex.pdf"), width = 6, height = 5)

#####################################################################

data("human_pwms_v2")
motifs <- human_pwms_v2

## <<< REVISION (R2/R3): match the full FOXP family, not FOXP2/FOXP3 only ------
## Reviewers noted that FOXP-family motifs are highly similar and effectively
## indistinguishable, and that singling out FOXP2 is misleading (FOXP2 is not a
## canonical T-cell TF; FOXP3 is the canonical Treg factor). We therefore match
## ALL FOXP-family motifs and display them as a family in the browser track,
## rather than restricting to FOXP2/FOXP3. Individual members are still captured
## for completeness but interpreted at the family level in the text.
foxp <- names(motifs)[grepl("FOXP", names(motifs), ignore.case = TRUE)]
message("FOXP-family motifs matched: ", paste(foxp, collapse = ", "))

motif_positions <- motifmatchr::matchMotifs(
  pwms = motifs[foxp],
  subject = getPeakSet(project),
  genome = "hg38",
  out = "matches"
)
diff_peaks_gr <- getMarkers(markerTest, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", returnGR = TRUE)$Tex
saveRDS(diff_peaks_gr, file = paste0(outputDir, "diff_peaks_tex.rds"))

# Collapse all FOXP-family motif matches into a single family-level peak set,
# so the browser track shows "FOXP family" rather than one arbitrary member.
foxp_cols <- seq_len(ncol(motif_positions))
foxp_any  <- Reduce(`|`, lapply(foxp_cols, function(j) assay(motif_positions)[, j]))
foxp_family_peaks <- getPeakSet(project)[foxp_any]

peaks_matched <- list(
  cCREs       = getPeakSet(project),
  FOXP_family = foxp_family_peaks,  # all FOXP members, union of matches
  DiffPeaks   = diff_peaks_gr
)
## <<< END REVISION -----------------------------------------------------------

genes <- c("CTLA4","HAVCR2")  # removed FOXP2 as a highlighted "gene"; family shown in track
plot <- plotBrowserTrack(
  ArchRProj = project,
  features = peaks_matched,
  minCells = 50,
  tileSize = 500,
  geneSymbol = genes,
  ylim = c(0.001, 0.999),
  upstream = 75000,
  downstream = 75000,
  groupBy = "TexClusters",
  pal = c("#B22222", "#2C3E50", "#27AE60", "#8E44AD","#E67E22", "#F1C40F")
)

plotPDF(
  plotList = plot,
  name = "Plot-Motif-HIV-CD8T-FOXPfamily.pdf",
  ArchRProj = project,
  addDOC = FALSE, width = 5, height = 5
)

#####################################################################
# TF-activity (chromVAR) UMAPs for the supplement
# (Fabian, comments 7/9/11: include the TF activities we previously only
#  checked manually in a supplementary panel). We compute per-cell motif
#  deviations and project the FOXP-family and exhaustion-associated TF
#  activities onto the T-cell UMAP so they can be shown directly.
#####################################################################
if (build_project) project <- addDeviationsMatrix(ArchRProj = project, peakAnnotation = "Motif", force = TRUE)

## getFeatures() on a MotifMatrix returns BOTH the bias-corrected z-scores
## ("z:TF_###") and the raw deviations ("deviations:TF_###"). The previous grep
## matched both, so half the panels were raw-deviation UMAPs that render
## near-flat/blank. We keep only the z-score features (the standard chromVAR
## readout) and match each TF symbol exactly (anchored after "z:", so "TOX"
## does not also grab TOX2/TOX3).
motif_features <- getFeatures(project, useMatrix = "MotifMatrix")
motif_z <- grep("^z:", motif_features, value = TRUE)

pick <- function(tf) {
  hits <- grep(paste0("^z:", tf, "(_|$)"), motif_z, value = TRUE, ignore.case = TRUE)
  if (length(hits) == 0L) message("  [TFactivity] no z-score motif for ", tf)
  hits
}
# FOXP family (the response focus) + canonical CD8 exhaustion / memory TFs
tf_wanted <- c("FOXP1", "FOXP2", "FOXP3", "TOX", "TBX21", "EOMES", "TCF7", "NR4A1")
tf_show   <- unique(unlist(lapply(tf_wanted, pick)))
message("TF-activity UMAP motifs: ", paste(tf_show, collapse = ", "))
stopifnot(length(tf_show) > 0)

p_tf <- plotEmbedding(
  ArchRProj     = project,
  colorBy       = "MotifMatrix",
  name          = tf_show,
  embedding     = "UMAP_HIV",
  imputeWeights = getImputeWeights(project)
)
if (!is.list(p_tf)) p_tf <- list(p_tf)

# Same compact colour scale as the gene-activity panels, so the z-score range shows
p_tf2 <- lapply(p_tf, style_umap_gene)
pdf(file = paste0(fig_dir, "TFactivity_UMAP_tcells.pdf"), width = 12, height = 10)
do.call(cowplot::plot_grid, c(list(ncol = 3), p_tf2))
dev.off()

# ArchR's native multi-page version kept for completeness
plotPDF(p_tf,
  name = "TFactivity_UMAP_tcells_archr.pdf",
  ArchRProj = project, addDOC = FALSE, width = 5, height = 5
)

#####################################################################
# chromVAR TF-activity heatmap across HIV clusters (C1-C6), Altius archetypes
# ArchR-native marker heatmap of the most differentially active (variable)
# Altius archetypes across the six HIV CD8+ T-cell clusters, on chromVAR
# z-scores. Coloured with the ChrAccR diverging scheme (`.default.div`,
# teal-to-brown) so it matches the cvar TF-deviation heatmap.
#####################################################################
# Add the Altius archetype annotation + deviations matrix on this T-cell peak set
if (build_project) {
  # subsetCells() inherits peak annotations from the full project, whose match
  # matrices are sized to the OLD peak set. After building a new T-cell peak set
  # a stale 'altius' annotation makes addDeviationsMatrix fail with
  # 'subscript out of bounds', so drop it and recompute against the new peaks.
  if ("altius" %in% names(project@peakAnnotation)) {
    project@peakAnnotation <- project@peakAnnotation[
      setdiff(names(project@peakAnnotation), "altius")
    ]
  }
  altius <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/altius_tf_bindsites.rds")
  project <- addPeakAnnotationsNew(ArchRProj = project,
                                   regions = altius, name = "altius")
  project <- addDeviationsMatrix(ArchRProj = project, peakAnnotation = "altius", force = TRUE)

  # Cache the fully-built project (all matrices) so later runs skip the heavy steps
  saveRDS(project, file = hiv_rds)
  message("Saved HIV T-cell project to ", hiv_rds)
}

# Top-30 most VARIABLE Altius archetypes (chromVAR combined variability), shown
# as mean per-cluster z-scores. The differential marker test returned too few
# archetypes on this small T-cell peak set, so we rank by variability instead.
vd <- as.data.frame(getVarDeviations(project, name = "altiusMatrix", plot = FALSE))
all_feats <- if ("name" %in% colnames(vd)) vd$name else rownames(vd)

n_top   <- 50
top_var <- head(all_feats, n_top)

# Always include the FOXP/FOX-family archetype(s), even if not among the top
# variable set, since they are the focus of the HIV/Tex analysis.
fox_feats <- grep("fox", all_feats, value = TRUE, ignore.case = TRUE)
if (length(fox_feats)) message("FOX-family archetypes added: ", paste(fox_feats, collapse = ", "))
top_var <- unique(c(top_var, fox_feats))

# Per-cell z-scores averaged per HIV cluster, keeping the variability/append order
mm   <- getMatrixFromProject(project, useMatrix = "altiusMatrix")
zmat <- assays(mm)[["z"]]
rownames(zmat) <- rowData(mm)$name
top_var <- top_var[top_var %in% rownames(zmat)]
zmat <- zmat[top_var, , drop = FALSE]

clust          <- project$ClustersHIV
cluster_levels <- paste0("C", 1:6)
mean_z <- sapply(cluster_levels, function(cl)
  rowMeans(zmat[, clust == cl, drop = FALSE], na.rm = TRUE))
mean_z <- t(scale(t(mean_z)))          # row-scale for cluster-to-cluster contrast

# Same diverging palette as the cvar TF-deviation heatmap (ChrAccR .default.div,
# teal-to-brown); fall back to a BrBG teal-brown ramp if the config is absent.
colors.cv <- tryCatch(
  ChrAccR::getConfigElement("colorSchemesCont")[[".default.div"]],
  error = function(e) rev(RColorBrewer::brewer.pal(11, "BrBG"))
)
col_fun <- circlize::colorRamp2(
  seq(-2, 2, length.out = length(colors.cv)), colors.cv)

hm_motifs <- ComplexHeatmap::Heatmap(
  mean_z,
  name            = "Row z\n(mean chromVAR)",
  col             = col_fun,
  cluster_columns = FALSE,
  column_order    = cluster_levels,
  show_row_names  = TRUE,
  row_names_gp    = grid::gpar(fontsize = 6),
  column_names_gp = grid::gpar(fontsize = 10)
)

pdf(file = paste0(fig_dir, "hiv_chromVAR_cluster_heatmap_altius.pdf"), width = 6, height = 11)
ComplexHeatmap::draw(hm_motifs)
dev.off()

#####################################################################
## <<< REVISION: lineage + exhaustion marker activity across ALL T-cell types --
## The panel above is restricted to the HIV CD8 subclusters. This is the same
## marker set (CD8A, CD8B, CD4, CTLA4, HAVCR2) across the T-cell compartment of
## the whole atlas, with the annotated cell type on the x axis instead of the HIV
## subcluster, so the exhaustion markers can be read against the lineage markers
## in every T-cell population rather than only in the HIV subset.
##
## T_mix is excluded: it is the mixed/unresolved T-cell cluster, so its marker
## distributions are a blend of the populations either side of it and would be
## read as intermediate biology rather than as an annotation artefact.
#####################################################################

tcell_types <- c("T_naive", "T_mem_CD4", "T_mem_CD8", "T_mait") # T_mix excluded

# Canonical atlas cell-type colours (same values as 13_cvar_analysis.R and
# 14_l2fc.R), so this panel is comparable with the rest of the figures.
tcell_palette <- c(
  "T_naive"   = "#C7E9B4",
  "T_mem_CD4" = "#4292c6",
  "T_mem_CD8" = "#0074cc",
  "T_mait"    = "#41B6C4"
)

# `project` above is the HIV-only subset, so reload the full atlas if the build
# branch did not already put it in the session.
if (!exists("echo_full")) {
  echo_full <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)
}

cd_all <- as.data.frame(getCellColData(echo_full, select = "ClusterCellTypes"))
t_cells <- rownames(cd_all)[cd_all$ClusterCellTypes %in% tcell_types]
message("T-cell compartment: ", length(t_cells), " cells across ",
        length(tcell_types), " cell types (T_mix excluded)")
tproj <- ArchR::subsetCells(echo_full, cellNames = t_cells)

# Subset the cells FIRST, so only the T-cell compartment is pulled into memory.
gsm_t   <- getMatrixFromProject(tproj, useMatrix = "GeneScoreMatrix")
gs_t_rn <- rowData(gsm_t)$name
# Reuse the marker set from the HIV panel so the two figures stay identical;
# fall back to the literal list if only this section is being re-run.
marker_genes_t <- if (exists("comp_genes_wanted")) comp_genes_wanted else
  c("CD8A", "CD8B", "CD4", "CTLA4", "HAVCR2")
genes_t <- intersect(marker_genes_t, gs_t_rn)
missing_t <- setdiff(marker_genes_t, genes_t)
if (length(missing_t)) {
  message("  not in the gene-score matrix, dropped: ",
          paste(missing_t, collapse = ", "))
}
mat_t <- assay(gsm_t)[match(genes_t, gs_t_rn), , drop = FALSE]
rownames(mat_t) <- genes_t

# align the gene-activity columns to the cell-type labels
ct_by_cell <- as.character(
  getCellColData(tproj, select = "ClusterCellTypes", drop = TRUE)
)[match(colnames(gsm_t), tproj$cellNames)]

tcell_df <- do.call(rbind, lapply(genes_t, function(g) {
  data.frame(Gene = g, CellType = ct_by_cell, Activity = mat_t[g, ])
})) %>%
  dplyr::mutate(
    CellType = factor(CellType, levels = tcell_types),
    Gene     = factor(Gene, levels = genes_t)
  )

## Gene-activity scores are zero-inflated with long right tails: nearly all the
## mass sits at zero and a few cells reach 50. On a linear axis with scale =
## "width", the drawn shape is set by those few extreme cells and everything
## informative is compressed into the bottom of the panel. Two readable views:

## ---- (a) Dot plot: detection rate and mean level ----------------------------
## The standard readout for zero-inflated single-cell data. Dot SIZE is the
## percentage of cells with non-zero activity, dot COLOUR is the mean activity
## scaled within each gene (so genes on different absolute scales are
## comparable). This separates "how many cells have the marker at all" from
## "how strong it is where present", which a violin conflates.
tcell_dot <- tcell_df %>%
  dplyr::group_by(Gene, CellType) %>%
  dplyr::summarise(
    n_cells      = dplyr::n(),
    pct_detected = 100 * mean(Activity > 0, na.rm = TRUE),
    mean_all     = mean(Activity, na.rm = TRUE),
    mean_detected = mean(Activity[Activity > 0], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(Gene) %>%
  dplyr::mutate(mean_scaled = as.numeric(scale(mean_all))) %>%
  dplyr::ungroup()

p_dot <- ggplot(tcell_dot, aes(x = CellType, y = Gene)) +
  geom_point(aes(size = pct_detected, colour = mean_scaled)) +
  scale_size_continuous(range = c(1, 9), name = "% cells\ndetected") +
  scale_colour_gradient2(low = "#4292c6", mid = "grey92", high = "#B2182B",
    midpoint = 0, name = "Mean activity\n(z within gene)") +
  scale_y_discrete(limits = rev(levels(tcell_df$Gene))) +
  labs(x = "T-cell type", y = NULL) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_line(colour = "grey94", linewidth = 0.3))

ggsave(p_dot,
       file = paste0(fig_dir, "tcell_lineage_exhaustion_geneactivity_dotplot.pdf"),
       width = 6.5, height = 4)

## ---- (b) Violin on a log scale ----------------------------------------------
## Same distributions as before, but log1p-transformed so the zero-inflated bulk
## is resolved instead of being flattened by the tail. Kept as the companion to
## the dot plot for anyone who wants the full distribution rather than a summary.
# Linear gene-activity units, matching every other gene-activity panel in the
# manuscript, with a free y scale so each gene gets a range suited to its own
# distribution.
#
# A per-gene trim is still needed: within a single gene the bulk of cells sit
# near zero while a handful reach 30-50, and on a linear axis those few cells set
# the range and flatten everything else (this is what the first version of this
# panel looked like). The violin layer is therefore drawn from the cells at or
# below each gene's 99.5th percentile, while the median and IQR are computed from
# EVERY cell -- so no reported statistic is affected by the trim, only the shape
# of the drawn tail.
viol_caps <- tcell_df %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(cap = stats::quantile(Activity, 0.995, na.rm = TRUE),
    .groups = "drop")
viol_df <- tcell_df %>%
  dplyr::left_join(viol_caps, by = "Gene") %>%
  dplyr::filter(Activity <= cap)
message("Violin panel: trimmed ",
  nrow(tcell_df) - nrow(viol_df), " of ", nrow(tcell_df),
  " cell-gene values above the per-gene 99.5th percentile (drawing only)")

tcell_violin <- ggplot(mapping = aes(x = CellType, y = Activity)) +
  geom_violin(data = viol_df, aes(fill = CellType),
    scale = "width", trim = TRUE, linewidth = 0.2) +
  # A narrow crossbar rather than a boxplot: the old white box was wider than
  # most of the violins it sat inside, so the box read as the data.
  stat_summary(data = tcell_df, fun = median, geom = "crossbar",
    width = 0.45, linewidth = 0.25, colour = "grey15", fatten = 0) +
  stat_summary(data = tcell_df,
    fun.min = function(z) stats::quantile(z, 0.25, na.rm = TRUE),
    fun.max = function(z) stats::quantile(z, 0.75, na.rm = TRUE),
    fun = median, geom = "linerange", linewidth = 0.5, colour = "grey15") +
  facet_wrap(~ Gene, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = tcell_palette) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  labs(x = "T-cell type", y = "Gene activity",
    caption = paste("Bar, median; line, interquartile range, both over all cells.",
      "Violin outlines omit the top 0.5% of cells per gene so the axis is not set by a few extreme cells.")) +
  theme_classic(base_size = 10) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.caption = element_text(size = 7, hjust = 0))

ggsave(tcell_violin,
       file = paste0(fig_dir, "tcell_lineage_exhaustion_geneactivity_violin.pdf"),
       width = 12, height = 4)

## ---- Numbers behind both panels ---------------------------------------------
## NB the unit here is the CELL. Every cell is one observation, so these are
## descriptive only -- a formal comparison between cell types would be
## pseudoreplicated across the cells of a donor and must be done at donor level.
write.csv(
  tcell_dot %>%
    dplyr::mutate(dplyr::across(c(pct_detected, mean_all, mean_detected,
      mean_scaled), ~ round(.x, 3))),
  file = paste0(fig_dir, "tcell_lineage_exhaustion_geneactivity_summary.csv"),
  row.names = FALSE
)

# wide per-cell-type means, kept for the figure legend
tcell_means <- tcell_df %>%
  dplyr::group_by(Gene, CellType) %>%
  dplyr::summarise(mean_activity = round(mean(Activity, na.rm = TRUE), 3),
                   n_cells = dplyr::n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = CellType,
                     values_from = c(mean_activity, n_cells))
write.csv(tcell_means,
          file = paste0(fig_dir, "tcell_lineage_exhaustion_geneactivity_means.csv"),
          row.names = FALSE)
## <<< END REVISION -----------------------------------------------------------

#####################################################################

