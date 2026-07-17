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
})

addArchRThreads(threads = 30) 
fig_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/figures/"
outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/"
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
p <- plotEmbedding(
  ArchRProj = project,
  colorBy = "GeneScoreMatrix",
  name = markerGenes,
  embedding = "UMAP_HIV",
  quantCut = c(0.01, 0.95)
)

p2 <- lapply(p, function(x) {
  x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank()
    )
})
pdf(file = paste0(fig_dir, "exhaustion_markers_umap.pdf"), width = 10, height = 10)
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
p_lin2 <- lapply(p_lin, function(x) {
  x + guides(color = FALSE, fill = FALSE) +
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
      axis.text.x = element_blank(), axis.ticks.x = element_blank(),
      axis.text.y = element_blank(), axis.ticks.y = element_blank()
    )
})
pdf(file = paste0(fig_dir, "lineage_treg_markers_umap.pdf"), width = 9, height = 6)
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


#####################################################################
# Add peak matrix and motif annotations
#####################################################################
project$TexClusters <- ifelse(project$ClustersHIV %in% c("C1","C2"), "Tex", project$ClustersHIV)

project <- addGroupCoverages(ArchRProj = project, groupBy = "TexClusters", threads = 30, force = TRUE)
pathToMacs2 <- findMacs2()
project <- addReproduciblePeakSet(
    ArchRProj = project, 
    groupBy = "TexClusters", maxPeaks = 300000,
    pathToMacs2 = pathToMacs2, threads = 30
)
project <- addPeakMatrix(project, threads = 30, force = TRUE)
project <- addBgdPeaks(project, force = TRUE)

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
project <- addMotifAnnotations(ArchRProj = project, motifSet = "cisbp", name = "Motif",force = TRUE)

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