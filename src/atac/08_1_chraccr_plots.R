#!/usr/bin/env Rscript

#####################################################################
# 08_1_chraccr_plots.R
# created on 2023-10-02 by Irem Gunduz
# ChrAccR plots for COVID samples
#####################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(ChrAccR)
  library(ggplot2)
  library(data.table)
  library(stringr)
  library(muRtools)
  library(ggplot2)
  library(circlize)
  library(ComplexHeatmap)
  library(ggrepel)
  library(GenomicRanges)
  library(muLogR)
  library(clue)
})
set.seed(12)
source("/scratch/icbb/igunduz/irem_github/exposure_atlas_manuscript/utils/lola.R")
source("/scratch/icbb/igunduz/irem_github/exposure_atlas_manuscript/utils/chraccr_plots.R")
source("/icbb/projects/igunduz/sc_epigenome_exp/utils/helpers.R")

color_mapping <- c(
  "C19_ctrl" = "#6CA6CD",
  "C19_mild" = "#8FBC8F",
  "C19_mod" = "#2C948F",
  "C19_sev" = "#006400",
  "HIV_acu" = "#893368",
  "HIV_chr" = "#825CA3",
  "HIV_ctrl" = "#4F619D",
  "Influenza_ctrl" = "#4F619D",
  "Influenza_d3" = "#FFD34E",
  "Influenza_d6" = "#EBB332",
  # day 28 is the true final timepoint; "d30" is a legacy mislabel
  "Influenza_d30" = "#E1861A",
  "Influenza_d28" = "#E1861A",
  "OP_low" = "#BADBF4",
  "OP_med" = "#94A9D3",
  "OP_high" = "#16528A"
)

diffCompNames <- c(
  "HIV_ctrl vs HIV_chr",
  "HIV_ctrl vs HIV_acu",
  "Influenza_d3 vs Influenza_ctrl",
  "Influenza_d30 vs Influenza_ctrl",
  "Influenza_d6 vs Influenza_ctrl",
  "OP_high vs OP_low",
  "OP_high vs OP_med",
  "OP_med vs OP_low"
)
cells <- c("B_mem", "B_naive", "Mono_CD14", "Mono_CD16", "NK_CD16", "T_mem_CD8", "T_mem_CD4", "T_naive", "T_mix", "T_naive", "T_mait")
outputDir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_2023-10-02/"
plotDir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/plots/"
if (!dir.exists(plotDir)) {
  dir.create(plotDir)
}
for (cell in cells) {
  # Iterate over each comparison in the comparison table
  comparisonTable <- readRDS(paste0(outputDir, cell, "/reports/differential_data/comparisonTable.rds"))
  comparisonTable$comps <- paste0(comparisonTable$grp1, " vs ", comparisonTable$grp2)
  not_present <- diffCompNames[which(!diffCompNames %in% comparisonTable$comps)]

  for (i in 1:nrow(comparisonTable)) {
    # Get the relevant information for the current comparison
    grp1 <- comparisonTable$grp1[i]
    grp2 <- comparisonTable$grp2[i]
    comp <- comparisonTable$comps[i]
    pp <- plotMAwithChrAccR(cell, outputDir, i, "archrPeaks")
    filename <- paste0(plotDir, "MA_plot_", cell, "_", grp1, " vs ", grp2, ".pdf")
    ggsave(filename = filename, plot = pp, width = 20, height = 20)
  }
}

outputDir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/"
diffCompNames <- c(
  "C19_mild vs C19_ctrl",
  "C19_mod vs C19_ctrl",
  "C19_sev vs C19_ctrl"
)
cells <- c("Mono_CD14", "NK_CD16", "T_mem_CD8", "T_mem_CD4", "T_naive")
plotDir <- "/icbb/projects/igunduz/finalize_echo_050824/MA_plots_251124/"
if (!dir.exists(plotDir)) {
  dir.create(plotDir)
}
for (cell in cells) {
  # Iterate over each comparison in the comparison table
  comparisonTable <- readRDS(paste0(outputDir, cell, "/reports/differential_data/comparisonTable.rds"))
  comparisonTable$comps <- paste0(comparisonTable$grp1, " vs ", comparisonTable$grp2)
  not_present <- diffCompNames[which(!diffCompNames %in% comparisonTable$comps)]

  for (i in 1:nrow(comparisonTable)) {
    # Get the relevant information for the current comparison
    grp1 <- comparisonTable$grp1[i]
    grp2 <- comparisonTable$grp2[i]
    comp <- comparisonTable$comps[i]
    pp <- plotMAwithChrAccR(cell, outputDir, i, "archrPeaks")
    filename <- paste0(plotDir, "MA_plot_", cell, "_", grp1, " vs ", grp2, ".pdf")
    ggsave(filename = filename, plot = pp, width = 20, height = 20)
  }
}

###############################################################################
# l2fc plots for COVID-MONO

# archr_peaks dataframe to the annotations
global_peaklist <- readRDS("/icbb/projects/igunduz/DARPA_analysis/artemis_031023/rawData/regionsGR.rds") %>%
  data.table::as.data.table() %>%
  dplyr::select(!width & !strand) %>%
  as.data.frame()

outputDir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/"
cell <- "Mono_CD14"
plot_dir <- "/icbb/projects/igunduz/finalize_echo_050824/C19/"
if (!dir.exists(plot_dir)) {
  dir.create(plot_dir)
}
cor_list <- list()
# ATAC get the differentially accessble regions
atac <- prepareDARforPlot(cell, outputDir, 3, "archrPeaks") # get the C19_severe
for (comp in diffCompNames[1:2]) {
  if (!file.exists(paste0(plot_dir, "L2FC_", cell, "_", comp, ".pdf"))) {
    j <- ifelse(comp == "C19_mild vs C19_ctrl", 1, 2)

    # Methylation, get the differntiall methylated regions
    atac2 <- prepareDARforPlot(cell, outputDir, j, "archrPeaks")
    # dmt <-  prepareDMRforPlot(cell,comp,path,global_peaklist,"sites")

    # filter methylation by overlapped regions of atac
    datatable <- prepareScatterMethAtac(atac, atac2, 4)
    cor_list[[comp]] <- cor.test(datatable$log2FC_1, datatable$log2FC_2)
    # plot the scatter plot
    scplot <- plotScatterL2FC(datatable, y_lab = paste0("Covid-19 ", j), x_lab = "Covid-19 Severe", comb = paste0("Covid-19 "), group1 = "log2FC_1", group2 = "log2FC_2")

    # save the plot
    ggsave(plot = scplot, paste0(plot_dir, "L2FC_", cell, "_", comp, ".pdf"), width = 20, height = 20)
    ggsave(plot = scplot, paste0(plot_dir, "L2FC_", cell, "_", comp, ".png"), width = 20, height = 20)
  }
}
saveRDS(cor_list, file = paste0(plot_dir, "cor_list.rds"))

###############################################################################

# load loladb database
lolaDb <- "/icbb/projects/share/annotations/lolaDB/hg38/"
lolaDb <- RnBeads::loadLolaDbs(lolaDb)
# lola_path <- "/icbb/projects/igunduz/DARPA_analysis/LOLA/"
lola_path <- "/icbb/projects/igunduz/finalize_echo_050824/C19/"
if (!dir.exists(lola_path)) {
  dir.create(lola_path)
}

# volcano plot for atac
p <- lolaVolcanoPlotC19("Mono_CD14", lolaDb, outputDir, pValCut = 1.5, region = "archrPeaks", database = "TF_motif_clusters", signifCol = "qValue", cnames = c("loss", "gain"), n = 3, colorpanel = c())$plot
ggsave(plot = p, paste0(lola_path, "Mono_atac_TF_motif_clusters.pdf"), width = 20, height = 20)

p <- lolaVolcanoPlotC19("Mono_CD14", lolaDb, outputDir, pValCut = 1.5, region = "archrPeaks", database = "TF_motifs", signifCol = "qValue", cnames = c("loss", "gain"), n = 3, colorpanel = c())$plot
ggsave(plot = p, paste0(lola_path, "Mono_atac_TF_motifs.pdf"), width = 20, height = 20)

p <- lolaVolcanoPlotC19("Mono_CD14", lolaDb, outputDir, pValCut = 1.5, region = "archrPeaks", database = "TF_motifs_vert", signifCol = "qValue", cnames = c("loss", "gain"), n = 3, colorpanel = c())$plot
ggsave(plot = p, paste0(lola_path, "Mono_atac_TF_motifs_vert.pdf"), width = 20, height = 20)

###############################################################################
# differential heatmap for COVID-MONO
logger.start("Starting differential heatmap for COVID-MONO")
diff <- prepareDARforPlot(cell, outputDir, 3, "archrPeaks") # , padj = 0.05) # get the C19_severe
# diff <- diff[isDiff == TRUE,] #filter the differentially accessible regions
diff$Start <- diff$Start + 1
diff$id <- paste0(diff$Chromosome, ":", diff$Start, "-", diff$End)
diffm <- diff
diff <- dplyr::select(diff, c("id", "isDiff"))

ds <- ChrAccR::loadDsAcc(paste0(outputDir, cell, "/data/dsATAC_processed/"))
logger.info("Getting fragment GR")
counts <- getCounts(ds, "archr_peaks")
# Get the counts matrix
coords <- data.table::as.data.table(getCoord(ds, "archr_peaks")) %>%
  dplyr::mutate(id = paste0(seqnames, ":", start, "-", end)) %>%
  dplyr::select(id)
# Apply k-means clustering
counts <- data.table::as.data.table(computeZScore(counts))

counts <- cbind(counts, coords)
counts <- merge(counts, diff)
counts <- counts[isDiff == TRUE, ]
counts <- dplyr::select(counts, !id & !isDiff)
assignments <- consensus_kmeans(counts, 5, 1)
diffm[diffm$isDiff == TRUE, "cluster"] <- assignments
diffm <- diffm[diffm$isDiff == TRUE, ]

# Reorder columns in the counts matrix and column annotations
counts_ordered <- setcolorder(counts, sort(names(counts)))

# Create column annotations based on the sample names
col_ann <- data.frame(Exposure = gsub("^([^_]*_[^_]*).*", "\\1", colnames(counts_ordered)))
rownames(col_ann) <- colnames(counts)
colnames(col_ann) <- "Exposure"

# Create row annotations based on cluster assignments
row_ann <- data.frame(Cluster = as.factor(assignments))
row_ann_ordered <- row_ann[order(row_ann$Cluster), ]

# counts_ordered$cluster_assignments <- assignments # add cluster info
setDT(counts_ordered)
rownames(counts_ordered) <- diffm$id

col_ha <- HeatmapAnnotation(
  Exposure = col_ann$Exposure,
  col = list(Exposure = c(
    "C19_ctrl" = "#4F609C",
    "C19_mild" = "#C43E96",
    "C19_mod" = "#06948E",
    "C19_sev" = "#C03830"
  ))
)

# Define the order of the columns based on their groups
group_order <- c("C19_ctrl", "C19_mild", "C19_mod", "C19_sev")
column_order <- colnames(counts_ordered)
column_order <- column_order[order(match(column_order, group_order))]

dend_cols <- cluster_within_group(counts_ordered, col_ann$Exposure)

# Create heatmap
ht <- Heatmap(
  counts_ordered,
  column_order = column_order,
  cluster_columns = FALSE,
  col = ArchR::paletteContinuous("solarExtra"),
  name = "Accessibility Z-score",
  cluster_rows = TRUE,
  show_row_names = FALSE,
  show_column_names = FALSE,
  heatmap_legend_param = list(at = c(-max(round(counts_ordered)), 0, max(round(counts_ordered)))),
  column_title_side = "top",
  top_annotation = col_ha
)


pdf(paste0("/icbb/projects/igunduz/pheatmap_", cell, "_samplewise_C19_diff_all.pdf"), width = 20, height = 20)
draw(ht)
dev.off()


###############################################################################
# Add label to peaks
###############################################################################

ctrl_columns <- grep("C19_ctrl", colnames(counts_ordered), value = TRUE)
sev_columns <- grep("C19_sev", colnames(counts_ordered), value = TRUE)
ctrl_means <- rowMeans(counts_ordered[, ..ctrl_columns])
sev_means <- rowMeans(counts_ordered[, ..sev_columns])
labels <- ifelse(ctrl_means > sev_means, "control", "severe")

labels <- as.data.frame(labels)
rownames(labels) <- rownames(counts_ordered)
saveRDS(labels, file = "/icbb/projects/igunduz/mono_projection_070824/labels.rds")


###############################################################################
# R3.4: moderate vs severe DAR overlap, sign concordance, one-group-only
# Uses prepareDARforPlot + isDiff filter (same DAR definition as Fig 3).
###############################################################################

# severe = index 3, moderate = index 2 (matches the j<-ifelse logic above)
sev_tab <- prepareDARforPlot(cell, outputDir, 3, "archrPeaks")   # C19_sev vs ctrl
mod_tab <- prepareDARforPlot(cell, outputDir, 2, "archrPeaks")   # C19_mod vs ctrl

# match the ID convention used elsewhere in this script (Start + 1)
mk_id <- function(d) {
  d$Start <- d$Start + 1
  paste0(d$Chromosome, ":", d$Start, "-", d$End)
}
sev_tab$id <- mk_id(sev_tab)
mod_tab$id <- mk_id(mod_tab)

# DARs = isDiff == TRUE (same filter as the heatmap section)
sevDAR <- sev_tab[sev_tab$isDiff == TRUE, ]
modDAR <- mod_tab[mod_tab$isDiff == TRUE, ]

cat("moderate DARs:", nrow(modDAR), " (expect ~373)\n")
cat("severe DARs:  ", nrow(sevDAR), " (expect ~923)\n")

# overlap and one-group-only
overlap  <- intersect(modDAR$id, sevDAR$id)
mod_only <- setdiff(modDAR$id, sevDAR$id)
sev_only <- setdiff(sevDAR$id, modDAR$id)

# sign concordance among shared DARs
lfc <- "log2FC"
m_lfc <- modDAR[[lfc]][match(overlap, modDAR$id)]
s_lfc <- sevDAR[[lfc]][match(overlap, sevDAR$id)]
concordant <- sign(m_lfc) == sign(s_lfc)

summary_r34 <- data.frame(
  moderate_DARs       = nrow(modDAR),
  severe_DARs         = nrow(sevDAR),
  shared              = length(overlap),
  moderate_only       = length(mod_only),
  severe_only         = length(sev_only),
  pct_sign_concordant = round(100 * mean(concordant), 1),
  shared_lfc_pearson  = round(cor(m_lfc, s_lfc, use = "complete.obs"), 3)
)
print(summary_r34)
write.csv(summary_r34,
          file = paste0(plot_dir, "moderate_vs_severe_DAR_summary.csv"),
          row.names = FALSE)

###############################################################################
# R3.4 (revision): moderate-vs-severe DAR overlap, split by direction.
# The direct moderate-vs-severe comparison is computed in
# 07_2_run_ChrAccR_C19.R and read from mod_vs_sev_DAR.rds.
# Direction convention follows ChrAccR: log2FC > 0 = higher accessibility in the
# disease group of each comparison.
###############################################################################
fig_dir <- "/scratch/icbb/igunduz/irem_github/exposure_atlas_manuscript/figures/"

# Select the ChrAccR comparisons BY NAME (robust to diffTab ordering) rather
# than by a hard-coded positional index; moderate-vs-severe comes from the
# direct comparison in 07_2_run_ChrAccR_C19.R (mod_vs_sev_DAR.rds).
comparisonTable <- readRDS(paste0(outputDir, cell, "/reports/differential_data/comparisonTable.rds"))
comp_index <- function(grp1, grp2) {
  i <- which(comparisonTable$grp1 == grp1 & comparisonTable$grp2 == grp2)
  if (length(i) != 1) stop("comparison not found in comparisonTable: ", grp1, " vs ", grp2)
  i
}

get_dar <- function(index) {
  tab <- prepareDARforPlot(cell, outputDir, index, "archrPeaks")
  tab$Start <- tab$Start + 1
  tab$id <- paste0(tab$Chromosome, ":", tab$Start, "-", tab$End)
  tab[tab$isDiff == TRUE, c("id", "log2FC")]
}

dar_cm <- get_dar(comp_index("C19_mod", "C19_ctrl")) # moderate vs control
dar_cs <- get_dar(comp_index("C19_sev", "C19_ctrl")) # severe vs control

ms_tab <- readRDS(paste0(plot_dir, "mod_vs_sev_DAR.rds")) # from 07_2 (R3.4)
dar_ms <- ms_tab[ms_tab$isDiff == TRUE, c("id", "log2FC")] # severe vs moderate

# report set sizes split by direction so empty sets are explicit
message(sprintf(
  "DAR counts (up/down)  ctrl-mod: %d/%d | ctrl-sev: %d/%d | mod-sev: %d/%d",
  sum(dar_cm$log2FC > 0), sum(dar_cm$log2FC < 0),
  sum(dar_cs$log2FC > 0), sum(dar_cs$log2FC < 0),
  sum(dar_ms$log2FC > 0), sum(dar_ms$log2FC < 0)
))

###############################################################################
# Single compact panel replacing the two UpSet plots.
# One stacked bar carries everything the upsets did: the three overlap
# categories (moderate-only / shared / severe-only), the hyper- vs hypo-
# accessible split within each, and -- in the subtitle -- the effect-size
# concordance of the shared core plus the null result of the direct
# moderate-vs-severe test. Small enough to sit as an inset in Figure 3
# (Fabian's suggestion) rather than taking a full panel.
###############################################################################
ov       <- intersect(dar_cm$id, dar_cs$id)
mod_only <- setdiff(dar_cm$id, dar_cs$id)
sev_only <- setdiff(dar_cs$id, dar_cm$id)

# direction taken from the group in which the region is differential
# (severe for the shared set; the shared set is direction-concordant anyway)
dir_of <- function(ids, tab) {
  ifelse(tab$log2FC[match(ids, tab$id)] > 0, "Hyper-accessible", "Hypo-accessible")
}
cat_df <- rbind(
  data.frame(Category = "Moderate\nonly", Direction = dir_of(mod_only, dar_cm)),
  data.frame(Category = "Shared",         Direction = dir_of(ov,       dar_cs)),
  data.frame(Category = "Severe\nonly",   Direction = dir_of(sev_only, dar_cs))
)
cat_df$Category <- factor(cat_df$Category,
  levels = c("Moderate\nonly", "Shared", "Severe\nonly")
)

# effect-size concordance within the shared core
m_lfc_ov <- dar_cm$log2FC[match(ov, dar_cm$id)]
s_lfc_ov <- dar_cs$log2FC[match(ov, dar_cs$id)]
r_shared <- cor(m_lfc_ov, s_lfc_ov, use = "complete.obs")
pct_conc <- 100 * mean(sign(m_lfc_ov) == sign(s_lfc_ov))

p_overlap <- ggplot(cat_df, aes(x = Category, fill = Direction)) +
  geom_bar(colour = "grey20", width = 0.7) +
  geom_text(stat = "count", aes(label = after_stat(count)),
    position = position_stack(vjust = 0.5), size = 3.2, colour = "white"
  ) +
  scale_fill_manual(values = c(
    "Hyper-accessible" = "#E64B35", "Hypo-accessible" = "#3C5488"
  )) +
  labs(
    title = "CD14+ monocyte DARs: moderate vs severe COVID-19",
    subtitle = sprintf(
      "Shared: %.0f%% direction-concordant, effect sizes r = %.2f | direct moderate-vs-severe test: %d DARs",
      pct_conc, r_shared, nrow(dar_ms)
    ),
    x = NULL, y = "Number of DARs", fill = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(legend.position = "top", plot.subtitle = element_text(size = 8))

ggsave(paste0(fig_dir, "DAR_moderate_vs_severe_overlap.pdf"), p_overlap,
  width = 5.2, height = 4.2
)
message("Wrote compact DAR-overlap panel to ", fig_dir)

