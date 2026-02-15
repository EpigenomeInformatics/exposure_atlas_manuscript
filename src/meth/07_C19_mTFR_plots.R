#!/usr/bin/env Rscript

##########################################################################################################################################
# sub_hm_090425.R
# Created by Irem Gunduz on 27-04-2025
# This script generates heatmaps and scatter plots for the methylation data per monocyte subtypes
##########################################################################################################################################

suppressPackageStartupMessages({
  library(methylTFR)
  library(circlize)
  library(ComplexHeatmap)
  library(ggplot2)
  library(ggrepel)
  library(muLogR)
  library(data.table)
  library(rlist)
})

set.seed(42)

# --- Configuration & Directories ---
mtfr_dir <- "/icbb/projects/igunduz/Figure_5_040425/covid_mtfr_070425/jaspar2020_distal/"
motifset <- "jaspar2020"
id <- "C19_mono2_vs_sev_090425"
condition <- c("C19_ctrl", "C19_sev")
plot_dir <- "/icbb/projects/igunduz/exposure_atlas_manuscript/src/meth_processing/Fig_5_sub_080425/plots/"
ds_dir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/Mono_CD14/data/dsATAC_filtered/"

if (!dir.exists(plot_dir)) {
  dir.create(plot_dir, recursive = TRUE)
}

logger.start("Starting Analysis Pipeline")

# --- Helper Functions ---
outliers <- c(
  "Monocyte_CoV_S_S8_D1.bedGraph.bed",
  "Monocyte_CoV_S_S15_D7.bedGraph.bed",
  "Monocyte_Ctrl_1_F_White_45yo.bedGraph.bed",
  "Monocyte_Ctrl_12_F_White_32yo.bedGraph.bed",
  "Monocyte_Ctrl_11_M_White_57yo.bedGraph.bed"
)

mono_1 <- c(
  "Monocyte_CoV_S_S15_D1.bedGraph.bed",
  "Monocyte_CoV_S_S11_D3.bedGraph.bed",
  "Monocyte_CoV_S_S7_D1.bedGraph.bed",
  "Monocyte_CoV_S_S11_D1.bedGraph.bed"
)

get_group <- function(filename) {
  if (filename %in% mono_1) {
    return("C19_sev_mono1")
  } else if (grepl("CoV", filename)) {
    return("C19_sev_mono2")
  } else if (grepl("Ctrl", filename)) {
    return("C19_ctrl")
  } else {
    return("Unknown")
  }
}

##################################################################################
# 1. Plot PCA
##################################################################################
logger.info("Generating PCA Plots")

# Read methylTFR deviations
mtfr_devs <- list.files(mtfr_dir, pattern = "Monocyte", full.names = TRUE)
mtfr_devs <- list.cbind(lapply(condition, function(x) {
  methylTFR::deviations(readRDS(mtfr_devs[grepl(x, mtfr_devs)]))
}))
mtfr_devs <- as.data.frame(mtfr_devs)

# Remove outliers and compute Z-score
mtfr_devs <- mtfr_devs[, !colnames(mtfr_devs) %in% outliers]
mtfr_devs_z <- methylTFR:::computeRowZScore(as.matrix(mtfr_devs))

# Get groups
groups <- unlist(lapply(colnames(mtfr_devs_z), get_group))

# Clean column names
colnames(mtfr_devs_z) <- gsub(".bedGraph.bed", "", colnames(mtfr_devs_z))
colnames(mtfr_devs_z) <- gsub("Monocyte_", "", colnames(mtfr_devs_z))

# --- PCA 1: All Groups ---
tdf <- as.data.frame(t(mtfr_devs_z))
pca_result <- prcomp(tdf, center = TRUE, scale. = TRUE)
tdf$groups <- groups
sample_names <- rownames(tdf)

pdf(paste0(plot_dir, "PCA_mTFR_allmotifs_no_outliers.pdf"), width = 10, height = 10)
autoplot(pca_result, data = tdf, colour = "groups", main = "PCA", size = 5) +
  geom_text_repel(aes(label = sample_names), size = 3) + 
  theme_classic() +
  theme(legend.position = "bottom")
dev.off()

# --- PCA 2: Mono1 vs Ctrl Only ---
subset_idx <- groups %in% c("C19_ctrl", "C19_sev_mono1")
mtfr_devs_z_sub <- mtfr_devs_z[, subset_idx]
groups_sub <- groups[subset_idx]

tdf <- as.data.frame(t(mtfr_devs_z_sub))
pca_result <- prcomp(tdf, center = TRUE, scale. = TRUE)
tdf$groups <- groups_sub
sample_names <- rownames(tdf)

pdf(paste0(plot_dir, "PCA_mTFR_allmotifs_no_outliers_mono1.pdf"), width = 10, height = 10)
autoplot(pca_result, data = tdf, colour = "groups", main = "PCA", size = 5) +
  geom_text_repel(aes(label = sample_names), size = 3) + 
  theme_classic() +
  theme(legend.position = "bottom")
dev.off()


##################################################################################
# 2. Plot Heatmaps
##################################################################################
logger.info("Generating Heatmaps")

# Reload Data for Heatmaps (Ensuring Clean State)
mtfr_devs <- list.files(mtfr_dir, pattern = "Monocyte", full.names = TRUE)
mtfr_devs <- list.cbind(lapply(condition, function(x) {
  methylTFR::deviations(readRDS(mtfr_devs[grepl(x, mtfr_devs)]))
}))
mtfr_devs <- as.data.frame(mtfr_devs)

# Filter Outliers and Compute Z-scores
mtfr_devs <- mtfr_devs[, !colnames(mtfr_devs) %in% outliers]
mtfr_devs_z <- methylTFR:::computeRowZScore(as.matrix(mtfr_devs))
motifs <- readRDS(paste0(ds_dir, "cvar_motifs_diff_080425.rds"))
mtfr_devs_z <- mtfr_devs_z[motifs, ]

groups <- unlist(lapply(colnames(mtfr_devs_z), get_group))
valid_groups <- c("C19_ctrl", "C19_sev_mono1", "C19_sev_mono2")
mtfr_devs_z <- mtfr_devs_z[, groups %in% valid_groups]
groups <- groups[groups %in% valid_groups]

# --- Heatmap 1: All 3 Groups (Ctrl, Mono1, Mono2) ---
ann <- data.frame(Exposure = groups)
rownames(ann) <- colnames(mtfr_devs_z)
group_colors <- c("C19_ctrl" = "#0072B2", "C19_sev_mono1" = "#D55E00", "C19_sev_mono2" = "#009E73")

column_ha <- HeatmapAnnotation(
  df = ann,
  col = list(Exposure = group_colors),
  annotation_legend_param = list(title = "Exposure")
)
c <- muRtools::colpal.cont(nrow(ann), "cptcity.arendal_temperature") 

pdf(paste0(plot_dir, "mTFR_allmotifs_no_outliers_mono2_CLUSTERED.pdf"), width = 20, height = 20)
draw(Heatmap(
  mtfr_devs_z,
  cluster_rows = TRUE,
  cluster_columns = TRUE, 
  column_split = ann$Exposure, 
  top_annotation = column_ha, 
  col = c, 
  heatmap_legend_param = list(title = "methylTFR"),
  show_column_names = FALSE
))
dev.off()

# --- Heatmap 2: Ctrl vs Mono1 Only ---
# Subset again for specific comparison
subset_idx <- groups %in% c("C19_ctrl", "C19_sev_mono1")
mtfr_devs_z_sub <- mtfr_devs_z[, subset_idx]
groups_sub <- groups[subset_idx]

colnames(mtfr_devs_z_sub) <- gsub(".bedGraph.bed", "", colnames(mtfr_devs_z_sub))
colnames(mtfr_devs_z_sub) <- gsub("Monocyte_", "", colnames(mtfr_devs_z_sub))

ann <- data.frame(Exposure = groups_sub)
rownames(ann) <- colnames(mtfr_devs_z_sub)
group_colors <- c("C19_ctrl" = "#0072B2", "C19_sev_mono1" = "#D55E00")

column_ha <- HeatmapAnnotation(
  df = ann,
  col = list(Exposure = group_colors),
  annotation_legend_param = list(title = "Exposure")
)
c <- muRtools::colpal.cont(nrow(ann), "cptcity.arendal_temperature") 

heatmap_obj <- Heatmap(
  mtfr_devs_z_sub,
  cluster_rows = TRUE,
  cluster_columns = TRUE, 
  column_split = ann$Exposure, 
  top_annotation = column_ha, 
  col = c, 
  heatmap_legend_param = list(title = "methylTFR"),
  show_column_names = TRUE
)

pdf(paste0(plot_dir, "mTFR_allmotifs_no_outliers_mono1_CLUSTERED.pdf"), width = 20, height = 20)
draw(heatmap_obj)
dev.off()

# Save Row Order
row_order <- row_order(heatmap_obj)
motifs_ordered <- rownames(mtfr_devs_z_sub)[row_order]
saveRDS(motifs_ordered, paste0(plot_dir, "motifs_ordered.rds"))


# --- Heatmap 3: chromVAR (using mTFR order) ---
logger.info("Plotting chromVAR heatmap")
chromvar_mat <- readRDS(paste0(ds_dir, "chromvar_jaspar2020_101224_corrected_zscores.rds"))
rownames(chromvar_mat) <- sub(".*_", "", rownames(chromvar_mat))
zmat <- chromvar_mat

# Clamp Z-scores
zmat[zmat > 5] <- 5
zmat[zmat < -5] <- -5

# Subset to same motifs and order
zmat <- zmat[motifs_ordered, ]
ann <- data.table(Exposure = ifelse(grepl("ctrl", colnames(zmat)), "C19_ctrl", "C19_sev"))
rownames(ann) <- colnames(zmat)

column_ha <- HeatmapAnnotation(
  Exposure = ann$Exposure, 
  col = list(Exposure = c("C19_ctrl" = "blue", "C19_sev" = "red")), 
  annotation_legend_param = list(title = "Exposure")
)

colors.cv <- ChrAccR::getConfigElement("colorSchemesCont")[[".default.div"]]
c <- grDevices::colorRampPalette(colors.cv)(nrow(ann))

pdf(paste0(plot_dir, "chromVAR_motifs_080425.pdf"), width = 20, height = 20)
draw(Heatmap(zmat,
  cluster_rows = FALSE, # Preserve mTFR order
  cluster_columns = TRUE, 
  column_split = ann$Exposure, 
  top_annotation = column_ha, 
  col = c,
  heatmap_legend_param = list(title = "chromVAR"),
  show_column_names = FALSE
))
dev.off()


##################################################################################
# 3. Scatter Plots: mTFR vs chromVAR
##################################################################################
logger.start("Plotting mTFR vs chromVAR scatter plots")

# Reload Data for Heatmaps (Ensuring Clean State)
mtfr_devs <- list.files(mtfr_dir, pattern = "Monocyte", full.names = TRUE)
mtfr_devs <- list.cbind(lapply(condition, function(x) {
  methylTFR::deviations(readRDS(mtfr_devs[grepl(x, mtfr_devs)]))
}))
mtfr_devs <- as.data.frame(mtfr_devs)

# Filter Outliers and Compute Z-scores
mtfr_devs <- mtfr_devs[, !colnames(mtfr_devs) %in% outliers]
mtfr_devs_z <- methylTFR:::computeRowZScore(as.matrix(mtfr_devs))

chromvar_mat <- readRDS(paste0(ds_dir, "chromvar_jaspar2020_101224_corrected_zscores.rds"))
rownames(chromvar_mat) <- sub(".*_", "", rownames(chromvar_mat))

# Set max min to 5 and -5 for scatter
chromvar_mat <- pmin(pmax(chromvar_mat, -5), 5)

# --- Data Prep for Scatter ---
# Re-aligning data (using full datasets for scatter calculation)
# Note: Reuse mtfr_devs_z_sub (mono1 vs ctrl) or calculate fresh if "groups" variable changed
# Recalculating strict Mono1 vs Ctrl environment for scatter
mtfr_scatter_sub <- mtfr_devs_z[, groups == "C19_ctrl" | groups == "C19_sev_mono1"]
groups_scatter <- groups[groups == "C19_ctrl" | groups == "C19_sev_mono1"]

# Ensure ChromVAR matches
common_motifs <- intersect(rownames(mtfr_scatter_sub), rownames(chromvar_mat))
mtfr_sub <- mtfr_scatter_sub[common_motifs, ]
cvar_sub <- chromvar_mat[common_motifs, ]


# Setup Annotation for Calculation
cvar_groups <- ifelse(grepl("ctrl", colnames(cvar_sub)), "C19_ctrl", "C19_sev")
mtfr_groups <- groups_scatter

# --- Calculate Means ---
mtfr_ctrl_means <- rowMeans(mtfr_sub[, mtfr_groups == "C19_ctrl", drop = FALSE], na.rm = TRUE)
mtfr_sev_means  <- rowMeans(mtfr_sub[, grepl("C19_sev", mtfr_groups), drop = FALSE], na.rm = TRUE) 

cvar_ctrl_means <- rowMeans(cvar_sub[, cvar_groups == "C19_ctrl", drop = FALSE], na.rm = TRUE)
cvar_sev_means  <- rowMeans(cvar_sub[, cvar_groups == "C19_sev", drop = FALSE], na.rm = TRUE)

df_scatter <- data.frame(
  Motif = common_motifs,
  mTFR_diff = mtfr_ctrl_means - mtfr_sev_means,
  cVAR_diff = cvar_ctrl_means - cvar_sev_means
)

# --- Run Wilcoxon Tests ---
run_wilcox <- function(mat, group_vec) {
  apply(mat, 1, function(x) {
    ctrl_vals <- x[group_vec == "C19_ctrl"]
    sev_vals  <- x[grepl("C19_sev", group_vec)]
    
    if(length(na.omit(ctrl_vals)) >= 3 && length(na.omit(sev_vals)) >= 3) {
      res <- wilcox.test(ctrl_vals, sev_vals, exact = FALSE)
      return(res$p.value)
    } else {
      return(NA)
    }
  })
}

df_scatter$mTFR_pval <- run_wilcox(mtfr_sub, mtfr_groups)
df_scatter$cVAR_pval <- run_wilcox(cvar_sub, cvar_groups)

# --- Significance Status ---
p_cutoff <- 0.05
mtfr_is_sig <- !is.na(df_scatter$mTFR_pval) & df_scatter$mTFR_pval < p_cutoff
cvar_is_sig <- !is.na(df_scatter$cVAR_pval) & df_scatter$cVAR_pval < p_cutoff

df_scatter$Status <- "Not Differential"
df_scatter$Status[mtfr_is_sig & cvar_is_sig]   <- "Differential in both"
df_scatter$Status[!mtfr_is_sig & cvar_is_sig]  <- "Differential in cVAR"
df_scatter$Status[mtfr_is_sig & !cvar_is_sig]  <- "Differential in mTFR"

df_scatter$Status <- factor(df_scatter$Status, 
                            levels = c("Differential in both", "Differential in cVAR", 
                                       "Differential in mTFR", "Not Differential"))

plot_cor <- cor(df_scatter$mTFR_diff, df_scatter$cVAR_diff, use = "complete.obs", method = "pearson")

# --- Ranking Logic (P-value Based) ---
# Combined P-value for "Both" group (max of the pair)
df_scatter$max_pval <- pmax(df_scatter$mTFR_pval, df_scatter$cVAR_pval, na.rm = TRUE)

# Select Top 25 based on P-VALUES (Smallest is best)
top_both <- head(df_scatter[df_scatter$Status == "Differential in both", ][order(df_scatter[df_scatter$Status == "Differential in both", ]$max_pval), "Motif"], 25)
top_cvar <- head(df_scatter[df_scatter$Status == "Differential in cVAR", ][order(df_scatter[df_scatter$Status == "Differential in cVAR", ]$cVAR_pval), "Motif"], 25)
top_mtfr <- head(df_scatter[df_scatter$Status == "Differential in mTFR", ][order(df_scatter[df_scatter$Status == "Differential in mTFR", ]$mTFR_pval), "Motif"], 25)

candidate_motifs <- c(top_both, top_cvar, top_mtfr)

# --- Define User Target List ---
user_target_list <- c(
  "BHLHE23", "SPDEF", "SOX10", "EVX1", "MSGN1", "HOXC11", "HOXC12", "SOX12", "GLI3", "KLF3",
  "MEF2B", "HOXD8", "HOXA7", "HOXB8", "ZBTB7C", "TFAP2A", "RUNX2", "HOXD10", "HOXD12",
  "HOXD11", "HOXC9", "HOXB9", "HOXC10", "TFAP2A(var.3)", "ZNF263", "HOXD4", "HOXB4", "HOXB6",
  "ZFP57", "RELB", "PBX1", "REL", "RELA", "BACH2(var.2)", "CEBPD", "CEBPG(var.2)", "JUN",
  "FOS::JUN(var.2)", "ZNF528", "FOSL1::JUND(var.2)", "ETV1", "ELF1", "ELF5", "GABPA", "ELF3",
  "FOSL2::JUND(var.2)", "JUNB(var.2)", "CREM", "JUN::JUNB(var.2)", "JDP2(var.2)", "ATF7",
  "FOSL2::JUN(var.2)", "FOSL2::JUNB(var.2)", "FOSB::JUN", "FOSB::JUNB(var.2)", "FOSL1::JUN(var.2)",
  "JUND(var.2)", "ETV4", "EWSR1−FLI1", "ETV6", "E2F1", "E2F4", "GFI1", "CUX2", "CUX1",
  "NFIA", "TEAD1", "TEAD3", "TEAD2", "TEAD4", "MAZ", "E2F2", "PAX7", "ESX1", "ATF2",
  "HOXA9", "ZBTB32", "EBF3", "ZBED1", "IKZF1", "EHF", "NFKB1", "NFKB2", "POU1F1", "DMRTC2",
  "TP53", "DMRT3", "DMRTA2", "PLAGL2", "IRF1", "ZKSCAN5", "SPI1", "SPIB", "SPIC",
  "ZNF652", "GBX1", "OSR2", "NR1I3", "HNF4G", "HNF4A"
)

# Extract clean name
df_scatter$Clean_Name <- sub(".*_", "", df_scatter$Motif)

# --- Create Label Columns ---

# Label 1: Filtered (Intersection of Candidates AND User List)
df_scatter$Label_Filtered <- ifelse(
  (df_scatter$Motif %in% candidate_motifs) & (df_scatter$Clean_Name %in% user_target_list), 
  df_scatter$Clean_Name, 
  ""
)

# Label 2: Top Only (Candidates only, ignoring User List)
df_scatter$Label_Top <- ifelse(
  df_scatter$Motif %in% candidate_motifs, 
  df_scatter$Clean_Name, 
  ""
)

custom_colors <- c(
  "Differential in both" = "#A95A52",
  "Differential in cVAR" = "#649A6C",
  "Differential in mTFR" = "#4A588A",
  "Not Differential"     = "#DCE0E5" 
)

# --- Plot Generation ---

# Plot 1 = Diffs from heatmap + Filtered Labels (Intersection of Candidates AND User List)
p1 <- ggplot(df_scatter, aes(x = mTFR_diff, y = cVAR_diff)) +
  geom_point(aes(color = Status), alpha = 0.8, size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_color_manual(values = custom_colors) +
  geom_text_repel(aes(label = Label_Filtered), 
                  size = 3, 
                  max.overlaps = Inf, 
                  box.padding = 0.5,
                  point.padding = 0.5,
                  segment.color = NA, 
                  force = 2) +        
  annotate("text", x = max(df_scatter$mTFR_diff, na.rm=TRUE) * 0.75, 
           y = max(df_scatter$cVAR_diff, na.rm=TRUE) * 0.95, 
           label = sprintf("Correlation: %.2f", plot_cor), size = 5) +
  labs(
    title = "Monocytes (Targeted Labels)",
    x = "mTFR z-score difference (Control - Severe)",
    y = "chromVAR z-score difference (Control - Severe)",
    color = NULL 
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(hjust = 0.05, face = "bold", margin = margin(b = 15))
  )

pdf(paste0(plot_dir, "mTFR_vs_chromVAR_scatter_wilcox_pval_filtered.pdf"), width = 8, height = 8)
print(p1)
dev.off()

#  PLOT 2: Top Candidates (Unfiltered) - Labels based on Candidate List only, ignoring User List
p2 <- ggplot(df_scatter, aes(x = mTFR_diff, y = cVAR_diff)) +
  geom_point(aes(color = Status), alpha = 0.8, size = 2.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  scale_color_manual(values = custom_colors) +
  geom_text_repel(aes(label = Label_Top), 
                  size = 3, 
                  max.overlaps = Inf, 
                  box.padding = 0.5,
                  point.padding = 0.5,
                  segment.color = NA, 
                  force = 2) +        
  annotate("text", x = max(df_scatter$mTFR_diff, na.rm=TRUE) * 0.75, 
           y = max(df_scatter$cVAR_diff, na.rm=TRUE) * 0.95, 
           label = sprintf("Correlation: %.2f", plot_cor), size = 5) +
  labs(
    title = "Monocytes (Top Significant Labels)",
    x = "mTFR z-score difference (Control - Severe)",
    y = "chromVAR z-score difference (Control - Severe)",
    color = NULL 
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(hjust = 0.05, face = "bold", margin = margin(b = 15))
  )

pdf(paste0(plot_dir, "mTFR_vs_chromVAR_scatter_wilcox_pval_top_unfiltered.pdf"), width = 8, height = 8)
print(p2)
dev.off()

logger.completed()

############################################################
# SELEX based coloring
###################################################################

selex <- read.csv("/icbb/projects/igunduz/exposure_atlas_manuscript/sample_annots/Selex_data.csv", skip = 20, header = 21, sep = ";")[, c(1, 2, 3, 4, 6)]
motifs <- rownames(mtfr_devs_z) # Use the same motifs as in the heatmap (mono1 vs ctrl)
plot_cor <- cor(df_scatter$mTFR_diff, df_scatter$cVAR_diff, use = "complete.obs", method = "spearman") 

# Add a new annotation column based on selex$Call
# Ensure the order of `motifs` matches the heatmap rows
selex_annotation <- data.frame(
  Call = ifelse(motifs %in% selex$TF.name,
    selex$Call[match(motifs, selex$TF.name)],
    "Inconclusive"
  ) # Replace NA with "Inconclusive"
)
rownames(selex_annotation) <- motifs
selex_annotation$Call[is.na(selex_annotation$Call)] <- "Inconclusive"
selex_annotation$TF <- rownames(selex_annotation)

clean_selex <- function(x) {
  # 1. Fix known typos
  x <- gsub("Mehyl", "Methyl", x)
  
  # 2. Convert to lowercase for easier matching
  x_lower <- tolower(x)
  
  # 3. Default to Inconclusive
  res <- rep("Inconclusive", length(x))
  
  # 4. Assign categories (Order matters for priority)
  
  # 'Little effect' -> Little effect
  res[grep("little effect", x_lower)] <- "Little effect"
  
  # 'MethylMinus' -> MethylMinus (overwrites "Little effect and MethylMinus")
  res[grep("methylminus", x_lower)] <- "MethylMinus"
  
  # 'MethylPlus' -> MethylPlus (overwrites "Little effect and MethylPlus")
  res[grep("methylplus", x_lower)] <- "MethylPlus"
  
  # 'Mixed' -> If both Plus and Minus are present
  mixed_idx <- grep("methylplus", x_lower)
  mixed_idx <- intersect(mixed_idx, grep("methylminus", x_lower))
  res[mixed_idx] <- "Mixed"
  
  # 'Inconclusive' explicitly
  res[grep("inconclusive", x_lower)] <- "Inconclusive"
  
  return(res)
}

df_scatter$Selex_Call <- selex_annotation$Call[match(df_scatter$Motif, selex_annotation$TF)]
df_scatter$Selex_Call_Clean <- clean_selex(df_scatter$Selex_Call)
selex_colors <- c(
  "MethylPlus"    = "#D55E00",  # Red/Orange
  "MethylMinus"   = "#0072B2",  # Blue
  "Mixed"         = "#CC79A7",  # Purple
  "Little effect" = "#009E73",  # Green
  "Inconclusive"  = "grey85"    # Light Grey
)

p3 <- ggplot(df_scatter, aes(x = mTFR_diff, y = cVAR_diff)) +
  # Points colored by cleaned SELEX call
  geom_point(aes(color = Selex_Call_Clean), alpha = 0.8, size = 2.5) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
  
  scale_color_manual(values = selex_colors) +
  
  # Label top candidates (Unfiltered)
  geom_text_repel(aes(label = Label_Top), 
                  size = 3, 
                  max.overlaps = Inf, 
                  box.padding = 0.5,
                  point.padding = 0.5,
                  segment.color = NA, 
                  force = 2) +        
  
  annotate("text", x = max(df_scatter$mTFR_diff, na.rm=TRUE) * 0.75, 
           y = max(df_scatter$cVAR_diff, na.rm=TRUE) * 0.95, 
           label = sprintf("Correlation: %.2f", plot_cor), size = 5) +
  
  labs(
    title = "Monocytes (SELEX Coloring)",
    x = "mTFR z-score difference (Control - Severe)",
    y = "chromVAR z-score difference (Control - Severe)",
    color = "SELEX Call" 
  ) +
  
  theme_classic(base_size = 16) +
  theme(
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", color = NA),
    legend.key = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(hjust = 0.05, face = "bold", margin = margin(b = 15))
  )

pdf(paste0(plot_dir, "mTFR_vs_chromVAR_scatter_wilcox_selex.pdf"), width = 10, height = 10)
print(p3)
dev.off()
