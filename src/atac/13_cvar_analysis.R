#!/usr/bin/env Rscript

#####################################################################
# 13_cvar_analysis.R
# Created by IBG on 03-11-2025
# Corrected and Updated
#####################################################################

suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(ggplot2)
  library(chromVAR)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(stringr)
  library(grid)
  library(RColorBrewer)
})
set.seed(12)

# -------------------------------------------------------------------
# 1. Setup and Load Project
# -------------------------------------------------------------------
outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/"
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)
save_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/sample_annots/cvar_diffs"
if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

# Define Colors and Valid Cell Types 
cell_type_colors <- c(
  "B_mem" = "#AE017E",
  "B_naive" = "#F768A1",
  "DC" = "#67000D",
  "Mono_CD14" = "#FE9929",
  "Mono_CD16" = "#CC4C02",
  "NK_CD16" = "#A65628",
  "Plasma" = "#A106BD",
  "T_mait" = "#41B6C4",
  "T_mem_CD4" = "#4292c6",
  "T_mem_CD8" = "#0074cc",
  "T_mix" = "#888FB5",
  "T_naive" = "#C7E9B4"
)
names(cell_type_colors) <- names(cell_type_colors)

# -------------------------------------------------------------------
# 2. Create Pseudo-bulk and Deviations Object
# -------------------------------------------------------------------
# Re-create pseudo-bulk samples
exposure_cell <- paste0(project$sample_exposure_group, "_", project$ClusterCellTypes, "_", project$Sample)

project <- addCellColData(
  ArchRProj = project, data = exposure_cell,
  name = "exposure_cell", cells = project$cellNames,
  force = TRUE
)

# Get the chromVAR matrix
pseudo_chrom <- ArchR::getGroupSE(project, "altiusMatrix", groupBy = "exposure_cell", divideN = TRUE)

# Select z-scores and create chromVAR object
seZ <- pseudo_chrom[rowData(pseudo_chrom)$seqnames == "z", ]
zmat <- assay(seZ)
rownames(zmat) <- rowData(seZ)$name
devmat <- zmat

se_dev <- SummarizedExperiment(
  assays = SimpleList(deviations = devmat, z = zmat),
  rowData = DataFrame(
    name = rowData(seZ)$name,
    fractionMatches = rep(NA, nrow(seZ)),
    fractionBackgroundOverlap = rep(NA, nrow(seZ))
  ),
  colData = colData(seZ)
)
dev <- as(se_dev, "chromVARDeviations")

# -------------------------------------------------------------------
# 3. ROBUST METADATA PARSING (Fixing the NA issue)
# -------------------------------------------------------------------
cn <- colnames(dev)

# A. Extract Condition (everything before the first underscore usually, but handled specifically)
colData(dev)$Condition <- str_extract(cn, "^[^_]+_[^_]+")
colData(dev)$Exposure <- str_extract(cn, "^[^_]+")

# B. Clean the remaining string to extract CellType and Sample
# Remove the Condition prefix
remainder <- str_remove(cn, paste0(colData(dev)$Condition, "_"))
# Remove file extension
remainder <- str_remove(remainder, "\\.tsv\\.gz$")
# Remove the _Bei tag if it exists in the middle
remainder <- str_replace(remainder, "_Bei_", "_")

# C. Match Cell Types Deterministically
# Instead of regex splitting, we check which valid cell type appears at the start of the remainder
colData(dev)$Cell_Type <- NA_character_
colData(dev)$Sample <- NA_character_

for (ct in names(cell_type_colors)) {
  # Pattern: Start of string, followed by cell type, followed by underscore
  # This distinguishes "T_mem" from "T_mem_CD4" by ensuring the longest match logic or exact boundaries
  pattern <- paste0("^", ct, "_")
  
  # Find indices where this cell type matches
  matches <- str_detect(remainder, pattern)
  
  # Assign Cell Type
  colData(dev)$Cell_Type[matches] <- ct
  
  # Extract Sample (remove the cell type + underscore from the remainder)
  colData(dev)$Sample[matches] <- str_remove(remainder[matches], pattern)
}

# Verify no NAs remain
if (sum(is.na(colData(dev)$Cell_Type)) > 0) {
  warning("Some Cell Types resulted in NA. Check 'names(cell_type_colors)' list against these names:")
  print(head(remainder[is.na(colData(dev)$Cell_Type)]))
}

# Save the cleaned object
saveRDS(dev, file = file.path(save_dir, "/atac_chromvar_deviations_exposure_celltype.rds"))
dev <- readRDS(file.path(save_dir, "/atac_chromvar_deviations_exposure_celltype.rds"))

# -------------------------------------------------------------------
# 4. Differential Analysis (Per Comparison AND Per Cell Type)
# -------------------------------------------------------------------

raw_comps <- c(
  "HIV_ctrl vs HIV_chr",
  "HIV_ctrl vs HIV_acu",
  "Influenza_d3 vs Influenza_ctrl",
  "Influenza_d30 vs Influenza_ctrl",
  "Influenza_d6 vs Influenza_ctrl",
  "OP_high vs OP_low",
  "OP_high vs OP_med",
  "OP_med vs OP_low",
  "C19_mild vs C19_ctrl",
  "C19_mod vs C19_ctrl",
  "C19_sev vs C19_ctrl"
)

results_list <- list()

# Outer Loop: Comparisons
for(comp_str in raw_comps) {
  
  parts <- str_split(comp_str, " vs ")[[1]]
  group1 <- parts[1]
  group2 <- parts[2]
  
  print(paste0("Processing comparison: ", comp_str))
  
  # Inner Loop: Cell Types
  for(ct in names(cell_type_colors)) {
    
    # Subset: Specific Exposure Groups AND Specific Cell Type
    subset_mask <- (colData(dev)$Condition %in% c(group1, group2)) & 
                   (colData(dev)$Cell_Type == ct)
    
    dev_sub <- dev[, subset_mask]
    
    # Check if we have data for both groups in this cell type
    # We need at least 2 samples per group generally, but at least 1 to run without immediate error
    # (differentialDeviations usually requires replicates, check your sample size)
    n_g1 <- sum(colData(dev_sub)$Condition == group1)
    n_g2 <- sum(colData(dev_sub)$Condition == group2)
    
    if(n_g1 > 1 && n_g2 > 1) {
      
      # Run Differential Deviations
      tryCatch({
        diff_res <- differentialDeviations(dev_sub, groups = "Condition")
        
        # Add Metadata
        diff_res$Motif_ID <- rownames(diff_res)
        diff_res$Comparison <- paste0(group1, "_vs_", group2)
        diff_res$Group1 <- group1
        diff_res$Group2 <- group2
        diff_res$Cell_Type <- ct # Crucial: Label the result with the cell type
        
        # Store in list
        list_name <- paste0(group1, "_vs_", group2, "_", ct)
        results_list[[list_name]] <- diff_res
        
      }, error = function(e) {
        message(paste("Skipping", comp_str, "for", ct, "- Error or insufficient variance."))
      })
      
    } else {
      message(paste("Skipping", comp_str, "for", ct, "- Insufficient samples (", n_g1, "vs", n_g2, ")"))
    }
  }
}

# Merge and Save
if(length(results_list) > 0) {
  final_merged_results <- do.call(rbind, results_list)
  saveRDS(final_merged_results, file.path(save_dir, "all_exposures_celltype_differential_deviations_merged.rds"))
} else {
  warning("No differential results were generated.")
}

# -------------------------------------------------------------------
# 5. Global PCA Plot
# -------------------------------------------------------------------
z_scores <- deviationScores(dev)
z_scores <- z_scores[complete.cases(z_scores), ]

# Run PCA
pca_res <- prcomp(t(z_scores), center = TRUE, scale. = TRUE)

# Prepare Data
pca_df <- as.data.frame(pca_res$x[, 1:2])
# Ensure alignment
pca_df$Cell_Type <- colData(dev)$Cell_Type[match(rownames(pca_df), colnames(dev))]
pca_df$Exposure <- colData(dev)$Exposure[match(rownames(pca_df), colnames(dev))]

# Plot
pca_plot <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cell_Type, shape = Exposure)) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_color_manual(values = cell_type_colors) +
  scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 4, 10)[1:length(unique(pca_df$Exposure))]) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "right"
  ) +
  labs(title = "Global ChromVAR Deviations PCA", x = "PC1", y = "PC2")

ggsave(filename = file.path(save_dir, "global_deviations_pca.pdf"), plot = pca_plot, width = 10, height = 8)
print("Analysis Complete.")

# -------------------------------------------------------------------
# 6. Heatmap of Top Differential Motifs
# -------------------------------------------------------------------
# Filter for significance and high deviation difference
if (exists("final_merged_results")) {
  top_hits <- final_merged_results %>%
    as.data.frame() %>%
    filter(p_value_adjusted < 0.05) %>%
    pull(Motif_ID) %>%
    unique()
  
  if (length(top_hits) < 5) {
    stop("Fewer than 5 significant motifs found. Adjust p-value threshold or check data.")
  }
} else {
  stop("Result object 'final_merged_results' not found.")
}

condition_cell <- paste0(project$sample_exposure_group, "_", project$ClusterCellTypes)

project <- addCellColData(
  ArchRProj = project, data = condition_cell,
  name = "condition_cell", cells = project$cellNames,
  force = TRUE
)

# Get the chromVAR matrix
pseudo_chrom <- ArchR::getGroupSE(project, "altiusMatrix", groupBy = "condition_cell", divideN = TRUE)

# Select z-scores and create chromVAR object
seZ <- pseudo_chrom[rowData(pseudo_chrom)$seqnames == "z", ]
zmat <- assay(seZ)
rownames(zmat) <- rowData(seZ)$name
devmat <- zmat[top_hits, ]


# -------------------------------------------------------------------
# 6. Heatmap with Figure 1A Color Matching
# -------------------------------------------------------------------

condition_colors <- c(
  # COVID-19 (Greens)
  "C19_ctrl"      = "#2b8cbe",  # Control (Blue-ish based on Figure 1A legend)
  "C19_mild"      = "#a1d99b",  # Light Green
  "C19_mod"       = "#41ab5d",  # Medium Green
  "C19_sev"       = "#006d2c",  # Dark Green

  # Influenza (Oranges)
  "Influenza_ctrl"= "#2b8cbe",  # Control (Matching COVID control if they are the same baseline)
  "Influenza_d3"  = "#fec44f",  # Day 3 (Yellow/Orange)
  "Influenza_d6"  = "#fe9929",  # Day 6 (Orange)
  "Influenza_d30" = "#d95f0e",  # Day 28/30 (Dark Orange/Red)

  # HIV (Purples)
  "HIV_ctrl"      = "#4F619D",  # Control (Lavender/Pink)
  "HIV_chr"       = "#825CA3",  # Chronic (Medium Purple)
  "HIV_acu"       = "#893368",  # Acute (Dark Magenta)

  # Organophosphate (Blues)
  "OP_low"        = "#9ecae1",  # Low (Light Blue)
  "OP_med"        = "#4292c6",  # Medium (Medium Blue)
  "OP_high"       = "#08519c"   # High (Dark Blue)
)

# Parse Metadata
cn <- colnames(devmat) 
anno_df <- data.frame(ID = cn, stringsAsFactors = FALSE)
sorted_conds <- names(condition_colors)[order(nchar(names(condition_colors)), decreasing = TRUE)]
cond_pattern <- paste0("^(", paste(sorted_conds, collapse = "|"), ")")
anno_df$Condition <- str_extract(cn, cond_pattern)
sorted_types <- names(cell_type_colors)[order(nchar(names(cell_type_colors)), decreasing = TRUE)]
ct_pattern <- paste0("_(", paste(sorted_types, collapse = "|"), ")(_|$)" )
anno_df$Cell_Type <- str_extract(cn, ct_pattern)
anno_df$Cell_Type <- str_remove_all(anno_df$Cell_Type, "^_|_$")


# Heatmap Z-score Colors (Red-White-Blue)
colors.cv <- tryCatch({
  ChrAccR::getConfigElement("colorSchemesCont")[[".default.div"]]
}, error = function(e) {
  c("#2166AC", "white", "#B2182B") 
})
col_fun <- colorRamp2(seq(-5, 5, length.out = length(colors.cv)), colors.cv)

# 3. Create Heatmap
col_ha <- HeatmapAnnotation(
  df = anno_df[, c("Cell_Type", "Condition")],
  col = list(
    Cell_Type = cell_type_colors,
    Condition = condition_colors
  ),
  annotation_name_side = "left",
  simple_anno_size = unit(0.4, "cm")
)

pdf(file.path(save_dir, "differential_motifs_heatmap_fig1_colors.pdf"), width = 14, height = 10)

Heatmap(
  matrix = devmat,
  name = "z-score",
  col = col_fun, 
  top_annotation = col_ha,
  column_split = anno_df$Cell_Type,
  cluster_column_slices = FALSE, # Keep Cell Types in order
  column_title_rot = 45,
  column_gap = unit(1, "mm"),
  show_column_names = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  column_title = "Differential Deviations",
  use_raster = TRUE
)

dev.off()

# -------------------------------------------------------------------
# Individual cell-type heatmaps
# -------------------------------------------------------------------
dend_global <- hclust(dist(devmat), method = "ward.D2")
global_row_order <- dend_global$order

# --- 3. Setup Colors ---
colors.cv <- tryCatch({
  ChrAccR::getConfigElement("colorSchemesCont")[[".default.div"]]
}, error = function(e) {
  c("#2166AC", "white", "#B2182B") 
})
col_fun <- colorRamp2(seq(-5, 5, length.out = length(colors.cv)), colors.cv)

# --- 4. Define Grid Layout Parameters ---
unique_types <- sort(unique(anno_df$Cell_Type))
n_plots <- length(unique_types)

# Define grid dimensions (e.g., 3 columns wide)
ncol_grid <- 3 
nrow_grid <- ceiling(n_plots / ncol_grid)

# Create PDF with enough size for the grid
pdf(file.path(save_dir, "differential_motifs_heatmap_grid_layout.pdf"), width = 20, height = 5 * nrow_grid)

# Create the Grid Page
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow_grid, ncol_grid)))

# --- 5. Loop and Draw ---
for (i in seq_along(unique_types)) {
  
  ct <- unique_types[i]
  
  # A. Subset Data
  idx <- which(anno_df$Cell_Type == ct)
  sub_mat <- devmat[, idx, drop = FALSE]
  sub_anno_df <- anno_df[idx, , drop = FALSE]
  
  # B. Create Annotation
  col_ha_sub <- HeatmapAnnotation(
    df = sub_anno_df[, "Condition", drop = FALSE],
    col = list(Condition = condition_colors),
    show_legend = FALSE, # Hide annotation legend to save space
    simple_anno_size = unit(0.3, "cm"),
    show_annotation_name = FALSE
  )
  
  # C. Create Heatmap Object
  ht <- Heatmap(
    matrix = sub_mat,
    name = "z-score",
    col = col_fun,
    top_annotation = col_ha_sub,
    
    # 1. Clustering Columns WITHIN groups
    # Split by condition, KEEP condition order (e.g. Ctrl then Disease), CLUSTER samples inside
    column_split = sub_anno_df$Condition, 
    cluster_column_slices = FALSE, 
    cluster_columns = TRUE,        
    
    # 2. Consistent Row Ordering
    cluster_rows = FALSE,          # Disable internal clustering
    row_order = global_row_order,  # Enforce the global order
    
    # 3. Aesthetics
    column_title = ct,
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    show_column_names = FALSE,
    
    show_row_names = TRUE, # We need row names on all since they are in a grid
    row_names_gp = gpar(fontsize = 6),
    
    show_heatmap_legend = (i == 1),
    use_raster = TRUE
  )
  
  row_idx <- ceiling(i / ncol_grid)
  col_idx <- (i - 1) %% ncol_grid + 1
  
  pushViewport(viewport(layout.pos.row = row_idx, layout.pos.col = col_idx))
  draw(ht, newpage = FALSE) 
  popViewport()
}

lgd_cond = Legend(title = "Condition", at = names(condition_colors), legend_gp = gpar(fill = condition_colors))
draw(lgd_cond, x = unit(0.95, "npc"), y = unit(0.05, "npc"), just = c("right", "bottom"))

dev.off()
