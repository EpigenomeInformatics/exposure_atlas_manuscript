library(ComplexHeatmap)
library(circlize)
library(chromVAR)
library(dplyr)
library(RColorBrewer)

# --- 1. Data Setup (Assumes 'dev' and 'final_merged_results' are loaded) ---

# Select Significant Motifs (Adjust p-value threshold if too many/few rows)
sig_motifs <- final_merged_results %>%
  filter(p_value_adjusted < 0.05) %>%
  pull(Motif_ID) %>%
  unique()

if(length(sig_motifs) < 2) stop("Not enough significant motifs to plot.")

# Get Z-scores
z_mat <- deviationScores(dev)
z_mat_sig <- z_mat[sig_motifs, , drop=FALSE]
z_mat_sig <- z_mat_sig[complete.cases(z_mat_sig), ] # Remove NAs

# --- 2. Metadata & Colors ---

# Align metadata
col_meta <- colData(dev) %>% as.data.frame()
col_meta <- col_meta[colnames(z_mat_sig), ] 

# Colors (Your custom palette)
cell_type_colors <- c(
  "B_mem" = "#AE017E", "B_naive" = "#F768A1", "DC" = "#67000D",
  "Mono_CD14" = "#FE9929", "Mono_CD16" = "#CC4C02", "NK_CD16" = "#A65628",
  "Plasma" = "#A106BD", "T_mait" = "#41B6C4", "T_mem_CD4" = "#4292c6",
  "T_mem_CD8" = "#0074cc", "T_mix" = "#888FB5", "T_naive" = "#C7E9B4"
)

# Exposure Colors
unique_exposures <- unique(col_meta$Exposure)
# Using Dark2 palette, but you can assign specific colors like "HIV"="red", "Control"="grey"
exp_colors <- setNames(brewer.pal(n = max(3, length(unique_exposures)), name = "Dark2")[1:length(unique_exposures)], unique_exposures)

# Annotation Bar
ha <- HeatmapAnnotation(
  Cell_Type = col_meta$Cell_Type,
  Exposure = col_meta$Exposure,
  col = list(
    Cell_Type = cell_type_colors,
    Exposure = exp_colors
  ),
  annotation_name_side = "left"
)

# --- 3. Define the Split Logic ---
# We create a data frame for splitting. 
# The order of columns determines the hierarchy: Cell_Type first, then Exposure.
split_df <- data.frame(
  Cell_Type = col_meta$Cell_Type,
  Exposure = col_meta$Exposure
)

# Optional: Enforce specific order for Exposures (e.g., put Control first)
# split_df$Exposure <- factor(split_df$Exposure, levels = c("Control", "HIV", "Influenza", "OP", "C19")) 

# --- 4. Generate Heatmap ---

# Cap scores for visualization
z_capped <- t(scale(t(z_mat_sig)))
z_capped[z_capped > 3] <- 3
z_capped[z_capped < -3] <- -3

col_fun <- colorRamp2(c(-3, 0, 3), c("#2166ac", "white", "#b2182b"))

pdf(file.path(save_dir, "Figure_Global_Heatmap_Clustered_by_Exposure.pdf"), width = 20, height = 14)

Heatmap(z_capped,
        name = "Z-Score",
        col = col_fun,
        top_annotation = ha,
        
        # --- THE KEY CHANGE IS HERE ---
        # Passing the dataframe splits hierarchically (CellType -> Exposure)
        column_split = split_df,
        
        # Keep slices fixed (alphabetical or factor order) so "Control" is always left of "HIV"
        cluster_column_slices = FALSE, 
        
        # Cluster samples strictly WITHIN the Exposure groups
        cluster_columns = TRUE,
        
        # Aesthetics
        column_title_gp = gpar(fontsize = 8), # Make split titles smaller to fit
        column_title_rot = 45, # Rotate titles if they overlap
        show_column_names = FALSE,
        show_row_names = length(sig_motifs) < 100, 
        row_names_gp = gpar(fontsize = 6),
        
        # Row Clustering (Find gene modules)
        row_km = 5, 
        border = TRUE
)

dev.off()
print("Cluster-by-Exposure Heatmap saved.")