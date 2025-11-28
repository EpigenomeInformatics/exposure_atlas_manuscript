library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)

# 1. Get Data
z_mat <- deviationScores(dev)
sig_motifs <- final_merged_results %>% filter(p_value_adjusted < 0.05) %>% pull(Motif_ID) %>% unique()
z_sig <- z_mat[sig_motifs, ]

# 2. Calculate Average Z-Score per Group (Cell_Type + Exposure)
# Transpose to dataframe
df_z <- as.data.frame(t(z_sig))
df_z$Group <- paste0(colData(dev)$Cell_Type, "::", colData(dev)$Exposure)

# Aggregate (Mean)
avg_df <- df_z %>%
  group_by(Group) %>%
  summarise(across(everything(), mean, na.rm=TRUE)) %>%
  column_to_rownames("Group")

# 3. Prepare Matrix for Plotting
avg_mat <- t(as.matrix(avg_df))
avg_mat <- t(scale(t(avg_mat))) # Row-scale to emphasize relative differences

# 4. Split Names for Annotation
group_names <- colnames(avg_mat)
cell_types <- str_split_fixed(group_names, "::", 2)[,1]
exposures  <- str_split_fixed(group_names, "::", 2)[,2]

# 5. Plot
ha <- HeatmapAnnotation(
  Exposure = exposures,
  CellType = cell_types,
  col = list(CellType = cell_type_colors, Exposure = exp_colors)
)

pdf(file.path(save_dir, "Figure_Summary_Average_Heatmap.pdf"), width = 12, height = 10)
Heatmap(avg_mat,
        name = "Avg Z-Score",
        top_annotation = ha,
        column_split = cell_types, # Keep cell types separate
        cluster_columns = TRUE,    # Cluster exposures to see which are similar!
        show_column_names = FALSE,
        row_km = 4,                # Cluster TFs
        col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
        border = TRUE)
dev.off()

library(umap)
library(ggplot2)

# 1. Use the Deviation Matrix (Motifs x Cells/Samples)
# We transpose it: Rows = Samples, Columns = Motifs
# Use only significant motifs to reduce noise
umap_input <- t(z_sig) 

# 2. Run UMAP
set.seed(123)
umap_res <- umap(umap_input)
umap_df <- as.data.frame(umap_res$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")

# 3. Add Metadata
umap_df$Cell_Type <- colData(dev)$Cell_Type
umap_df$Exposure <- colData(dev)$Exposure

# 4. Plot
# Facet by Cell Type to see if Exposures shift the state
p <- ggplot(umap_df, aes(x=UMAP1, y=UMAP2, color=Exposure)) +
  geom_point(size=3, alpha=0.8) +
  scale_color_manual(values=exp_colors) +
  theme_bw() +
  facet_wrap(~Cell_Type, scales = "free") + # Crucial: Look at each lineage independently
  labs(title = "Regulatory Landscape (Motif Deviations)",
       subtitle = "Do exposures drive distinct regulatory states?")

ggsave(file.path(save_dir, "Figure_Regulatory_UMAP.pdf"), p, width = 14, height = 10)