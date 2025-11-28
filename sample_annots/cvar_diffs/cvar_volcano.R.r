library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)

# 1. Load your merged results (ensure this object exists from the previous calculation step)
# res <- readRDS(file.path(save_dir, "all_exposures_celltype_differential_deviations_merged.rds"))

# 2. Create specific directory for these plots
volcano_dir <- file.path(save_dir, "volcano_plots_individual")
if (!dir.exists(volcano_dir)) dir.create(volcano_dir)

# 3. Get lists of Comparisons and Cell Types actually present in the data
comparisons_found <- unique(final_merged_results$Comparison)
cell_types_found <- unique(final_merged_results$Cell_Type)

# 4. Nested Loop to Plot and Save Individually
for (comp in comparisons_found) {
  
  print(paste0("Plotting comparison: ", comp))
  
  for (ct in cell_types_found) {
    
    # Filter data for this specific combination
    plot_df <- final_merged_results %>% 
      filter(Comparison == comp, Cell_Type == ct)
    
    # Skip if no data found for this combo (e.g., if a cell type isn't in that batch)
    if(nrow(plot_df) == 0) next
    
    # --- PREPARE DATA ---
    padj_cutoff <- 0.05
    
    # Create Significance Labels
    plot_df <- plot_df %>%
      mutate(
        Significance = case_when(
          p_value_adjusted < padj_cutoff & meanDiff > 0 ~ "Enriched in Group 1",
          p_value_adjusted < padj_cutoff & meanDiff < 0 ~ "Enriched in Group 2",
          TRUE ~ "NS"
        )
      )
    
    # Pick top labels (Top 10 significant by P-value)
    top_labels <- plot_df %>%
      filter(Significance != "NS") %>%
      arrange(p_value_adjusted) %>%
      slice_head(n = 10)
    
    # Dynamic Titles
    g1_name <- unique(plot_df$Group1)
    g2_name <- unique(plot_df$Group2)
    
    # --- GENERATE PLOT ---
    p <- ggplot(plot_df, aes(x = meanDiff, y = -log10(p_value_adjusted))) +
      # Base points (Not Significant)
      geom_point(data = subset(plot_df, Significance == "NS"), 
                 color = "grey85", size = 1, alpha = 0.5) +
      
      # Significant points
      geom_point(data = subset(plot_df, Significance != "NS"), 
                 aes(color = Significance), size = 2, alpha = 0.8) +
      
      # Colors: Red (Group 1), Blue (Group 2)
      scale_color_manual(values = c("Enriched in Group 1" = "#B2182B", 
                                    "Enriched in Group 2" = "#2166AC")) +
      
      # Text Labels
      geom_text_repel(data = top_labels, aes(label = Motif_ID),
                      size = 3, max.overlaps = 20, box.padding = 0.5,
                      segment.color = "grey50") +
      
      # Aesthetics
      theme_bw() +
      theme(
        legend.position = "top",
        panel.grid.minor = element_blank()
      ) +
      
      # Lines and Limits
      geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
      geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey50") +
      
      # Labels
      labs(
        title = paste0(g1_name, " vs ", g2_name),
        subtitle = paste0("Cell Type: ", ct),
        x = paste0("Mean Deviation Diff\n<-- Higher in ", g2_name, " | Higher in ", g1_name, " -->"),
        y = "-log10(Adjusted P-value)",
        color = ""
      )

    # --- SAVE PLOT ---
    # We use ggsave which automatically closes the device after saving
    # File name: Volcano_HIV_ctrl_vs_HIV_chr_B_mem.pdf
    safe_comp_name <- str_replace_all(comp, " ", "_") # Remove spaces if any
    filename <- paste0("Volcano_", safe_comp_name, "_", ct, ".pdf")
    
    tryCatch({
      ggsave(filename = file.path(volcano_dir, filename), 
             plot = p, width = 6, height = 6, dpi = 300)
    }, error = function(e) {
      message(paste("Failed to save:", filename))
    })
  }
}

print("All volcano plots generated.")