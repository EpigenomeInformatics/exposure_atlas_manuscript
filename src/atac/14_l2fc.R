suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(pheatmap)  
  library(RColorBrewer)
})

# --- Configuration ---
mainDir  <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_2023-10-02/"
covidDir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/"
plotDir  <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/plots/correlations_matrix/"

if (!dir.exists(plotDir)) dir.create(plotDir, recursive = TRUE)

sensitivity_limit <- 0.5 

cells <- c("B_mem", "B_naive", "Mono_CD14", "Mono_CD16", "NK_CD16", 
           "T_mem_CD8", "T_mem_CD4", "T_naive", "T_mait")

main_comps <- c(
"HIV_ctrl vs HIV_chr",
"HIV_ctrl vs HIV_acu",
"Influenza_d3 vs Influenza_ctrl",
"Influenza_d6 vs Influenza_ctrl" ,
"Influenza_d30 vs Influenza_ctrl",
"OP_high vs OP_low",
"OP_med vs OP_low"

)

covid_comps <- c(
  "C19_mild vs C19_ctrl",
  "C19_mod vs C19_ctrl",
  "C19_sev vs C19_ctrl"
)

# --- Colors ---
condition_colors <- c(
  "C19_ctrl" = "#2b8cbe", "C19_mild" = "#a1d99b", "C19_mod" = "#41ab5d", "C19_sev" = "#006d2c",
  "Flu_ctrl" = "#2b8cbe", "Flu_d3" = "#fec44f", "Flu_d6" = "#fe9929", "Flu_d30" = "#d95f0e",
  "HIV_ctrl" = "#4F619D", "HIV_chr" = "#825CA3", "HIV_acu" = "#893368",
  "OP_low" = "#9ecae1", "OP_med" = "#4292c6", "OP_high" = "#08519c"
)

comp_group_colors <- c(
  "HIV" = "#88419D", 
  "Flu" = "#D95F0E", 
  "C19" = "#238B45", 
  "OP"  = "#084594", 
  "Other" = "grey50"
)

cell_colors <- c(
  "B_mem" = "#AE017E", "B_naive" = "#F768A1", "DC" = "#67000D",
  "Mono_CD14" = "#FE9929", "Mono_CD16" = "#CC4C02", "NK_CD16" = "#A65628",
  "Plasma" = "#A106BD", "T_mait" = "#41B6C4", "T_mem_CD4" = "#4292c6",
  "T_mem_CD8" = "#0074cc", "T_mix" = "#888FB5", "T_naive" = "#C7E9B4"
)

# --- Helper Functions ---
make_label <- function(cell, comp_name) {
  n <- str_remove(comp_name, " \\[.*\\]")
  n <- str_replace_all(n, "Influenza", "Flu")
  n <- str_replace_all(n, "C19_", "C19_") 
  n <- str_replace(n, " vs ", "_v_")
  return(paste0(cell, ":", n))
}

load_lfc <- function(base_dir, cell, target_comps) {
  res <- list()
  compTablePath <- file.path(base_dir, cell, "reports/differential_data/comparisonTable.rds")
  if (!file.exists(compTablePath)) return(list())
  
  compTable <- readRDS(compTablePath)
  compTable$plain_name <- paste0(compTable$grp1, " vs ", compTable$grp2)
  matched_idx <- which(compTable$plain_name %in% target_comps)
  
  for (i in matched_idx) {
    f <- file.path(base_dir, cell, "reports/differential_data", paste0("diffTab_", i, "_archrPeaks.tsv"))
    if (file.exists(f)) {
      dt <- fread(f, select = c("chrom", "chromStart", "chromEnd", "log2FoldChange"))
      dt[, peakID := paste(chrom, chromStart, chromEnd, sep = "_")]
      global_name <- make_label(cell, compTable$plain_name[i])
      tmp <- dt[, .(peakID, log2FoldChange)]
      setnames(tmp, "log2FoldChange", global_name)
      res[[global_name]] <- tmp
    }
  }
  return(res)
}

# --- Main Execution ---
message("Starting Global Data Load...")
global_lfc_list <- list()

for (cell in cells) {
  global_lfc_list <- c(global_lfc_list, load_lfc(mainDir, cell, main_comps))
  global_lfc_list <- c(global_lfc_list, load_lfc(covidDir, cell, covid_comps))
}

num_loaded <- length(global_lfc_list)
message(paste0("Total Loaded: ", num_loaded))

if (num_loaded > 1) {
  # Merge & Filter
  merged_df <- Reduce(function(x, y) merge(x, y, by = "peakID", all = TRUE), global_lfc_list)
  mat_data <- as.matrix(merged_df[, -1, with = FALSE])
  rownames(mat_data) <- merged_df$peakID
  mat_data[is.na(mat_data)] <- 0
  
  row_vars <- apply(mat_data, 1, var)
  keep_idx <- order(row_vars, decreasing = TRUE)[1:min(25000, length(row_vars))]
  cor_mat <- cor(mat_data[keep_idx, ], method = "pearson")

  # Annotations
  all_labels <- colnames(cor_mat)
  annotation_df <- data.frame(row.names = all_labels)
  annotation_df$Cell_Type <- sapply(all_labels, function(x) str_split(x, ":")[[1]][1])
  
  annotation_df$Group <- sapply(all_labels, function(x) {
    if (grepl("HIV", x)) return("HIV")
    if (grepl("Flu", x)) return("Flu")
    if (grepl("C19", x)) return("C19")
    if (grepl("OP", x)) return("OP")
    return("Other")
  })
  
  annotation_df$Condition <- sapply(all_labels, function(lbl) {
    comp_part <- str_split(lbl, ":")[[1]][2]
    targets <- c("C19_sev", "C19_mod", "C19_mild", "Flu_d30", "Flu_d6", "Flu_d3", 
                 "HIV_acu", "HIV_chr", "OP_high", "OP_med", "OP_low")
    found <- NA
    for (t in targets) { if (grepl(t, comp_part)) { found <- t; break } }
    
    if (is.na(found)) {
      if (grepl("C19_ctrl", comp_part)) return("C19_ctrl")
      if (grepl("Flu_ctrl", comp_part)) return("Flu_ctrl")
      if (grepl("HIV_ctrl", comp_part)) return("HIV_ctrl")
      if (grepl("OP_low", comp_part))   return("OP_low") 
    }
    return(ifelse(is.na(found), "Other", found))
  })
  
  ann_colors <- list(
    Cell_Type = cell_colors[names(cell_colors) %in% unique(annotation_df$Cell_Type)],
    Group = comp_group_colors,
    Condition = condition_colors[names(condition_colors) %in% unique(annotation_df$Condition)]
  )

  # Constrained Clustering (Within Group)
  ordered_groups <- c("HIV", "Flu", "C19", "OP")
  ordered_groups <- ordered_groups[ordered_groups %in% unique(annotation_df$Group)]
  
  final_order <- c()
  gap_indices <- c()
  current_idx <- 0
  
  for (g in ordered_groups) {
    group_items <- rownames(annotation_df)[annotation_df$Group == g]
    if (length(group_items) > 1) {
      sub_mat <- cor_mat[group_items, group_items]
      hc <- hclust(as.dist(1 - sub_mat), method = "ward.D2")
      final_order <- c(final_order, group_items[hc$order])
    } else {
      final_order <- c(final_order, group_items)
    }
    current_idx <- current_idx + length(group_items)
    gap_indices <- c(gap_indices, current_idx)
  }
  
  gap_indices <- gap_indices[-length(gap_indices)]
  cor_mat_ordered <- cor_mat[final_order, final_order]
  annotation_df_ordered <- annotation_df[final_order, ]

  # Plot
  h_size <- max(12, num_loaded * 0.4) 
  save_path <- file.path(plotDir, "Global_Correlation_Constrained.pdf")
  
  pdf(save_path, width = h_size, height = h_size)
  pheatmap(cor_mat_ordered,
           main = "Global Correlation (Constrained by Group)",
           annotation_col = annotation_df_ordered,
           annotation_row = annotation_df_ordered,
           annotation_colors = ann_colors,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           gaps_row = gap_indices,
           gaps_col = gap_indices,
           display_numbers = FALSE, 
           color = colorRampPalette(c("#053061", "#2166ac", "#f7f7f7", "#b2182b", "#67001f"))(100),
           breaks = seq(-sensitivity_limit, sensitivity_limit, length.out = 101),
           fontsize = 8,
           border_color = NA 
  )
  dev.off()
  message(paste("Saved:", save_path))
} else {
  message("Error: Not enough data found.")
}
