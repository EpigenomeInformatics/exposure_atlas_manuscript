suppressPackageStartupMessages({
  library(dplyr)
  library(data.table)
  library(ggplot2)
  library(stringr)
  library(cowplot)
  library(ggpubr)
})

# --- Configuration ---

# 1. Main ArchR Project Directory (Where the standard cell folders are)
archrDir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_2023-10-02/"

# 2. COVID Results Directory 
# (Update this path to point to the specific folder containing your COVID .tsv files)
covidDir <- "/path/to/your/separate/covid_results_folder/" 

# 3. Output Plot Directory
plotDir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/plots/correlations_scatter/"

if (!dir.exists(plotDir)) {
  dir.create(plotDir, recursive = TRUE)
}

# 4. Cells and Comparisons
cells <- c("B_mem", "B_naive", "Mono_CD14", "Mono_CD16", "NK_CD16", 
           "T_mem_CD8", "T_mem_CD4", "T_naive", "T_mix", "T_mait")

target_comparisons <- c(
  "HIV_ctrl vs HIV_chr", 
  "HIV_ctrl vs HIV_acu",
  "Influenza_d3 vs Influenza_ctrl",
  "Influenza_d30 vs Influenza_ctrl",
  "OP_high vs OP_low"
)

# --- Helper Function ---
format_comp_name <- function(name) {
  new_name <- str_replace_all(name, "Influenza", "Flu")
  new_name <- str_replace_all(new_name, " vs ", "\nvs\n") 
  return(new_name)
}

# --- Main Loop ---

for (cell in cells) {
  message(paste0("Processing: ", cell))
  
  lfc_list <- list()
  
  # ==========================================
  # PART A: Load Standard ArchR Results
  # ==========================================
  compTablePath <- file.path(archrDir, cell, "reports/differential_data/comparisonTable.rds")
  
  if (file.exists(compTablePath)) {
    comparisonTable <- readRDS(compTablePath)
    comparisonTable$comps <- paste0(comparisonTable$grp1, " vs ", comparisonTable$grp2)
    matched_indices <- which(comparisonTable$comps %in% target_comparisons)
    
    for (i in matched_indices) {
      raw_comp_name <- comparisonTable$comps[i]
      diff_file <- file.path(archrDir, cell, "reports/differential_data", paste0("diffTab_", i, "_archrPeaks.tsv"))
      
      if (file.exists(diff_file)) {
        dt <- fread(diff_file, select = c("chrom", "chromStart", "chromEnd", "log2FoldChange"))
        dt[, peakID := paste(chrom, chromStart, chromEnd, sep = "_")]
        
        nice_name <- format_comp_name(raw_comp_name)
        tmp_lfc <- dt[, .(peakID, log2FoldChange)]
        setnames(tmp_lfc, "log2FoldChange", nice_name)
        
        lfc_list[[nice_name]] <- tmp_lfc
      }
    }
  }

  # ==========================================
  # PART B: Load Separate COVID Results
  # ==========================================
  
  # Find a file in the covidDir that contains the cell name
  # e.g., "Mono_CD14_C19_vs_Ctrl.tsv"
  covid_files_found <- list.files(covidDir, pattern = cell, full.names = TRUE)
  
  # Filter for .tsv or .csv if necessary
  covid_files_found <- covid_files_found[grep("\\.tsv$|\\.csv$", covid_files_found)]
  
  if (length(covid_files_found) > 0) {
    # Take the first match
    c19_path <- covid_files_found[1] 
    message(paste0("  Found COVID file: ", basename(c19_path)))
    
    # Read file (Adapting to potential column names)
    c19_dt <- fread(c19_path)
    
    # Ensure we have the coordinate columns to build peakID
    # Adjust "chrom", "start", "end" strings below if your COVID file headers differ
    if (all(c("chrom", "chromStart", "chromEnd", "log2FoldChange") %in% colnames(c19_dt))) {
      
      c19_dt[, peakID := paste(chrom, chromStart, chromEnd, sep = "_")]
      c19_name <- "COVID\nvs\nCtrl"
      
      tmp_c19 <- c19_dt[, .(peakID, log2FoldChange)]
      setnames(tmp_c19, "log2FoldChange", c19_name)
      
      lfc_list[[c19_name]] <- tmp_c19
      
    } else {
      message("  Error: COVID file columns must include: chrom, chromStart, chromEnd, log2FoldChange")
    }
  } else {
    message(paste0("  Warning: No COVID file found for ", cell))
  }

  # ==========================================
  # PART C: Merge, Filter, and Plot
  # ==========================================
  
  # Need at least 2 comparisons to make a scatter plot
  if (length(lfc_list) < 2) next
  
  # Merge all data.tables in the list by peakID
  merged_lfc <- Reduce(function(x, y) merge(x, y, by = "peakID", all = TRUE), lfc_list)
  
  # Create Matrix for Variance Calculation
  mat <- as.matrix(merged_lfc[, -1, with = FALSE])
  row_vars <- apply(mat, 1, var, na.rm = TRUE)
  
  # Select Top 5000 Variable Peaks
  n_keep <- min(5000, length(row_vars))
  keep_idx <- order(row_vars, decreasing = TRUE)[1:n_keep]
  plot_data <- merged_lfc[keep_idx, ]
  
  # Generate Plots
  plot_list <- list()
  col_names <- colnames(plot_data)[-1] 
  pairs <- combn(col_names, 2, simplify = FALSE)
  
  for (p in pairs) {
    x_var <- p[1]
    y_var <- p[2]
    
    all_vals <- c(plot_data[[x_var]], plot_data[[y_var]])
    limit_range <- range(all_vals, na.rm = TRUE)
    
    g <- ggplot(plot_data, aes_string(x = paste0("`", x_var, "`"), y = paste0("`", y_var, "`"))) +
      geom_point(alpha = 0.3, size = 0.8, color = "#2c3e50") +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red", size = 0.5) +
      geom_smooth(method = "lm", color = "blue", size = 0.5, se = FALSE) +
      stat_cor(method = "pearson",
               aes(label = ..r.label..),
               label.x.npc = "left", label.y.npc = "top",
               size = 5, color = "red") +
      theme_bw() +
      coord_fixed(xlim = limit_range, ylim = limit_range) +
      labs(x = paste0("LFC: ", str_replace_all(x_var, "\n", " ")),
           y = paste0("LFC: ", str_replace_all(y_var, "\n", " "))) +
      theme(axis.title = element_text(size = 9))
    
    plot_list[[paste(x_var, y_var)]] <- g
  }
  
  # Save Grid
  if (length(plot_list) > 0) {
    final_grid <- plot_grid(plotlist = plot_list, ncol = 3, align = 'hv')
    save_path <- paste0(plotDir, "Scatter_Dots_", cell, ".pdf")
    n_rows <- ceiling(length(plot_list) / 3)
    ggsave(save_path, final_grid, width = 12, height = 4 * n_rows, limitsize = FALSE)
    message(paste("Saved:", save_path))
  }
}