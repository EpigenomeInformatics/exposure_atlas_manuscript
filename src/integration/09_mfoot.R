#!/usr/bin/env Rscript

#####################################################################
# 09_mfoot.R
# created on 10-11-2024 by Irem Gunduz
# Compare TF footprints between T cells using pseudobulks
#####################################################################

## Load Libraries
set.seed(42)
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(methylTFR)
  library(muLogR)
  library(muRtools)
  library(methylTFRAnnotationHg38)
  library(ggplot2)
  library(GenomicRanges)
})


# Set the paths
fig_dir <- "/icbb/projects/igunduz/methylTFR_081124/Tcell_fp_250625_divided/"
if (!dir.exists(fig_dir)) dir.create(fig_dir)
motifset <- "jaspar2020_distal"
tfset <- if (motifset == "jaspar2020_distal") {
  "jaspar2020"
} else {
  motifset
}
gcfreqs <- getGCfreq(motifSet = motifset)
tf_bindsites <- getTFbindsites(motifSet = tfset)
gc_dist <- getGenomeGC()
enhancer <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
filtered_names <- names(tf_bindsites)
flankNorm <- 50

# Define the paths to the methylome data for Th and Tc samples
th_mem_sample <- "/icbb/projects/igunduz/new_scmeth_pseudobulks_070224/Th-Mem.bedGraph"
th_naive_sample <- "/icbb/projects/igunduz/new_scmeth_pseudobulks_070224/Th-Naive.bedGraph"
tc_naive_sample <- "/icbb/projects/igunduz/new_scmeth_pseudobulks_070224/Tc-Naive.bedGraph"
tc_mem_sample <- "/icbb/projects/igunduz/new_scmeth_pseudobulks_070224/Tc-Mem.bedGraph"

# Define sample lists and names
groups <- list(
  Th = list(
    Th_Mem = list(sample = th_mem_sample, name = "Th-Mem"),
    Th_Naive = list(sample = th_naive_sample, name = "Th-Naive")
  ),
  Tc = list(
    Tc_Mem = list(sample = tc_mem_sample, name = "Tc-Mem"),
    Tc_Naive = list(sample = tc_naive_sample, name = "Tc-Naive")
  )
)

# Load the methylome data for each sample in groups
groups <- lapply(groups, function(group) {
  lapply(group, function(s) {
    logger.info(paste("Loading sample", s$name))
    s$sample <- read_methylome(s$sample, type = "bismarkcytosine", 5)
    return(s) # Return the updated list item with the loaded sample data
  })
})

# Function to generate and save a single plot for all groups (Th and Tc) with observed substracted expected
plot_and_save_difference <- function(groups, save_dir, obs_colors) {
  for (motif in filtered_names) {
    if (!file.exists(file.path(save_dir, paste0("TF_footprint_diff_", motif, "_all_groups.pdf")))) {
      logger.start(paste("Processing motif", motif, "for all groups"))

      # Generate plot data for all groups and calculate the difference
      combined_data <- rbindlist(lapply(names(groups), function(group_name) {
        group_samples <- groups[[group_name]]
        rbindlist(lapply(names(group_samples), function(name) {
          plot_data <- plotExpectedFootprint(
            motif, tf_bindsites, group_samples[[name]]$sample,
            sample_name = group_samples[[name]]$name,
            gc_dist = gc_dist, gcfreqs = gcfreqs,
            enhancer = enhancer, returnPlotData = TRUE
          )

          # Extract observed and expected data directly from plotDF
          difference_data <- plot_data$plotDF[, .(avg_methyl = avg_methyl[type == "Observed"] / avg_methyl[type == "Expected"]), by = x]
          difference_data[, type := paste("Observed divided Expected", group_samples[[name]]$name)]

          # Now normalize the ratio by flanking region
          flank <- max(abs(difference_data$x), na.rm = TRUE)
          idx <- abs(difference_data$x) >= flank - flankNorm
          norm_factor <- mean(difference_data$avg_methyl[idx], na.rm = TRUE)
          difference_data[, avg_methyl := avg_methyl / norm_factor]

          return(difference_data)
        }))
      }))
      logger.completed()

      logger.info("Plotting combined difference data for all groups")

      # Plot the combined difference data for all groups
      p_combined <- ggplot(combined_data, aes(x = x, y = avg_methyl, color = type)) +
        geom_line() + # Add lines for each type
        # geom_point() +  # Add points for each type
        xlab("Distance from motif center") +
        ylab("Methylation difference (Observed /Expected)") +
        theme_classic() +
        ggtitle(paste("TF footprint difference for", motif, "in all groups")) +
        scale_color_manual(values = obs_colors) + # Use the observed color vector
        theme(legend.position = "bottom")

      # Save the plot as a PDF in the specified directory
      ggsave(
        filename = file.path(save_dir, paste0("TF_footprint_diff_", motif, "_all_groupsVS2.pdf")),
        plot = p_combined, width = 12, height = 8
      )
    }
  }
}


# Define observed colors for all groups (since expected is subtracted, we don't need exp_colors)
obs_colors <- c(
  "Observed divided Expected Th-Mem" = "#41B6C4",
  "Observed divided Expected Th-Naive" = "#C7E9B4",
  "Observed divided Expected Tc-Mem" = "#4292C6",
  "Observed divided Expected Tc-Naive" = "#888FB5"
)

# Generate and save plots for all groups with observed substracted expected
plot_and_save_difference(groups, fig_dir, obs_colors)

#####################################################################
