#!/usr/bin/env Rscript

#####################################################################
# c19_mfoot_080425.R
# Created on 08-04-2025 by Irem Gunduz
# Compare TF footprints in COVID-19 samples for Jaspar2020 motifs
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
plot_dir <- "/icbb/projects/igunduz/C19_mFOOT_140725_mono1_jaspar2020_v2/"
if (!dir.exists(plot_dir)) dir.create(plot_dir)
motifset <- "jaspar2020_distal"

# Load annotations
gcfreqs <- getGCfreq(motifSet = motifset)
tfset <- if (motifset == "jaspar2020_distal") {
  "jaspar2020"
} else {
  motifset
}
tf_bindsites <- getTFbindsites(motifSet = tfset)
gc_dist <- readRDS("/icbb/projects/igunduz/annotation/genome_wide_GC.RDS")
enhancer <- readRDS("/icbb/projects/share/annotations/methylTFRAnnotationHg38/inst/extdata/distal_regions.RDS")
filtered_names <- names(tf_bindsites)
flankNorm <- 50

# Define the paths to the methylome data for COVID-19 samples
covid_samples <- list(
  CommercialControl_healthy = list(sample = "/icbb/projects/igunduz/new_covid19_pseudobulks_180924/CommercialControl_healthy.bedGraph", name = "CommercialControl_healthy"),
  COVID_severe_mono_1 = list(sample = "/icbb/projects/igunduz/new_covid19_pseudobulks_080425/severe_mono1.bedGraph", name = "COVID_severe_mono_1") # ,
)

# Load the methylome data for each sample
covid_samples <- lapply(covid_samples, function(s) {
  logger.info(paste("Loading sample", s$name))
  s$sample <- read_methylome(s$sample, type = "bismarkcytosine", 1)
  return(s)
})

publication <- c("SPIB","CEBPD","RELA","CREB1","CEBPA","JUN","BATF","REL","BATF::JUN","SPIC","FOSB::JUN")
# Function to generate and save a plot for all COVID-19 samples
plot_and_save_difference_covid <- function(samples, save_dir, obs_colors) {
  for (motif in filtered_names) {
    if (!file.exists(file.path(save_dir, paste0("TF_footprint_diff_", motif, "_covid.pdf")))) {
      logger.start(paste("Processing motif", motif, "for COVID-19 samples"))

      # Generate plot data and calculate the difference
      combined_data <- rbindlist(lapply(names(samples), function(name) {
        plot_data <- plotExpectedFootprint(
          motif, tf_bindsites, samples[[name]]$sample,
          sample_name = samples[[name]]$name,
          gc_dist = gc_dist, gcfreqs = gcfreqs,
          enhancer = enhancer, returnPlotData = TRUE
        )
        # Calculate observed/expected methylation
        difference_data <- plot_data$plotDF[, .(avg_methyl = avg_methyl[type == "Observed"] / avg_methyl[type == "Expected"]), by = x]
        difference_data[, type := paste("Observed divided Expected", samples[[name]]$name)]

        # Now normalize the ratio by flanking region
        flank <- max(abs(difference_data$x), na.rm = TRUE)
        idx <- abs(difference_data$x) >= flank - flankNorm
        norm_factor <- mean(difference_data$avg_methyl[idx], na.rm = TRUE)
        difference_data[, avg_methyl := avg_methyl / norm_factor]
        
        return(difference_data)
      }))
      logger.completed()

      # Plot the combined difference data for all samples
      logger.info("Plotting combined difference data for all COVID-19 samples")
      p_combined <- ggplot(combined_data, aes(x = x, y = avg_methyl, color = type)) +
        geom_line() + # Add lines for each type
        # geom_point() +  # Add points for each type
        xlab("Distance from motif center") +
        ylab("Methylation difference (Observed / Expected)") +
        theme_classic() +
        ggtitle(paste("TF footprint difference for", motif, "in COVID-19 samples")) +
        scale_color_manual(values = obs_colors) + # Use the observed color vector
        theme(legend.position = "bottom") +
        xlim(-200, 200) # Set the x-axis limits

      # Save the plot as a PDF in the specified directory
      ggsave(
        filename = file.path(save_dir, paste0("TF_footprint_diff_", motif, "_covid.pdf")),
        plot = p_combined, width = 12, height = 8
      )
    }
  }
}

# Define the observed colors for COVID-19 samples
obs_colors <- c(
  "Observed divided Expected CommercialControl_healthy" = "#A8CEE3",
  "Observed divided Expected COVID_severe_mono_1" = "#D95F02"
)

# Generate and save plots for all samples with observed divided expected
plot_and_save_difference_covid(covid_samples, plot_dir, obs_colors)

#####################################################################

