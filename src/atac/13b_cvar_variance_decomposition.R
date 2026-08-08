#!/usr/bin/env Rscript

#####################################################################
# 13b_cvar_variance_decomposition.R
# Add-on to 13_cvar_analysis.R
# Purpose (reviewer R3.3): move from descriptive, pairwise cross-exposure
# comparisons to a design-aware decomposition of TF-activity variance into
# shared (cell-type / lineage), exposure-specific, severity/stage, and
# donor components. 
#
# Output:
#   - varpart_per_motif.csv         : variance fractions per TF motif
#   - varpart_summary.csv           : median variance fraction per factor
#   - varpart_violin.pdf            : distribution of variance explained
#   - (optional) severity subset run for COVID-19 only
#####################################################################

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(variancePartition)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(chromVAR)
})
set.seed(12)

save_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript/sample_annots/cvar_diffs"

# -------------------------------------------------------------------
# 1. Load the pseudobulk chromVAR object built in 13_cvar_analysis.R
#    (condition x cell type x sample; Z-scores in assay 'z')
# -------------------------------------------------------------------
dev <- readRDS(file.path(save_dir, "atac_chromvar_deviations_exposure_celltype.rds"))

# Z-score matrix: rows = TF motifs, cols = pseudobulk samples
zmat <- SummarizedExperiment::assay(dev, "z")
meta <- as.data.frame(SummarizedExperiment::colData(dev))

# -------------------------------------------------------------------
# 2. Build a clean design table
#    Cell_Type  : lineage identity (shared axis)
#    Exposure   : cohort / stimulus (COVID/HIV/Influenza/OP)
#    Condition  : exposure + stage/severity (e.g. C19_sev, HIV_acu)
#    Sample     : donor  (random effect; controls pseudoreplication)
# -------------------------------------------------------------------
meta$Cell_Type <- factor(meta$Cell_Type)
meta$Exposure  <- factor(meta$Exposure)
meta$Condition <- factor(meta$Condition)
meta$Sample    <- factor(meta$Sample)

# Drop pseudobulks with missing labels or singleton factor levels
keep <- complete.cases(meta[, c("Cell_Type","Exposure","Condition","Sample")])
zmat <- zmat[, keep, drop = FALSE]
meta <- meta[keep, , drop = FALSE]

# Remove motifs with zero variance (variancePartition will error otherwise)
motif_var <- apply(zmat, 1, var, na.rm = TRUE)
zmat <- zmat[is.finite(motif_var) & motif_var > 0, , drop = FALSE]

message("Motifs: ", nrow(zmat), " | pseudobulks: ", ncol(zmat))
message("Cell types: ", nlevels(droplevels(meta$Cell_Type)),
        " | Exposures: ", nlevels(droplevels(meta$Exposure)),
        " | Donors: ", nlevels(droplevels(meta$Sample)))

# -------------------------------------------------------------------
# 3. Variance partition
#    All factors as random effects (categorical). This asks: of the
#    variance in each TF's activity across pseudobulks, how much is
#    attributable to lineage vs exposure vs donor?
#    NOTE: Condition and Exposure are nested (Condition contains Exposure).
#    Fit them in SEPARATE models to avoid aliasing:
#      Model A (shared vs exposure vs donor):  ~ Cell_Type + Exposure + Sample
#      Model B (adds stage/severity):          ~ Cell_Type + Condition + Sample
# -------------------------------------------------------------------

## ---- Model A: lineage vs exposure vs donor ----
formA <- ~ (1 | Cell_Type) + (1 | Exposure) + (1 | Sample)
vpA <- fitExtractVarPartModel(zmat, formA, meta)
vpA_df <- as.data.frame(vpA)
vpA_df$Motif <- rownames(vpA_df)
write.csv(vpA_df, file.path(save_dir, "varpart_per_motif_modelA.csv"), row.names = FALSE)

summ_A <- vpA_df %>%
  select(-Motif) %>%
  summarise(across(everything(), median)) %>%
  pivot_longer(everything(), names_to = "Factor", values_to = "Median_VarExplained")
write.csv(summ_A, file.path(save_dir, "varpart_summary_modelA.csv"), row.names = FALSE)
print(summ_A)

## ---- Model B: lineage vs condition(stage/severity) vs donor ----
formB <- ~ (1 | Cell_Type) + (1 | Condition) + (1 | Sample)
vpB <- fitExtractVarPartModel(zmat, formB, meta)
vpB_df <- as.data.frame(vpB)
vpB_df$Motif <- rownames(vpB_df)
write.csv(vpB_df, file.path(save_dir, "varpart_per_motif_modelB.csv"), row.names = FALSE)

# -------------------------------------------------------------------
# 4. Plots (violin of variance explained per factor)
# -------------------------------------------------------------------
plot_vp <- function(vp_df, title, file) {
  long <- vp_df %>%
    select(-Motif) %>%
    pivot_longer(everything(), names_to = "Factor", values_to = "VarExplained")
  # order factors by median
  ord <- long %>% group_by(Factor) %>%
    summarise(m = median(VarExplained), .groups = "drop") %>%
    arrange(desc(m)) %>% pull(Factor)
  long$Factor <- factor(long$Factor, levels = ord)
  p <- ggplot(long, aes(x = Factor, y = VarExplained, fill = Factor)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.8) +
    geom_boxplot(width = 0.12, outlier.size = 0.4, fill = "white") +
    theme_classic() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(title = title, x = NULL, y = "Fraction of variance explained")
  ggsave(file.path(save_dir, file), p, width = 7, height = 5)
}

plot_vp(vpA_df, "TF-activity variance: lineage vs exposure vs donor",
        "varpart_violin_modelA.pdf")
plot_vp(vpB_df, "TF-activity variance: lineage vs condition vs donor",
        "varpart_violin_modelB.pdf")

# Boxplot version of the same decomposition (Fabian: supplement figure)
plot_vp_box <- function(vp_df, title, file) {
  long <- vp_df %>%
    select(-Motif) %>%
    pivot_longer(everything(), names_to = "Factor", values_to = "VarExplained")
  ord <- long %>%
    group_by(Factor) %>%
    summarise(m = median(VarExplained), .groups = "drop") %>%
    arrange(desc(m)) %>%
    pull(Factor)
  long$Factor <- factor(long$Factor, levels = ord)
  p <- ggplot(long, aes(x = Factor, y = VarExplained, fill = Factor)) +
    geom_boxplot(outlier.size = 0.4, alpha = 0.85) +
    theme_classic() +
    theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(title = title, x = NULL, y = "Fraction of variance explained")
  ggsave(file.path(save_dir, file), p, width = 7, height = 5)
}
plot_vp_box(vpA_df, "TF-activity variance: lineage vs exposure vs donor",
            "varpart_boxplot_modelA.pdf")
plot_vp_box(vpB_df, "TF-activity variance: lineage vs condition vs donor",
            "varpart_boxplot_modelB.pdf")

# -------------------------------------------------------------------
# 5. Which motifs are the exception -- exposure-driven, not lineage-driven?
#    These are the candidate 'condition-specific' TFs the reviewer wants
#    identified systematically (rather than eyeballed from a heatmap).
# -------------------------------------------------------------------
exposure_driven <- vpA_df %>%
  filter(Exposure > Cell_Type) %>%           # exposure explains more than lineage
  arrange(desc(Exposure)) %>%
  select(Motif, Cell_Type, Exposure, Sample, Residuals)
write.csv(exposure_driven,
          file.path(save_dir, "exposure_driven_motifs_modelA.csv"),
          row.names = FALSE)
message("Exposure-driven motifs (Exposure > Cell_Type): ", nrow(exposure_driven))

# -------------------------------------------------------------------
# 6. OPTIONAL: severity-specific decomposition within COVID-19 only
#    Directly addresses 'associated with severity/stage' for the cohort
#    where severity is defined (COVID-19: ctrl/mild/mod/sev).
# -------------------------------------------------------------------
covid_mask <- grepl("^C19_", meta$Condition)
if (sum(covid_mask) > 10) {
  meta_c <- droplevels(meta[covid_mask, ])
  zmat_c <- zmat[, covid_mask, drop = FALSE]
  zvar_c <- apply(zmat_c, 1, var, na.rm = TRUE)
  zmat_c <- zmat_c[is.finite(zvar_c) & zvar_c > 0, , drop = FALSE]
  formC <- ~ (1 | Cell_Type) + (1 | Condition) + (1 | Sample)
  vpC <- fitExtractVarPartModel(zmat_c, formC, meta_c)
  vpC_df <- as.data.frame(vpC); vpC_df$Motif <- rownames(vpC_df)
  write.csv(vpC_df, file.path(save_dir, "varpart_per_motif_COVID_severity.csv"), row.names = FALSE)
  plot_vp(vpC_df, "COVID-19: lineage vs severity vs donor",
          "varpart_violin_COVID_severity.pdf")
}

# -------------------------------------------------------------------
# 7. Within-cell-type decomposition (comment 112: "do this within each
#    cell type rather than across all cell types"). For each cell type we
#    drop the lineage axis and ask how TF-activity variance splits between
#    exposure and donor WITHIN that lineage.
# -------------------------------------------------------------------
# NOTE: within a single cell type there is one pseudobulk per donor, so the
# donor random effect (1 | Sample) has as many levels as observations and is not
# estimable. We therefore drop Sample here and decompose exposure vs residual
# only (donor variance is captured in the across-cell-type Model A above).
cell_types <- levels(droplevels(meta$Cell_Type))
vp_list <- lapply(cell_types, function(ct) {
  idx <- meta$Cell_Type == ct
  m <- droplevels(meta[idx, ])
  z <- zmat[, idx, drop = FALSE]
  # need >= 2 exposures and more observations than exposure levels
  if (nlevels(m$Exposure) < 2 || ncol(z) <= nlevels(m$Exposure) + 1) {
    message("  skip ", ct, " (insufficient design)")
    return(NULL)
  }
  zv <- apply(z, 1, var, na.rm = TRUE)
  z <- z[is.finite(zv) & zv > 0, , drop = FALSE]
  vp <- as.data.frame(fitExtractVarPartModel(z, ~ (1 | Exposure), m))
  vp$Cell_Type <- ct
  vp$n_pseudobulk <- ncol(z)
  vp
})
vp_within <- dplyr::bind_rows(vp_list) # per-motif, per-cell-type

# median fraction per factor within each cell type
within_ct <- vp_within %>%
  group_by(Cell_Type) %>%
  summarise(
    n_pseudobulk = dplyr::first(n_pseudobulk),
    median_Exposure = median(Exposure),
    median_Residuals = median(Residuals),
    .groups = "drop"
  ) %>%
  arrange(desc(median_Exposure))
write.csv(within_ct, file.path(save_dir, "varpart_within_celltype.csv"), row.names = FALSE)
message("Within-cell-type variance partition (median fraction per factor):")
print(within_ct)

# Boxplot: distribution of exposure-associated variance per motif, by cell type
vp_within$Cell_Type <- factor(vp_within$Cell_Type, levels = within_ct$Cell_Type)
p_within <- ggplot(vp_within, aes(x = Cell_Type, y = Exposure, fill = Cell_Type)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.85) +
  theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(
    title = "Exposure-associated TF-activity variance within each cell type",
    x = NULL, y = "Fraction of variance explained by exposure"
  )
ggsave(file.path(save_dir, "varpart_within_celltype_boxplot.pdf"), p_within,
  width = 8, height = 5
)

message("variancePartition decomposition complete.")
