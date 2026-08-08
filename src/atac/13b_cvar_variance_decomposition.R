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
#    Sample     : fragment-file identifier (one per sample, NOT per person)
#    Donor      : subject identity (random effect; controls pseudoreplication)
# -------------------------------------------------------------------
meta$Cell_Type <- factor(meta$Cell_Type)
meta$Exposure  <- factor(meta$Exposure)
meta$Condition <- factor(meta$Condition)
meta$Sample    <- factor(meta$Sample)

# -------------------------------------------------------------------
# DONOR (subject) identity.
# meta$Sample is the fragment-file identifier, i.e. one level per SAMPLE, not
# per person. Several cohorts sample the same donor repeatedly:
#   HIV        3 timepoints per subject (hiv1..hiv12 -> 4 subjects)
#   Influenza  4 timepoints per subject (fluXX_1..fluXX_4)
#   COVID-19   some donors sampled twice (e.g. 555_1/555_2, 66D0/66D7)
# Using Sample as the "donor" term therefore does NOT control for repeated
# sampling of the same individual. We derive a true Donor ID below so the
# decomposition is genuinely donor-aware (reviewer R2 / R3.2).
# -------------------------------------------------------------------
hiv_subject <- c(
  hiv6 = "HIV_S1", hiv12 = "HIV_S1", hiv9 = "HIV_S1",
  hiv8 = "HIV_S2", hiv4  = "HIV_S2", hiv1 = "HIV_S2",
  hiv2 = "HIV_S3", hiv7  = "HIV_S3", hiv3 = "HIV_S3",
  hiv11 = "HIV_S4", hiv10 = "HIV_S4", hiv5 = "HIV_S4"
)

derive_donor <- function(sample_id) {
  s <- sub("_fragments\\.tsv\\.gz$", "", as.character(sample_id))
  s <- sub("^ATAC_", "", s)
  vapply(s, function(x) {
    if (!is.na(hiv_subject[x])) return(unname(hiv_subject[x])) # HIV timepoints
    if (grepl("^flu", x)) return(sub("_[0-9]+$", "", x))       # fluXX_1..4 -> fluXX
    # COVID-19 repeat sampling: 555_1/555_2, 564A/564B, 66D0/66D7, 132D0/132D7
    y <- sub("_[12]$", "", x)          # 555_1 -> 555
    y <- sub("([0-9])[AB]$", "\\1", y) # 564A  -> 564
    y <- sub("([0-9])D[0-9]+$", "\\1", y) # 66D0 -> 66, 132D0 -> 132
    y
  }, character(1), USE.NAMES = FALSE)
}
meta$Donor <- factor(derive_donor(meta$Sample))
message("Samples: ", nlevels(meta$Sample), " | Donors: ", nlevels(meta$Donor))
print(head(unique(data.frame(Sample = meta$Sample, Donor = meta$Donor)), 20))

# Drop pseudobulks with missing labels or singleton factor levels
keep <- complete.cases(meta[, c("Cell_Type","Exposure","Condition","Sample","Donor")])
zmat <- zmat[, keep, drop = FALSE]
meta <- meta[keep, , drop = FALSE]

# Remove motifs with zero variance (variancePartition will error otherwise)
motif_var <- apply(zmat, 1, var, na.rm = TRUE)
zmat <- zmat[is.finite(motif_var) & motif_var > 0, , drop = FALSE]

message("Motifs: ", nrow(zmat), " | pseudobulks: ", ncol(zmat))
message("Cell types: ", nlevels(droplevels(meta$Cell_Type)),
        " | Exposures: ", nlevels(droplevels(meta$Exposure)),
        " | Samples: ", nlevels(droplevels(meta$Sample)),
        " | Donors: ", nlevels(droplevels(meta$Donor)))

# -------------------------------------------------------------------
# 3. Variance partition
#    All factors as random effects (categorical). This asks: of the
#    variance in each TF's activity across pseudobulks, how much is
#    attributable to lineage vs exposure vs donor?
#    NOTE: Condition and Exposure are nested (Condition contains Exposure).
#    Fit them in SEPARATE models to avoid aliasing:
#      Model A (shared vs exposure vs donor):  ~ Cell_Type + Exposure + Donor
#      Model B (adds stage/severity):          ~ Cell_Type + Condition + Donor
# -------------------------------------------------------------------

## ---- Model A: lineage vs exposure vs donor ----
formA <- ~ (1 | Cell_Type) + (1 | Exposure) + (1 | Donor)
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
formB <- ~ (1 | Cell_Type) + (1 | Condition) + (1 | Donor)
vpB <- fitExtractVarPartModel(zmat, formB, meta)
vpB_df <- as.data.frame(vpB)
vpB_df$Motif <- rownames(vpB_df)
write.csv(vpB_df, file.path(save_dir, "varpart_per_motif_modelB.csv"), row.names = FALSE)

# -------------------------------------------------------------------
# 4. Plots (boxplots of variance explained per factor; scientific palette)
# -------------------------------------------------------------------ß
sci_pal <- c(
  "#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4",
  "#91D1C2", "#DC0000", "#7E6148", "#B09C85", "#3B4992", "#008B45",
  "#631879", "#008280", "#BB0021", "#5F559B"
)

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
    geom_boxplot(outlier.size = 0.4, alpha = 0.9, colour = "grey20") +
    scale_fill_manual(values = sci_pal) +
    theme_classic(base_size = 13) +
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
  filter(Exposure > Cell_Type) %>%
  arrange(desc(Exposure)) %>%
  select(Motif, Cell_Type, Exposure, Donor, Residuals)
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
  formC <- ~ (1 | Cell_Type) + (1 | Condition) + (1 | Donor)
  vpC <- fitExtractVarPartModel(zmat_c, formC, meta_c)
  vpC_df <- as.data.frame(vpC); vpC_df$Motif <- rownames(vpC_df)
  write.csv(vpC_df, file.path(save_dir, "varpart_per_motif_COVID_severity.csv"), row.names = FALSE)
  plot_vp_box(vpC_df, "COVID-19: lineage vs severity vs donor",
              "varpart_boxplot_COVID_severity.pdf")
}

# -------------------------------------------------------------------
# 7. Within-cell-type decomposition (comment 112: "do this within each
#    cell type rather than across all cell types"). For each cell type we
#    drop the lineage axis and ask how TF-activity variance splits between
#    exposure and donor WITHIN that lineage.
# -------------------------------------------------------------------
# Within a cell type there is one pseudobulk per SAMPLE, but several cohorts
# sample the same donor repeatedly (HIV 3 timepoints, influenza 4, some COVID 2),
# so the DONOR term has fewer levels than observations and remains estimable.
# Including it separates exposure-associated from inter-individual variation
# within each lineage; where the design does not support it we fall back to an
# exposure-only model.
cell_types <- levels(droplevels(meta$Cell_Type))
vp_list <- lapply(cell_types, function(ct) {
  idx <- meta$Cell_Type == ct
  m <- droplevels(meta[idx, ])
  z <- zmat[, idx, drop = FALSE]
  if (nlevels(m$Exposure) < 2 || ncol(z) <= nlevels(m$Exposure) + 1) {
    message("  skip ", ct, " (insufficient design)")
    return(NULL)
  }
  zv <- apply(z, 1, var, na.rm = TRUE)
  z <- z[is.finite(zv) & zv > 0, , drop = FALSE]
  use_donor <- nlevels(m$Donor) >= 2 && nlevels(m$Donor) < ncol(z)
  form <- if (use_donor) ~ (1 | Exposure) + (1 | Donor) else ~ (1 | Exposure)
  vp <- tryCatch(as.data.frame(fitExtractVarPartModel(z, form, m)),
    error = function(e) {
      message("  ", ct, ": donor term failed (", conditionMessage(e),
        "); falling back to exposure-only")
      use_donor <<- FALSE
      as.data.frame(fitExtractVarPartModel(z, ~ (1 | Exposure), m))
    }
  )
  if (is.null(vp$Donor)) vp$Donor <- NA_real_
  vp$Cell_Type <- ct
  vp$n_pseudobulk <- ncol(z)
  vp$donor_modelled <- use_donor
  vp
})
vp_within <- dplyr::bind_rows(vp_list) # per-motif, per-cell-type

# median fraction per factor within each cell type
within_ct <- vp_within %>%
  group_by(Cell_Type) %>%
  summarise(
    n_pseudobulk = dplyr::first(n_pseudobulk),
    donor_modelled = dplyr::first(donor_modelled),
    median_Exposure = median(Exposure),
    median_Donor = median(Donor),
    median_Residuals = median(Residuals),
    .groups = "drop"
  ) %>%
  arrange(desc(median_Exposure))
write.csv(within_ct, file.path(save_dir, "varpart_within_celltype.csv"), row.names = FALSE)
message("Within-cell-type variance partition (median fraction per factor):")
print(within_ct)

# Boxplot: distribution of exposure-associated variance per motif, by cell type
vp_within$Cell_Type <- factor(vp_within$Cell_Type, levels = within_ct$Cell_Type)
# Cell-type colours: identical to the UMAP / heatmap palette used throughout
# (see src/atac/13_cvar_analysis.R) so panels are directly comparable.
cell_type_colors <- c(
  "B_mem" = "#AE017E", "B_naive" = "#F768A1", "DC" = "#67000D",
  "Mono_CD14" = "#FE9929", "Mono_CD16" = "#CC4C02", "NK_CD16" = "#A65628",
  "Plasma" = "#A106BD", "T_mait" = "#41B6C4", "T_mem_CD4" = "#4292c6",
  "T_mem_CD8" = "#0074cc", "T_mix" = "#888FB5", "T_naive" = "#C7E9B4"
)
missing_cols <- setdiff(levels(vp_within$Cell_Type), names(cell_type_colors))
if (length(missing_cols)) {
  warning("cell types without an assigned colour: ", paste(missing_cols, collapse = ", "))
}

p_within <- ggplot(vp_within, aes(x = Cell_Type, y = Exposure, fill = Cell_Type)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.9, colour = "grey20") +
  scale_fill_manual(values = cell_type_colors) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(
    title = "Exposure-associated TF-activity variance within each cell type",
    x = NULL, y = "Fraction of variance explained by exposure"
  )
ggsave(file.path(save_dir, "varpart_within_celltype_boxplot.pdf"), p_within,
  width = 8, height = 5
)

message("variancePartition decomposition complete.")
