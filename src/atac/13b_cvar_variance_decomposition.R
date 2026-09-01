#!/usr/bin/env Rscript

#####################################################################
# 13b_cvar_variance_decomposition.R -- add-on to 13_cvar_analysis.R
# Reviewer R3.3: design-aware decomposition of TF-activity variance into
# lineage, exposure, severity/stage, donor and (for C19) processing-batch
# components. Outputs varpart_*.csv / .pdf into save_dir.
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

# Override the repo location with EXPOSURE_ATLAS_REPO=/path Rscript ...
repo_candidates <- c(
  Sys.getenv("EXPOSURE_ATLAS_REPO"),
  "/scratch/icbb/igunduz/irem_github/exposure_atlas_manuscript",
  "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript"
)
repo_candidates <- repo_candidates[nzchar(repo_candidates)]
repo_dir <- repo_candidates[which(dir.exists(repo_candidates))[1L]]
if (is.na(repo_dir)) stop("Set EXPOSURE_ATLAS_REPO to the repository path.")
annot_dir <- file.path(repo_dir, "sample_annots")
save_dir  <- file.path(annot_dir, "cvar_diffs")

# 1. Load the pseudobulk chromVAR object from 13_cvar_analysis.R (Z in assay 'z')
dev <- readRDS(file.path(save_dir, "atac_chromvar_deviations_exposure_celltype.rds"))

# Z-score matrix: rows = TF motifs, cols = pseudobulk samples
zmat <- SummarizedExperiment::assay(dev, "z")
meta <- as.data.frame(SummarizedExperiment::colData(dev))

# 2. Design table: Cell_Type (lineage), Exposure (cohort), Condition (exposure +
#    stage/severity), Sample (fragment file), Donor (subject)
meta$Cell_Type <- factor(meta$Cell_Type)
meta$Exposure  <- factor(meta$Exposure)
meta$Condition <- factor(meta$Condition)
meta$Sample    <- factor(meta$Sample)

# DONOR identity. meta$Sample is one level per SAMPLE, not per person: HIV has 3
# timepoints/subject, influenza 4, and some C19 donors were sampled twice.
hiv_subject <- c(
  hiv6 = "HIV_S1", hiv12 = "HIV_S1", hiv9 = "HIV_S1",
  hiv8 = "HIV_S2", hiv4  = "HIV_S2", hiv1 = "HIV_S2",
  hiv2 = "HIV_S3", hiv7  = "HIV_S3", hiv3 = "HIV_S3",
  hiv11 = "HIV_S4", hiv10 = "HIV_S4", hiv5 = "HIV_S4"
)

# Patterns match anywhere in the identifier, whatever prefix the column carries.
derive_donor <- function(sample_id) {
  x <- as.character(sample_id)
  out <- x # default: one donor per sample (OP cohort)

  # HIV: hiv<N> anywhere -> subject ID (3 timepoints per subject)
  has_hiv <- grepl("hiv[0-9]+", x)
  if (any(has_hiv)) {
    key <- sub(".*?(hiv[0-9]+).*", "\\1", x[has_hiv])
    mapped <- unname(hiv_subject[key])
    mapped[is.na(mapped)] <- key[is.na(mapped)] # unknown hiv id -> keep as-is
    out[has_hiv] <- mapped
  }

  # Influenza: flu<ID>_<timepoint> -> flu<ID> (4 timepoints per subject)
  has_flu <- grepl("flu[0-9]+_[0-9]+", x)
  if (any(has_flu)) {
    out[has_flu] <- sub(".*?(flu[0-9]+)_[0-9]+.*", "\\1", x[has_flu])
  }

  # C19: 555_1/555_2 -> 555 ; 564A/564B -> 564 ; 66D0/66D7 -> 66 ; 132D0 -> 132
  has_c19 <- grepl("ATAC_.+_fragments", x) & !has_flu
  if (any(has_c19)) {
    id <- sub(".*ATAC_(.+)_fragments.*", "\\1", x[has_c19])
    id <- sub("_[12]$", "", id)
    id <- sub("([0-9])[AB]$", "\\1", id)
    id <- sub("([0-9])D[0-9]+$", "\\1", id)
    out[has_c19] <- paste0("C19_", id)
  }
  out
}
meta$Donor <- factor(derive_donor(meta$Sample))

message("Samples: ", nlevels(meta$Sample), " | Donors: ", nlevels(meta$Donor))
if (nlevels(meta$Donor) == nlevels(meta$Sample)) {
  warning("Donor did not collapse any samples -- check the Sample ID format printed below")
}
collapsed <- unique(data.frame(
  Sample = as.character(meta$Sample), Donor = as.character(meta$Donor)
))
message("Donors contributing more than one sample:")
print(head(subset(
  merge(collapsed, as.data.frame(table(collapsed$Donor)),
    by.x = "Donor", by.y = "Var1"
  ), Freq > 1
), 30))

# PROCESSING BATCH -- C19 only (710/720), so section 6 uses it; A and B cannot.
batch_path <- file.path(annot_dir, "ATAC_metadata_covid.csv")
meta$Batch <- factor(rep(NA_character_, nrow(meta)))
if (file.exists(batch_path)) {
  amc <- read.csv(batch_path, stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(amc$arrow_name))
  sid <- as.character(meta$Sample)
  has <- grepl("ATAC_.+_fragments", sid)
  lab <- rep(NA_character_, length(sid))
  lab[has] <- paste0("ATAC_", sub(".*ATAC_(.+)_fragments.*", "\\1", sid[has]))
  meta$Batch <- factor(amc$processing_date[match(lab, amc$arrow_name)])
  message("Processing batch attached to ", sum(!is.na(meta$Batch)), " of ",
    nrow(meta), " pseudobulks; levels: ",
    paste(levels(meta$Batch), collapse = ", "))
} else {
  warning("ATAC_metadata_covid.csv not found -- batch will not be modelled.")
}

# Drop incomplete rows (Batch excluded: it is NA for three of four cohorts)
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

# 3. Variance partition, all factors as random effects.
#    Condition contains Exposure, so they are fitted in SEPARATE models:
#      A ~ Cell_Type + Exposure + Donor      B ~ Cell_Type + Condition + Donor

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

# 4. Plots: variance explained per factor
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

# 5. Motifs that are exposure-driven rather than lineage-driven
exposure_driven <- vpA_df %>%
  filter(Exposure > Cell_Type) %>%
  arrange(desc(Exposure)) %>%
  select(Motif, Cell_Type, Exposure, Donor, Residuals)
write.csv(exposure_driven,
          file.path(save_dir, "exposure_driven_motifs_modelA.csv"),
          row.names = FALSE)
message("Exposure-driven motifs (Exposure > Cell_Type): ", nrow(exposure_driven))

# 6. COVID-19 only: severity, and the one cohort where batch is recorded
covid_mask <- grepl("^C19_", meta$Condition)
if (sum(covid_mask) > 10) {
  meta_c <- droplevels(meta[covid_mask, ])
  zmat_c <- zmat[, covid_mask, drop = FALSE]
  zvar_c <- apply(zmat_c, 1, var, na.rm = TRUE)
  zmat_c <- zmat_c[is.finite(zvar_c) & zvar_c > 0, , drop = FALSE]

  ## ---- Is severity confounded with processing batch? ----
  has_batch <- !all(is.na(meta_c$Batch))
  if (has_batch) {
    don <- unique(data.frame(
      Donor = as.character(meta_c$Donor),
      Condition = as.character(meta_c$Condition),
      Batch = as.character(meta_c$Batch), stringsAsFactors = FALSE
    ))
    don <- don[!is.na(don$Batch), , drop = FALSE]
    bt <- table(Condition = don$Condition, Batch = don$Batch)
    message("COVID-19 severity arm vs processing batch (donor level):")
    print(bt)
    write.csv(as.data.frame(bt),
      file.path(save_dir, "COVID_severity_vs_batch.csv"), row.names = FALSE)
    if (all(dim(bt) >= 2)) {
      ft <- tryCatch(fisher.test(bt, simulate.p.value = TRUE, B = 1e5),
                     error = function(e) NULL)
      if (!is.null(ft)) {
        message(sprintf("Fisher exact test of severity vs batch: p = %.4g%s",
          ft$p.value,
          if (ft$p.value < 0.05)
            "  <-- NOT independent; the Condition fraction in model C is partly batch"
          else "  (balanced)"))
      }
    }
  }

  ## ---- Model C: lineage vs severity vs donor ----
  formC <- ~ (1 | Cell_Type) + (1 | Condition) + (1 | Donor)
  vpC <- fitExtractVarPartModel(zmat_c, formC, meta_c)
  vpC_df <- as.data.frame(vpC); vpC_df$Motif <- rownames(vpC_df)
  write.csv(vpC_df, file.path(save_dir, "varpart_per_motif_COVID_severity.csv"), row.names = FALSE)
  plot_vp_box(vpC_df, "COVID-19: lineage vs severity vs donor",
              "varpart_boxplot_COVID_severity.pdf")

  ## ---- Model C2: processing batch in place of donor ----
  # Batch cannot be ADDED to C: every sample of a donor shares a batch, so
  # (1|Donor) would absorb it. Swapping donor out identifies batch instead.
  # If Condition holds up in C2 the severity signal is not merely batch; if it
  # drops while Batch takes a comparable share, the two are not separable.
  if (has_batch && nlevels(droplevels(meta_c$Batch)) >= 2) {
    keep_b  <- !is.na(meta_c$Batch)
    meta_c2 <- droplevels(meta_c[keep_b, ])
    zmat_c2 <- zmat_c[, keep_b, drop = FALSE]
    zv2     <- apply(zmat_c2, 1, var, na.rm = TRUE)
    zmat_c2 <- zmat_c2[is.finite(zv2) & zv2 > 0, , drop = FALSE]

    formC2  <- ~ (1 | Cell_Type) + (1 | Condition) + (1 | Batch)
    vpC2    <- fitExtractVarPartModel(zmat_c2, formC2, meta_c2)
    vpC2_df <- as.data.frame(vpC2); vpC2_df$Motif <- rownames(vpC2_df)
    write.csv(vpC2_df,
      file.path(save_dir, "varpart_per_motif_COVID_severity_batch.csv"),
      row.names = FALSE)
    plot_vp_box(vpC2_df, "COVID-19: lineage vs severity vs processing batch",
                "varpart_boxplot_COVID_severity_batch.pdf")

    # compare on the motifs both models retained
    shared <- intersect(vpC_df$Motif, vpC2_df$Motif)
    i1 <- match(shared, vpC_df$Motif); i2 <- match(shared, vpC2_df$Motif)
    cmp <- data.frame(
      Factor = c("Cell_Type", "Condition", "Donor (C) / Batch (C2)", "Residuals"),
      Median_modelC = c(median(vpC_df$Cell_Type[i1]), median(vpC_df$Condition[i1]),
                        median(vpC_df$Donor[i1]),     median(vpC_df$Residuals[i1])),
      Median_modelC2 = c(median(vpC2_df$Cell_Type[i2]), median(vpC2_df$Condition[i2]),
                         median(vpC2_df$Batch[i2]),     median(vpC2_df$Residuals[i2]))
    )
    write.csv(cmp, file.path(save_dir, "varpart_COVID_severity_C_vs_C2.csv"),
      row.names = FALSE)
    message("COVID-19 median variance fractions -- model C (donor) vs C2 (batch), ",
      length(shared), " shared motifs:")
    print(cmp)
    message(sprintf(
      "Condition: %.3f with donor -> %.3f with batch (%+.3f). Batch median: %.3f",
      cmp$Median_modelC[2], cmp$Median_modelC2[2],
      cmp$Median_modelC2[2] - cmp$Median_modelC[2], cmp$Median_modelC2[3]))
  }
}

# 7. Within-cell-type decomposition (comment 112): drop the lineage axis and
#    split variance between exposure and donor within each lineage.
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
# cell-type colours as in 13_cvar_analysis.R, so panels are comparable
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

# facet by component so the cell-type colour scheme matches the UMAPs
vp_within_long <- vp_within %>%
  dplyr::select(Cell_Type, Exposure, Donor) %>%
  tidyr::pivot_longer(c(Exposure, Donor),
    names_to = "Component", values_to = "VarExplained"
  ) %>%
  dplyr::filter(!is.na(VarExplained)) %>%
  dplyr::mutate(Component = factor(Component,
    levels = c("Exposure", "Donor"),
    labels = c("Exposure", "Donor (inter-individual)")
  ))

p_within <- ggplot(vp_within_long, aes(x = Cell_Type, y = VarExplained, fill = Cell_Type)) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.9, colour = "grey20") +
  facet_wrap(~Component, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = cell_type_colors) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(
    title = "TF-activity variance within each cell type",
    x = NULL, y = "Fraction of variance explained"
  )
ggsave(file.path(save_dir, "varpart_within_celltype_boxplot.pdf"), p_within,
  width = 8, height = 7
)

# exposure-only version (single panel), in case the figure has room for one row
p_within_exp <- ggplot(
  subset(vp_within_long, Component == "Exposure"),
  aes(x = Cell_Type, y = VarExplained, fill = Cell_Type)
) +
  geom_boxplot(outlier.size = 0.4, alpha = 0.9, colour = "grey20") +
  scale_fill_manual(values = cell_type_colors) +
  theme_classic(base_size = 13) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1)) +
  labs(
    title = "Exposure-associated TF-activity variance within each cell type",
    x = NULL, y = "Fraction of variance explained by exposure"
  )
ggsave(file.path(save_dir, "varpart_within_celltype_exposureonly.pdf"), p_within_exp,
  width = 8, height = 5
)

message("variancePartition decomposition complete.")
