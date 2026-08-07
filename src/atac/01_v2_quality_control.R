#!/usr/bin/env Rscript

#####################################################################
# 01_v2_quality_control.R
# created on 2026-08-07 by Irem Gunduz
# Association of sample-level PCs with known covariates
# (Reviewer response: characterise variation w.r.t. potential confounders)
#
# Rationale
#  - We ask whether the low-dimensional structure of the ATAC data is
#    associated with cohort (exposure_type), Sex and Age.
#  - Age and Sex are recorded per SAMPLE and only for the OP and HIV
#    cohorts (see Table S1), so those two tests are run on that subset.
#  - The covariates are sample-level, therefore the test is performed at
#    the sample level: the cell-level embedding is aggregated to a
#    per-sample mean before PCA. Testing cell-level PCs against
#    sample-level covariates would be pseudoreplication and would
#    massively inflate significance.
#  - We run the test on BOTH the pre-Harmony IterativeLSI embedding and
#    the Harmony-corrected embedding (added in 02_cluster_and_batch.R).
#    Expectation: cohort is associated before correction (the bias we
#    account for via Harmony) and is attenuated afterwards, while Age and
#    Sex are not significantly associated with any PC in either case.
#####################################################################

## Load Libraries
suppressPackageStartupMessages({
  library(readxl)
  library(ArchR)
  library(dplyr)
  library(ggplot2)
})
set.seed(12) # set seed

## ---- Paths / parameters -----------------------------------------------------
repo_dir <- "/icbb/projects/igunduz/irem_github/exposure_atlas_manuscript"
suppTables <- file.path(repo_dir, "sample_annots/All_Supplementary_Tables.xlsx")
n_pc <- 10 # number of top PCs to test
embeddings_to_test <- c("IterativeLSI", "Harmony") # pre- and post-batch-correction
outputDir <- "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/"

# Load the ArchR project (already subsetted to non-BA samples)
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)

## ---- Sample-level covariate table (Table S1) --------------------------------
# NB: project$Sample is the arrow filename (e.g. "ATAC_055_fragments.tsv.gz"),
# which matches Table S1$arrow_name -- NOT Table S1$sampleId ("C19_mod_055").
# We therefore join the covariates on arrow_name.
covar <- readxl::read_excel(suppTables, sheet = "Table S1") %>%
  dplyr::select(sampleId, arrow_name, exposure_type, Age, Sex) %>%
  dplyr::mutate(
    exposure_type = factor(exposure_type),
    Sex = factor(Sex),
    Age = as.numeric(Age)
  )

## ---- Helper: association of one PC with one covariate ------------------------
# Returns variance explained (R^2), F-test p-value and the usable n.
# Works for both categorical (Sex, exposure_type) and continuous (Age) x.
pc_assoc <- function(y, x) {
  keep <- !is.na(y) & !is.na(x)
  if (is.factor(x)) x <- droplevels(x[keep]) else x <- x[keep]
  y <- y[keep]
  n <- length(y)
  # need >= 2 groups / non-constant predictor and enough residual df
  ok <- n >= 3 && (if (is.factor(x)) nlevels(x) >= 2 else stats::sd(x) > 0)
  if (!ok) {
    return(c(r2 = NA_real_, p = NA_real_, n = n))
  }
  fit <- stats::lm(y ~ x)
  fs <- summary(fit)$fstatistic
  p <- if (is.null(fs)) NA_real_ else stats::pf(fs[1L], fs[2L], fs[3L], lower.tail = FALSE)
  c(r2 = summary(fit)$r.squared, p = unname(p), n = n)
}

## ---- Helper: cohort-ADJUSTED association of one PC with one covariate --------
# Tests the covariate (x) while controlling for cohort, i.e. compares
#   lm(y ~ cohort)  vs  lm(y ~ cohort + x)
# via a nested F-test, and reports the PARTIAL R^2 (extra variance explained
# by x beyond cohort). This separates a genuine Age/Sex effect from the
# cohort confounding, since Age/Sex are only available for OP + HIV.
pc_assoc_adj <- function(y, x, cohort) {
  keep <- !is.na(y) & !is.na(x) & !is.na(cohort)
  y <- y[keep]
  cohort <- droplevels(cohort[keep])
  if (is.factor(x)) x <- droplevels(x[keep]) else x <- x[keep]
  n <- length(y)
  ok <- n >= 4 && (if (is.factor(x)) nlevels(x) >= 2 else stats::sd(x) > 0)
  if (!ok) {
    return(c(r2_partial = NA_real_, p = NA_real_, n = n))
  }
  # reduced model: cohort only (fall back to intercept if <2 cohorts present)
  reduced <- if (nlevels(cohort) >= 2) stats::lm(y ~ cohort) else stats::lm(y ~ 1)
  full <- stats::update(reduced, . ~ . + x)
  an <- stats::anova(reduced, full)
  p <- an[["Pr(>F)"]][2L]
  rss_r <- sum(stats::resid(reduced)^2)
  rss_f <- sum(stats::resid(full)^2)
  r2_partial <- (rss_r - rss_f) / rss_r
  c(r2_partial = r2_partial, p = unname(p), n = n)
}

## ---- Run per embedding ------------------------------------------------------
assoc_results <- list()

for (emb in embeddings_to_test) {
  # cell x dims embedding
  reddim <- getReducedDims(project, reducedDims = emb)

  # map each cell to its sample, then take the per-sample mean embedding
  cell_sample <- factor(getCellColData(project, select = "Sample", drop = TRUE))
  sums <- rowsum(reddim, group = cell_sample)
  n_per_sample <- as.numeric(table(cell_sample)[rownames(sums)])
  sample_emb <- sums / n_per_sample # samples x dims (pseudobulk embedding)

  # sample-level PCA
  pca <- stats::prcomp(sample_emb, center = TRUE, scale. = TRUE)
  var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
  k <- min(n_pc, ncol(pca$x))
  pcs <- pca$x[, seq_len(k), drop = FALSE]

  # align covariates to the PCA sample order (join on arrow_name)
  cv <- covar[match(rownames(pcs), covar$arrow_name), ]
  if (anyNA(cv$arrow_name)) {
    warning(sprintf(
      "%s: %d sample(s) in the embedding not found in Table S1",
      emb, sum(is.na(cv$arrow_name))
    ))
  }

  # test each PC against each covariate
  covariate_cols <- list(
    Cohort = cv$exposure_type, # all samples
    Sex = cv$Sex, # OP + HIV subset (rest NA)
    Age = cv$Age # OP + HIV subset (rest NA)
  )

  res <- do.call(rbind, lapply(seq_len(k), function(i) {
    do.call(rbind, lapply(names(covariate_cols), function(cn) {
      a <- pc_assoc(pcs[, i], covariate_cols[[cn]])
      # cohort-adjusted partial test (not meaningful for Cohort itself)
      adj <- if (cn == "Cohort") {
        c(r2_partial = NA_real_, p = NA_real_, n = a[["n"]])
      } else {
        pc_assoc_adj(pcs[, i], covariate_cols[[cn]], cv$exposure_type)
      }
      data.frame(
        embedding = emb,
        PC = paste0("PC", i),
        pc_index = i,
        var_explained = var_explained[i],
        covariate = cn,
        r2 = a[["r2"]], # marginal variance explained
        p = a[["p"]], # marginal p
        r2_partial_adjCohort = adj[["r2_partial"]], # variance beyond cohort
        p_adjCohort = adj[["p"]], # cohort-adjusted p
        n = a[["n"]],
        stringsAsFactors = FALSE
      )
    }))
  }))
  assoc_results[[emb]] <- res
}

assoc_df <- dplyr::bind_rows(assoc_results)

# BH-adjust p-values within each embedding x covariate family
# (both the marginal test and the cohort-adjusted test)
assoc_df <- assoc_df %>%
  dplyr::group_by(embedding, covariate) %>%
  dplyr::mutate(
    p_adj = p.adjust(p, method = "BH"),
    p_adjCohort_bh = p.adjust(p_adjCohort, method = "BH")
  ) %>%
  dplyr::ungroup()

## ---- Variance-weighted summary ----------------------------------------------
# Total sample-level variance associated with each covariate =
#   sum over PCs of (variance explained by covariate in that PC) x (PC's share
#   of total variance). Gives one interpretable number per covariate/embedding.
var_summary <- assoc_df %>%
  dplyr::group_by(embedding, covariate) %>%
  dplyr::summarise(
    total_var_marginal = sum(r2 * var_explained, na.rm = TRUE),
    total_var_adjCohort = sum(r2_partial_adjCohort * var_explained, na.rm = TRUE),
    .groups = "drop"
  )
message("Total sample-level variance associated with each covariate:")
print(as.data.frame(var_summary))
write.csv(var_summary,
  file = file.path(repo_dir, "figures/pc_covariate_variance_summary.csv"),
  row.names = FALSE
)

# write the results table
write.csv(assoc_df,
  file = file.path(repo_dir, "figures/pc_covariate_association.csv"),
  row.names = FALSE
)

## ---- Visualise: -log10(adjusted p) per PC x covariate -----------------------
p_assoc <- ggplot(
  assoc_df,
  aes(x = PC, y = covariate, fill = -log10(p_adj))
) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = ifelse(is.na(p_adj), "", sprintf("%.2f", r2))), size = 3) +
  facet_wrap(~embedding, ncol = 1) +
  scale_fill_gradient(
    low = "white", high = "#006400",
    na.value = "grey95", name = "-log10(adj. p)"
  ) +
  scale_x_discrete(limits = paste0("PC", seq_len(n_pc))) +
  labs(
    title = "Association of sample-level PCs with known covariates",
    subtitle = "tile label = variance explained (R^2); Sex/Age tested on OP+HIV subset",
    x = NULL, y = NULL
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(p_assoc,
  filename = file.path(repo_dir, "figures/pc_covariate_association.pdf"),
  width = 9, height = 6, units = "in", dpi = 300
)

# quick console summary of any significant confounder associations
message("Significant MARGINAL PC-covariate associations (adj. p < 0.05):")
print(subset(assoc_df, p_adj < 0.05,
  select = c(embedding, PC, covariate, r2, p_adj, n)
))

# Key check for the rebuttal: does Age/Sex survive once cohort is controlled for?
message("Significant COHORT-ADJUSTED associations for Age/Sex (adj. p < 0.05):")
print(subset(assoc_df, covariate %in% c("Age", "Sex") & p_adjCohort_bh < 0.05,
  select = c(embedding, PC, covariate, r2_partial_adjCohort, p_adjCohort_bh, n)
))

#####################################################################
