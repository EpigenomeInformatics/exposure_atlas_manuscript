# Dissecting the epigenome dynamics in human immune cells upon viral and chemical exposure by multimodal single-cell profiling

Supporting repository for the manuscript of the same name. 

### 📁 Project Directory Structure

```text
figures/              # Output figures from the analysis
sample_annots/        # Sample-level annotations and metadata
data/                 # methylTFR objects
src/
├── atac/             # scATAC-seq analysis scripts
├── meth/             # snmC-seq analysis scripts
├── integration/      # scATAC-seq + snmC-seq integrative analysis scripts

utils/                # Utility functions shared across pipelines
```

### 📦 Dependencies

| Package           | Description                                     |
|-------------------|-------------------------------------------------|
| R ≥ 4.1           | Minimum required R version                      |
| RnBeads           | DNA methylation analysis                        |
| ChrAccR           | Chromatin accessibility analysis                |
| ArchR             | Single-cell ATAC-seq analysis framework         |
| dplyr             | Data manipulation and transformation            |
| data.table        | Data manipulation and transformation            |
| methylTFR         | Methylation based TF activities                 |
| chromVAR         | Accessibility based TF activities                 |
| ggplot2           | Data visualization                              |
| ComplexHeatmap    | Complex heatmaps with annotations               |

## 🧬 Getting Started

We separated the analysis into three main categories:

---

## 1. 🔓 ATAC: Chromatin Accessibility Analysis

This section includes single-cell and pseudobulk-based ATAC-seq analysis. The workflow is modular and organized into the following steps:

| Step | Script                        | Description |
|------|-------------------------------|-------------|
| 01   | `01_quality_control.R`        | Perform quality control on raw scATAC data |
| 02   | `02_cluster_and_batch.R`      | Handle clustering and batch correction |
| 03   | `03_annotate.R`               | Annotate cell-types |
| 04.1 | `04_1_markers.R`              | Plot cell type markers  |
| 04.2 | `04_2_cellprops.R`            | Analyze cell-type proportions across exposures |
| 05   | `05_pseudobulk.R`             | Perform pseudobulk aggregation per cell-type |
| 07.1 | `07_1_run_ChrAccR.R`          | Run ChrAccR analysis |
| 07.2 | `07_2_run_ChrAccR_C19.R`      | Run ChrAccR analysis focused on COVID-19 samples |
| 08.1 | `08_1_chraccR_plots.R`        | Generate visualizations from ChrAccR outputs |
| 08.2 | `08_2_C19_trackplots.R`       | Create genome track plots for COVID-19 differential peaks in CD14+ Monocytes |
| 09   | `09_gene_exp_vs_activity.R`   | Correlate gene expression with chromatin accessibility in CD14+ Monocytes for protein coding genes|
| 10   | `10_tcells.R` | T-cell subset analysis for longitudinal HIV cohort |
---

## 2. 🧬 METH: Methylation Analysis

This section includes pseudobulk methylation analysis focused on ATAC peaks.

| Step | Script                        | Description |
|------|-------------------------------|-------------|
| 01   | `01_meth_pseudobulks.R`       | Create pseudobulks for methylation data |
| 02.1   | `02_1_run_RnBeads.R`        | Run RnBeads analysis using the pseudobulks |
| 02.2   | `02_2_run_RnBeads_C19.R`    | Run RnBeads analysis for COVID-19 monocytes using the pseudobulks |
| 03   | `03_RnBeads_plots.R`          | Generate visualizations from RnBeads outputs  |
| 04   | `04_C19_pseudobulks.R`        | Generate pseudobulk per condition in C19 for mTFR visualizations  |
| 05   | `05_C19_mfoot.R`              | Generate motif footprint plots for C19 |
| 06   | `06_C19_mTFR.R`               | Run methylTFR algorithm to create deviation scores for C19 |
| 07   | `07_C19_mTFR_plots.R`         | Generate visualizations methylTFR deviation scores for C19 |
---

## 3. 🧩 Integration

This section includes integration of single-cell methylation and chromatin accessibility data from overlapping samples, based on shared ATAC peaks.

| Step | Script                  | Description                                                   |
|------|-------------------------|---------------------------------------------------------------|
| 01   | `01_prepare_sampleannot.R` | Format sample annotation ready for aggregation               |
| 02   | `02_aggregate_meth.R`      | Aggregate scMeth over peak regions                           |
| 03   | `03_quality_check.R`       | Perform quality control on aggregated data                   |
| 04   | `04_lsi.R`                 | Apply Latent Semantic Indexing (LSI) for dimensionality reduction |
| 05   | `05_cca.R`                 | Run Canonical Correlation Analysis (CCA) for multi-omic alignment |
| 06   | `06_plot_cca.R`            | Visualize results of CCA                                     |
| 07   | `07_mTFR_run.R`            | Run methylTFR algorithm to create deviation scores           |
| 08   | `08_cor_analysis.R`        | Correlation analysis of mTFR and cVAR matricies              |
| 09  | `09_mfoot.R`                | Generate motif footprint plots for T-cells                    |
| 10 |`10_zdiff.R`                | Z-score difference plots for T-cells                    |
## 📫 Contact

For questions or contributions, feel free to reach out to the [maintainer](https://github.com/igunduz).
