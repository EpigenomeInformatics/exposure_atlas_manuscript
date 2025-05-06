# Dissecting the Epigenome Dynamics of Human Immune Cells Upon Viral and Chemical Exposure at Single-Cell Resolution

Supporting repository for the manuscript of the same name. 

### 📁 Project Directory Structure

```text
figures/              # Output figures from the analysis
sample_annots/        # Sample-level annotations and metadata

src/
├── atac/             # scATAC-seq analysis scripts
├── meth/             # snmC-seq analysis scripts
├── integration/      # scATAC-seq + snmC-seq integrative analysis scripts

utils/                # Utility functions shared across pipelines
```

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
| 06   | `06_C19_mTFR.R`               | Generate methylTFR deviation scores for C19 |
| 07   | `07_C19_mTFR_plots.R`         | Generate visualizations methylTFR deviation scores for C19 |

---

## 3. 🧩 Integration 

This section includes integration of single-cell methylation and chromatin accessibility data from overlapping samples, based on shared ATAC peaks.


## 📫 Contact

For questions or contributions, feel free to reach out to the maintainer.