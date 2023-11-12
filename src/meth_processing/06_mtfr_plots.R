suppressPackageStartupMessages({
  library(dplyr)
  library(methylTFR)
  library(circlize)
  library(ComplexHeatmap)
  library(chromVAR)
  library(ChrAccR)
  library(ggrepel)
  library(pheatmap)
  library(ggplot2)
  library(muLogR)
  library(rescueR)
 # library(JASPAR2020)
#  library(TFBSTools)
  library(ArchR)
  library(muRtools)
})
set.seed(42)

logger.info("Starting analysis for Monocyte in chromVAR and methylTFR")
mtfr_dir <- "/icbb/projects/igunduz/DARPA_analysis/methyltfr_041023/updated_jaspar2020_091123"
#mtfr_dir <- "/icbb/projects/igunduz/DARPA_analysis/methyltfr_041023/updatedClustResults_071023"
cell <- "Monocyte"
condition <- c("C19_ctrl", "C19_sev")
ds_dir <- "/icbb/projects/igunduz/DARPA_analysis/chracchr_run_011023/ChrAccRuns_covid_2023-10-02/Mono_CD14/data/"
source("/icbb/projects/igunduz/sc_epigenome_exp/utils/mtfr_plots_helpers.R")
source("/icbb/projects/igunduz/exposure_atlas_manuscript/src/utils/methyltfr_utils.R")
source("/icbb/projects/igunduz/exposure_atlas_manuscript/src/utils/mtfr_plots_helpers.R")
source("/icbb/projects/igunduz/sc_epigenome_exp/utils/chraccr.R")

logger.start("Reading methylTFR deviations")
mtfr_devs <- list.files(mtfr_dir, pattern = cell, full.names = TRUE)
# filter conditions from the list
mtfr_devs <- rlist::list.cbind(lapply(condition, function(x) {
  readRDS(mtfr_devs[grepl(x, mtfr_devs)])
}))

mtfr_devs <- as.data.frame(mtfr_devs)
logger.completed()

#tf_bindsites <- readRDS("/icbb/projects/igunduz/annotation/methylTFRAnnotationHg38/inst/extdata/jaspar2020_tf_bindsites.rds")
#rownames(mtfr_devs) <- names(tf_bindsites)
logger.info("Reading sample annotation")
sannot <- data.table::fread(paste0("/icbb/projects/igunduz/DARPA/Generated/methylTFR/bed/Monocyte/C19_sev_vs_Ctrl/sample_methylation_summary.tsv"))
groups <- sannot$C19_sev_vs_Ctrl
logger.start("Computing differential deviations")
diffm <- computeL2FCdevs_vs2(deviations = as.matrix(rev(mtfr_devs)), computezscore = TRUE, group = groups,parametric=T)
logger.completed()


logger.start("Starting chromVAR analysis")
ds <- ChrAccR::loadDsAcc(paste0(ds_dir, "dsATAC_raw/"))
samples <- ChrAccR::getSamples(ds)
samples <- samples[grep("C19_ctrl|C19_sev", samples)]
ds <- ds[samples] # subset the samples

# logger.start("Loading altius motifs")
# cvMot <- methylTFRAnnotationHg38::getTFbindsites()
# logger.completed()

# setTimeOut(10000)
# cvd <- ChrAccR:::computeDeviations_altius(ds, "archr_peaks")
# setConfigElement("chromVarMotifs","jaspar")
# cvMot <- getConfigElement("chromVarMotifs")
motifPFMatrixList <- getMatrixSet(
  x = JASPAR2020,
  opts = list(species = 9606, all_versions = FALSE, collection = "CORE")
)
# Get the motif names
motifPFMatrixList <- organizeMotifNames(motifPFMatrixList)

cvd <- ChrAccR:::getChromVarDev(ds, "archr_peaks", motifs = motifPFMatrixList)
logger.info(c("Saving chromvar object: ", paste0(ds_dir, "chromvar_diffs_jaspar2020.rds")))
# saveRDS(cvd,paste0(ds_dir,"chromvar_diffs_mono_clust.rds"))
saveRDS(cvd, paste0(ds_dir, "chromvar_diffs_jaspar2020_rawds.rds"))
logger.start("Running differential analysis")
chromvar_mat <- chromVAR::deviations(cvd)
chromvar_mat <- as.data.frame(chromvar_mat)
# rownames(chromvar_mat) <- sub(".*_", "", rownames(chromvar_mat))
groups <- ifelse(grepl("ctrl", colnames(chromvar_mat)), "C19_ctrl", "C19_sev")
diff <- computeL2FCdevs_vs2(deviations = as.matrix(chromvar_mat), computezscore = TRUE,parametric=TRUE,
 #grp1name = "C19_ctrl", grp2name="C19_sev",
 group = groups)#,chromvar_obj=cvd)


diff$isDiff_1 <- ifelse(diff$p_value_adjusted < 0.05, TRUE, FALSE)
diffm$isDiff_2 <- ifelse(diffm$p_value_adjusted < 0.05, TRUE, FALSE)

# merge the two matrices
diff <- diff[, c("motifs", "zDiff", "isDiff_1")]
colnames(diff) <- c("name", "zDiff_chromvar", "isDiff_1")
diffm <- diffm[, c("motifs", "zDiff", "isDiff_2")]
colnames(diffm) <- c("name", "zDiff_methylTFR", "isDiff_2")
merged <- merge(diff, diffm, by = "name")

merged_orig <- merged
# organise motifnames in merged
merged$name <- sub("_.*", "", merged$name)

#organize motif names
merged$name <- sub("\\.var.*", "", merged$name)
merged$name <- sub(".*\\.\\.", "", merged$name)

#set threshold for chromvar 
#merged$zDiff_chromvar[abs(merged$zDiff_chromvar) > ] <- NA
#merged <- merged[complete.cases(merged), ]

# Call the plotScatterL2FC function with the merged data
l2fcmc <- plotScatterL2FC(
  datatable = merged,
  y_lab = "zDiff_methylTFR",
  x_lab = "zDiff_chromvar",
  comb = "Comparison",
  group1 = "zDiff_chromvar",
  group2 = "zDiff_methylTFR",
  label = TRUE,
  max_overlaps = 14# ,
  # textsize= 20,
  # bins=50
)
cor(merged$zDiff_chromvar, merged$zDiff_methylTFR)
logger.info("Plotting scatter plot")
ggsave(paste0("/icbb/projects/igunduz/Figures/chromvar_mtfr_scatter_mono_jaspar2020.pdf"), l2fcmc, width = 10, height = 10)

merged <- merged_orig
#subset merged based on diffs
merged <- merged[merged$isDiff_1 | merged$isDiff_2, ]
#subset based on negative zdiff on chromvar positive on methylTFR
merged_an <- merged[merged$zDiff_chromvar < 0 & merged$zDiff_methylTFR > 0, ]
#subset based on positive zdiff on chromvar negative on methylTFR
merged_an2 <- merged[merged$zDiff_chromvar > 0 & merged$zDiff_methylTFR < 0, ]
#merge them together
merged_an <- rbind(merged_an, merged_an2)
motifs <- merged_an$name

#sannot <- data.table::fread(paste0("/icbb/projects/igunduz/DARPA/Generated/methylTFR/bed/Monocyte/C19_sev_vs_Ctrl/sample_methylation_summary.tsv"))
#groups <- sannot$C19_sev_vs_Ctrl
#diffm <- computeL2FCdevs_vs2(deviations = as.matrix(rev(mtfr_devs)), computezscore = TRUE, group = groups)
#diffm$isDiff_2 <- ifelse(diffm$p_value_adjusted < 0.05, TRUE, FALSE)
#motifs <- rownames(diffm)[diffm$isDiff_2]

logger.start("Plotting the heatmap for methylTFR")
# mtfr_devs <- mtfr_devs[motifs,]
mtfr_devs <- methylTFR:::computeZScore(as.matrix(mtfr_devs))
mtfr_devs <- as.data.frame(mtfr_devs)
# rownames(mtfr_devs) <- names(tf_bindsites)
logger.info("Create a data frame for the samples' conditions")
ann <- data.frame(Condition = ifelse(grepl("Ctrl", colnames(mtfr_devs)), "Control", "Severe"))
rownames(ann) <- colnames(mtfr_devs)
mtfr_devs <- mtfr_devs[motifs, ]
hm_mtfr <- deviationHeatmap(mtfr_devs, ann, cluster_rows = TRUE)
logger.completed()

logger.info("Plotting the heatmap for methylTFR") 
pdf(paste0(ds_dir, "mtfr_diff_pheatmap.pdf"), width = 20, height = 20)
draw(hm_mtfr$hm)
dev.off()



logger.start("Plotting the heatmap for chromVAR")
logger.info("Subset chromVAR matrix for differentials")
chromvar_mat <- chromVAR::deviationScores(cvd)#computeZScore(chromVAR::deviations(cvd))
chromvar_mat <- chromvar_mat[motifs, ]
chromvar_mat <- as.data.frame(chromvar_mat)
logger.info("Create a data frame for the samples' conditions")
# Create a data frame for the samples' conditions
ann <- data.frame(Condition = ifelse(grepl("ctrl", colnames(chromvar_mat)), "Control", "Severe"))
rownames(ann) <- colnames(chromvar_mat)
hm <- deviationHeatmap(chromvar_mat, ann, package = "chromVAR", colors = "cb.BrBG", row_order = hm_mtfr$row_order, cluster_rows = FALSE)

# plot the chromvar heatmap
pdf(paste0(ds_dir, "chromvar_diff_pheatmap.pdf"), width = 20, height = 20)
draw(hm$hm)
dev.off()
logger.completed()


logger.start("Running analysis on cell-type level with chromVAR")
conditions <- c("HIV_acu", "HIV_chr", "HIV_ctrl", "OP_high", "OP_low", "OP_med")
# Define color vectors
exposure_colors <- c(
  "HIV_acu" = "#DC143C",
  "HIV_chr" = "#800080",
  "HIV_ctrl" = "#FFC0CB",
  "OP_high" = "#A0522D",
  "OP_low" = "#5E3D23",
  "OP_med" = "#D2691E"
)
at_cell_cols <- c(
  "B_mem" = "#AE017E",
  "B_naive" = "#F768A1",
  "DC" = "#67000D",
  "Mono_CD14" = "#FE9929",
  "Mono_CD16" = "#CC4C02",
  "NK_CD16" = "#A65628",
  "Plasma" = "#A106BD",
  "T_mait" = "#41B6C4",
  "T_mem_CD4" = "#4292c6",
  "T_mem_CD8" = "#0074cc",
  "T_mix" = "#888FB5",
  "T_naive" = "#C7E9B4"
)

outputDir <- "/icbb/projects/igunduz/archr_project_011023/"
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)

# Subset arhr project for conditions
idxSample <- BiocGenerics::which(project$sample_exposure_group %in% conditions)
cellsSample <- project$cellNames[idxSample]
project <- project[cellsSample, ]

logger.start("Looking into variable z-scores using chromVAR")
pseudo_chrom <- ArchR::getGroupSE(project, "altiusMatrix", groupBy = "ClusterCellTypes", divideN = T)
seZ <- pseudo_chrom[rowData(pseudo_chrom)$seqnames == "z", ]
zmat <- assay(seZ)
# zmat <- computeZScore(zmat)
zmat <- as.matrix(zmat)
rownames(zmat) <- rowData(seZ)$name
# zmat <- zmat[motifs,]
logger.completed()


logger.info("Preparing annotation for chromVAR")
ann <- data.frame(Cell = colnames(zmat))
# ann <- unique(ann)
rownames(ann) <- ann$Cell
# ann <- dplyr::select(ann, -cell_sample)
# ann <- ann[colnames(zmat), ]

zmat[abs(round(zmat, 2)) < 0.1] <- NA
zmat[abs(round(zmat, 2)) > 5] <- NA
zmat <- zmat[complete.cases(zmat), ]
motifs <- rownames(zmat)

heatmap_result <- deviationW2annotHeatmap(
  mat = zmat,
  ann_df = ann,
  cluster_col = "Cell",
  fill_col_cell_type = at_cell_cols, # Assuming you have colors in ann$CellType
  # fill_col_condition = exposure_colors, # Assuming you have colors in ann$Condition
  colors = "cptcity.arendal_temperature",
  clustering_method = "between",
  cluster_rows = TRUE,
  package = "chromVAR",
  row_order = NULL
)


# Plot the heatmap
pdf("/icbb/projects/igunduz/chromvar_heatmap_zscores.pdf", width = 20, height = 20)
draw(heatmap_result$hm)
dev.off()
logger.completed()

# archr_mat <- getMatrixFromProject(project,"altiusMatrix")
ccd <- getCellColData(project)
ccd$cell_sample <- paste0(ccd$ClusterCellTypes, "_", ccd$Sample)

# add sample explosure group
project <- addCellColData(
  ArchRProj = project, data = ccd$cell_sample,
  name = "pseudo_group", cells = project$cellNames, force = T
)

pseudo_chrom <- ArchR::getGroupSE(project, "altiusMatrix", groupBy = "pseudo_group", divideN = F)
seZ <- pseudo_chrom[rowData(pseudo_chrom)$seqnames == "deviations", ]
mat <- assay(seZ)
rownames(mat) <- rowData(seZ)$name


logger.info("Preparing annotation for chromVAR")
ann <- data.frame(cell_sample = ccd$cell_sample, Condition = ccd$sample_exposure_group, Cell = ccd$ClusterCellTypes)
ann <- unique(ann)
rownames(ann) <- ann$cell_sample
ann <- dplyr::select(ann, -cell_sample)
ann <- ann[colnames(mat), ]

heatmap_result <- deviationW2annotHeatmap(
  mat = mat[motifs, ],
  ann_df = ann,
  cluster_col = "Cell",
  fill_col_cell_type = at_cell_cols, # Assuming you have colors in ann$CellType
  fill_col_condition = exposure_colors, # Assuming you have colors in ann$Condition
  colors = "cptcity.arendal_temperature",
  clustering_method = "within",
  cluster_rows = TRUE,
  package = "chromVAR",
  row_order = NULL
)


# Plot the heatmap
pdf("/icbb/projects/igunduz/chromvar_heatmap.pdf", width = 20, height = 20)
draw(heatmap_result$hm)
dev.off()
logger.completed()

logger.start("Running analysis on cell-type level with methylTFR")
cells <- c("B-cell", "Monocyte", "NK-cell", "Tc-Mem", "Tc-Naive", "Th-Mem", "Th-Naive")

# Replace forward slashes with underscores in conditions
conditions <- gsub("/", "_", conditions)

# Create all possible combinations of cells and conditions
combinations <- expand.grid(Cell = cells, Condition = conditions)

# Create the desired strings
result <- paste0(combinations$Cell, "_", combinations$Condition, "_")

mtfr_dir <- "/icbb/projects/igunduz/DARPA_analysis/methyltfr_041023/ClustResults"

# Initialize data frames to store results
cell <- condition <- list()

# Loop through combinations and read the data
mtfr_devs <- list()
for (x in seq_along(result)) {
  path <- paste0(mtfr_dir, "/", result[x], "deviations.RDS")
  cur_dev <- readRDS(path)
  mtfr_devs[[x]] <- cur_dev
  cell[[x]] <- rep(as.character(combinations$Cell[[x]]), ncol(cur_dev))
  condition[[x]] <- rep(as.character(combinations$Condition[[x]]), ncol(cur_dev))
}

# Combine the lists into data frames
mtfr_devs <- do.call(cbind, mtfr_devs)
cell <- unlist(cell)
condition <- unlist(condition)

logger.info("Successfully read all the methylTFR deviations")
ann <- data.frame(Condition = condition, Cell = cell)
rownames(ann) <- colnames(mtfr_devs)

logger.start("Adding color columns to ann based on Condition and Cell")

cell_type_colors <- c(
  "B-cell" = "#AE017E",
  "Monocyte" = "#CC4C02",
  "NK-cell" = "#A65628",
  "Th-Mem" = "#41B6C4",
  "Tc-Mem" = "#4292C6",
  "Tc-Naive" = "#888FB5",
  "Th-Naive" = "#C7E9B4"
)

# ann$ExposureColor <- exposure_colors[ann$Condition]
# ann$CellTypeColor <- cell_type_colors[ann$Cell]

heatmap_result <- deviationW2annotHeatmap(
  mat = mtfr_devs[motifs, ],
  ann_df = ann,
  cluster_col = "Cell",
  fill_col_cell_type = cell_type_colors, # Assuming you have colors in ann$CellType
  fill_col_condition = exposure_colors, # Assuming you have colors in ann$Condition
  colors = "cptcity.arendal_temperature",
  clustering_method = "within",
  cluster_rows = FALSE,
  row_order = motifs
)


# Plot the heatmap
pdf("/icbb/projects/igunduz/mtfr_heatmap.pdf", width = 20, height = 20)
draw(heatmap_result$hm)
dev.off()

logger.completed()

logger.start("Running z-score analysis using methylTFR")

#logger.info("Dividing the matrix by the number of cells in each cluster")
#nCells <- table(cell)
#groupMat <- t(mtfr_devs) / as.vector(nCells)
#groupMat <- as.data.frame(t(groupMat))

groupMat <- mtfr_devs
#groupMat <- computeZScore(as.matrix(groupMat))
groupMat <- as.data.frame(t(groupMat))

groupMat$Cell <- cell
groupMat <- aggregate(. ~ Cell, data = groupMat, FUN = sum)
groupMat <- as.data.frame(t(groupMat))
colnames(groupMat) <- groupMat[1, ]
groupMat <- groupMat[-1, ]
groupMat <- apply(groupMat, 2, as.numeric)
rownames(groupMat) <- rownames(mtfr_devs)
logger.completed()
groupMat <- computeZScore(as.matrix(groupMat))

logger.info("Successfully read all the methylTFR deviations")
ann <- data.frame(Cell = colnames(groupMat))
rownames(ann) <- colnames(groupMat)

heatmap_result <- deviationW2annotHeatmap(
  mat = groupMat[motifs, ],
  ann_df = ann,
  cluster_col = "Cell",
  fill_col_cell_type = cell_type_colors, # Assuming you have colors in ann$CellType
  # fill_col_condition = exposure_colors, # Assuming you have colors in ann$Condition
  colors = "cptcity.arendal_temperature",
  clustering_method = "within",
  cluster_rows = TRUE,
  row_order = motifs
)

# Plot the heatmap
pdf("/icbb/projects/igunduz/mtfr_heatmap_zscore.pdf", width = 20, height = 20)
draw(heatmap_result$hm)
dev.off()
logger.completed()

logger.start("Correlation analysis between chromVAR and methylTFR")
sub_mtfr <- groupMat # [motifs,]
chromvar <- data.frame(
  B_cell = zmat[, "B_mem"],#rowMeans(zmat[, c(1:2)]),
  Monocyte = rowMeans(zmat[, c(4:5)]),
  NK_cell = zmat[, "NK_CD16"],
  Tc_Mem = zmat[, "T_mem_CD8"],
  T_naive = zmat[, "T_naive"],
  Th_mem = zmat[, "T_mem_CD4"]
)
methyltfr <- data.frame(
  B_cell = sub_mtfr[, "B-cell"],
  Monocyte = sub_mtfr[, "Monocyte"],
  NK_cell = sub_mtfr[, "NK-cell"],
  Tc_Mem = sub_mtfr[, "Tc-Mem"],
  T_naive = rowMeans(data.frame(sub_mtfr[, 5], sub_mtfr[, 7])),
  Th_mem = sub_mtfr[, "Th-Mem"]
)
methyltfr <- as.data.frame(t(methyltfr))
chromvar <- as.data.frame(t(chromvar))
cormat <- cor(chromvar, methyltfr)
# cormat <- apply(cormat,2,as.numeric)


d <- data.frame(cor = diag(cormat))
rownames(d) <- colnames(chromvar)
neg_cor <- dplyr::filter(d, cor <= -0.5)
c <- apply(chromvar, 2, as.numeric)
m <- apply(methyltfr, 2, as.numeric)
c <- c[, rownames(neg_cor)]
m <- m[, rownames(neg_cor)]

pvals <- sapply(colnames(c), function(motif) p.adjust(cor.test(c[, motif], m[, motif])$p.value))
pvals <- data.frame(vals = pvals)
p <- ggplot(d, aes(cor)) +
  geom_histogram(color = "black", fill = "white")
ggsave(p, file = "/icbb/projects/igunduz/correlation_hist.pdf")

ml <- c
mr <- m
rownames(ml) <- rownames(mr) <- rownames(chromvar)
mr <- t(mr)
ml <- t(ml)

cell_type <- c(
  "B_cell" = "#AE017E",
  "Monocyte" = "#CC4C02",
  "NK_cell" = "#A65628",
  "Th_mem" = "#41B6C4",
  "Tc_Mem" = "#4292C6",
  "T_naive" = "#888FB5"
)

#cs <- colpal.cont(9, "cb.Reds")
#colors.cv <- ChrAccR::getConfigElement("colorSchemesCont")
#colors.cv <- colors.cv[[".default.div"]]
#csl <- circlize::colorRamp2(seq(round(min(ml))-1, round(max(ml))+1, length.out = length(colors.cv )), colors.cv )

cs_col <- colpal.cont(9,  "cptcity.jjg_misc_temperature")
#cs <- circlize::colorRamp2(seq(round(min(mr)), round(max(mr)), length.out = length(cs_col )), cs_col )
#csl <- circlize::colorRamp2(seq(round(min(mr)), round(max(mr)), length.out = length(cs_col )), cs_col )

cs <- circlize::colorRamp2(seq(round(min(mr)), round(max(mr)), length.out = length(cs_col )), cs_col )
csl <- circlize::colorRamp2(seq(round(min(mr)), round(max(mr)), length.out = length(cs_col )), cs_col )

ann <- data.frame(Cell = colnames(mr))
rownames(ann) <- colnames(mr)

heatmap_result <- deviationW2annotHeatmap(
  mat = mr,
  ann_df = ann,
  cluster_col = "Cell",
  fill_col_cell_type = cell_type, # Assuming you have colors in ann$CellType
  # fill_col_condition = exposure_colors, # Assuming you have colors in ann$Condition
  colors = "cptcity.arendal_temperature",
  clustering_method = "within",
  cluster_rows = FALSE#,
  #row_order = motifs
)

# Plot the heatmap
pdf("/icbb/projects/igunduz/mtfr_heatmap_zscore_vs2.pdf", width = 20, height = 20)
draw(heatmap_result$hm)
dev.off()

heatmap_result <- deviationW2annotHeatmap(
  mat = ml,
  ann_df = ann,
  cluster_col = "Cell",
  fill_col_cell_type = cell_type, # Assuming you have colors in ann$CellType
  # fill_col_condition = exposure_colors, # Assuming you have colors in ann$Condition
  colors = "cptcity.arendal_temperature",
  clustering_method = "within",
  cluster_rows = FALSE,
  package = "chromVAR"
  #row_order = motifs
)

# Plot the heatmap
pdf("/icbb/projects/igunduz/cvar_heatmap_zscore_vs2.pdf", width = 20, height = 20)
draw(heatmap_result$hm)
dev.off()

colAnnot <- HeatmapAnnotation(
  Cells = colnames(mr),
  col = list(
    Cells = cell_type
  ),
  show_annotation_name = FALSE
)

hmap_cm <- diagDivCellHeatmap(ml, mr,
  top_annotation = colAnnot,
  name.l = "chromVAR",
  name.r = "methylTFR",
  ,col.r= cs,col.l = csl)
#) # ,col.r = csr,col.l = cs)

pdf("/icbb/projects/igunduz/correlation_heatmap.pdf", width = 20, height = 20)
draw(hmap_cm)
dev.off()
logger.completed()

logger.start("Preparing figure 1 chromVAR heatmap...")
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)

logger.start("Looking into variable z-scores using chromVAR")
pseudo_chrom <- ArchR::getGroupSE(project, "altiusMatrix", groupBy = "ClusterCellTypes", divideN = TRUE)
seZ <- pseudo_chrom[rowData(pseudo_chrom)$seqnames == "z", ]
zmat <- assay(seZ)
zmat <- as.matrix(zmat)
rownames(zmat) <- rowData(seZ)$name
logger.completed()


logger.info("Preparing annotation for chromVAR")
ann <- data.frame(Cell = colnames(zmat))
rownames(ann) <- ann$Cell

zmat[abs(round(zmat, 2)) < 0.1] <- NA
zmat[abs(round(zmat, 2)) > 3] <- NA
zmat <- zmat[complete.cases(zmat), ]
motifs <- rownames(zmat)

cmat <- deviationW2annotHeatmap(
  mat = zmat,
  ann_df = ann,
  cluster_col = "Cell",
  fill_col_cell_type = at_cell_cols, # Assuming you have colors in ann$CellType
  # fill_col_condition = exposure_colors, # Assuming you have colors in ann$Condition
  colors = "cptcity.arendal_temperature",
  clustering_method = "within",
  cluster_rows = TRUE,
  package = "chromVAR",
  row_order = NULL
)


# Plot the heatmap
pdf("/icbb/projects/igunduz/figure_1_chromvar_heatmap_zscores.pdf", width = 20, height = 20)
draw(cmat$hm)
dev.off()
logger.completed()

logger.start("Looking for methylTFR values of same motifs")
conditions <- c("HIV_acu", "HIV_chr", "HIV_ctrl", "OP_high", "OP_low", "OP_med", "C19_ctrl", "C19_sev", "C19_mild", "Influenza_d30", "Influenza_ctrl")

# Replace forward slashes with underscores in conditions
conditions <- gsub("/", "_", conditions)

# Create all possible combinations of cells and conditions
combinations <- expand.grid(Cell = cells, Condition = conditions)

# Create the desired strings
result <- paste0(combinations$Cell, "_", combinations$Condition, "_")

# Initialize data frames to store results
cell <- condition <- list()

# Loop through combinations and read the data
mtfr_devs <- list()
for (x in seq_along(result)) {
  path <- paste0(mtfr_dir, "/", result[x], "deviations.RDS")
  cur_dev <- readRDS(path)
  mtfr_devs[[x]] <- cur_dev
  cell[[x]] <- rep(as.character(combinations$Cell[[x]]), ncol(cur_dev))
  condition[[x]] <- rep(as.character(combinations$Condition[[x]]), ncol(cur_dev))
}

# Combine the lists into data frames
mtfr_devs <- do.call(cbind, mtfr_devs)
cell <- unlist(cell)
condition <- unlist(condition)


groupMat <- mtfr_devs
groupMat <- computeZScore(as.matrix(groupMat))
groupMat <- as.data.frame(t(groupMat))

groupMat$Cell <- cell
groupMat <- aggregate(. ~ Cell, data = groupMat, FUN = median)
groupMat <- as.data.frame(t(groupMat))
colnames(groupMat) <- groupMat[1, ]
groupMat <- groupMat[-1, ]
groupMat <- apply(groupMat, 2, as.numeric)
rownames(groupMat) <- rownames(mtfr_devs)
groupMat <- groupMat[motifs, ]


logger.info("Successfully read all the methylTFR deviations")
ann <- data.frame(Cell = colnames(groupMat))
rownames(ann) <- colnames(groupMat)

heatmap_result <- deviationW2annotHeatmap(
  mat = groupMat[motifs, ],
  ann_df = ann,
  cluster_col = "Cell",
  fill_col_cell_type = cell_type_colors, # Assuming you have colors in ann$CellType
  # fill_col_condition = exposure_colors, # Assuming you have colors in ann$Condition
  colors = "cptcity.arendal_temperature",
  clustering_method = "within",
  cluster_rows = TRUE
)

# Plot the heatmap
pdf("/icbb/projects/igunduz/figure_1_mtfr_heatmap_zscore.pdf", width = 20, height = 20)
draw(heatmap_result$hm)
dev.off()

