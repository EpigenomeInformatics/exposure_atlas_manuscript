suppressPackageStartupMessages({
library(ArchR)
library(dplyr)
library(GenomicRanges)
library(SummarizedExperiment)
library(Seurat)
library(ggplot2)
library(rescueR)
})
set.seed(12) # set seed

addArchRThreads(threads = 30) # set the cores
addArchRGenome("hg38") # set the reference genom
# load the functions
source("/icbb/projects/igunduz/sc_epigenome_pathogen_exposure/utils/jaccard.R")
outputDir <- "/icbb/projects/igunduz/Tcell_subset/"
project <- ArchR::loadArchRProject(outputDir,showLogo=FALSE)


#wget https://changseq.s3.amazonaws.com/Jeff/10x_ScATAC/scATAC_TME_TCells_SummarizedExperiment.final.rds
tcell_se <- readRDS("/icbb/projects/igunduz/DARPA/scATAC_TME_TCells_SummarizedExperiment.final.rds")
tcell_se_clean <- tcell_se[complete.cases(rowData(tcell_se)$SYMBOL), ]
gene_names <- rowData(tcell_se_clean)$SYMBOL
rownames(tcell_se_clean) <- gene_names
tcell_matrix <- assay(tcell_se_clean)
seurat_obj <- CreateSeuratObject(counts = tcell_matrix)

remapClust <- c(
  "Cluster1" = "CD4T_naive",
  "Cluster2" = "Activated_CD4T",
  "Cluster3" = "Th1",
  "Cluster4" = "CD4T_mem",
  "Cluster5" = "Th17",
  "Cluster6" = "Tfh1",
  "Cluster7" = "Tfh2",
  "Cluster8" = "Treg1",
  "Cluster9" = "Treg2",
  "Cluster10" = "Treg3",
  "Cluster11" = "Treg4",
  "Cluster12" = "Effector_CD8T",
  "Cluster13" = "CD8T_naive",
  "Cluster14" = "CD8T_mem",
  "Cluster15" = "Early_TEx",
  "Cluster16" = "Intermediate_TEx",
  "Cluster17" = "Terminal_TEx",
  "Cluster18" = "Other_T",
  "Cluster19" = "Other_T"
)

# Add colData as metadata to the Seurat object
metadata <- colData(tcell_se_clean)

# Replace T_Cell_Cluster with new names
metadata$T_Cell_Cluster <- remapClust[metadata$T_Cell_Cluster]

# Convert metadata columns to vectors
metadata <- as.data.frame(metadata, stringsAsFactors = FALSE) # Ensure metadata is a data frame

# Add metadata to the Seurat object
seurat_obj <- AddMetaData(object = seurat_obj, metadata = metadata)


# add iterative LSI
project <- addIterativeLSI(
  ArchRProj = project,
  useMatrix = "TileMatrix",
  name = "IterativeLSI_tcell_final",
  iterations = 2,
  clusterParams = list(
    resolution = c(0.2),
    sampleCells = 10000,
    n.start = 10
  ),
  varFeatures = 25000,
  dimsToUse = 1:30,
  saveIterations = TRUE,
  nPlot = 100000,
  force=TRUE
)
# perform harmony
project <- addHarmony(
  ArchRProj = project,
  reducedDims = "IterativeLSI_tcell_final",
  name = "Harmony_TCELL_final",
  force = TRUE,
  groupBy = "sample_exposure_type")

# add clusters
project <- addClusters(
  input = project,
  reducedDims = "Harmony_TCELL_final",
  method = "Seurat",
  name = "tcellClust_final",
  resolution = 0.8,
  force = TRUE,
  maxClusters = 25
)

# add umap
project <- addUMAP(
  ArchRProj = project,
  reducedDims = "Harmony_TCELL_final",
  name = "UMAPTcell_final",
  nNeighbors = 40,
  minDist = 0.5,
    force = TRUE,
  metric = "cosine",
)
 
# integration
project <- addGeneIntegrationMatrix(
    ArchRProj = project, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "tcellGeneIntegrationMatrix",
    reducedDims = "IterativeLSI_tcell_final",
    seRNA = seurat_obj,
    addToArrow = FALSE,
    force = TRUE,
    groupRNA = "T_Cell_Cluster",
    nameCell = "predicted_tcell_cluster_final",
    nameGroup = "predictedGroup_tcell_final",
    nameScore = "predictedScore_tcell_final"
)


# Reassign the cluster names
cM <- confusionMatrix(project$tcellClust_final, project$predictedGroup_tcell_final)
jh <- Jaccard(confusionMatrix(project$predictedGroup_tcell, project$tcellClust)) # Jaccard heatmap
pdf(paste0("/icbb/projects/igunduz/Tcell_subset/Plots/jaccard_heatmap.pdf"), width = 10, height = 10)
jh
dev.off()

labelOld <- rownames(cM)
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
#labelNew2 <- mapLabels(labelNew, oldLabels = names(remapClust), newLabels = remapClust)
project$TcellSub <- mapLabels(project$tcellClust_final, newLabels = labelNew, oldLabels = labelOld)
table(project$TcellSub)


cell <- as.data.frame(project@cellColData)%>%
  dplyr::filter(!sample_exposure_group %in% c("BA_na","BA_vac")) 

cell <- table(cell$TcellSub,cell$sample_exposure_group)
cell <- as.matrix(cell)
cell <- scale(cell, center=FALSE, scale=colSums(cell))
cell <- reshape2::melt(cell)
colnames(cell) <- c("CellTypes","Exposure","FractionOfCells")
cell$Exposure <-forcats::fct_reorder(cell$Exposure, cell$FractionOfCells,mean, .desc = T)
#define the order
cell$Exposure <- factor(cell$Exposure, levels=c("C19_ctrl","C19_mild","C19_mod","C19_sev",
                            "HIV_ctrl","HIV_acu","HIV_chr","Influenza_ctrl",
                            "Influenza_d3","Influenza_d6","Influenza_d30",
                            "OP_low","OP_med","OP_high"))

#custom_palette <- c("#08519c", "#3182bd", "#6baed6", "#1b7837", "#5aae61", "#a6dba0", "#762a83", "#c2a5cf")
stacked <- ggplot(cell, aes(fill=CellTypes, y=FractionOfCells, x=Exposure)) + 
    geom_bar(position="fill", stat="identity")+
    scale_fill_manual(values=colpal.bgp2 )+
      labs(x = "Exposure Type",
       y = "Cell Fractions")+ 
     theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
if(!dir.exists("/icbb/projects/igunduz/Tcell_subset/Plots/")){
  dir.create("/icbb/projects/igunduz/Tcell_subset/Plots/")
}
pdf(paste0("/icbb/projects/igunduz/Tcell_subset/Plots/cellFraction_STACKED_plot2.pdf"), width = 15, height = 10)
stacked
dev.off()


df <- getEmbedding(project, embedding = "UMAPTcell_final", returnDF = TRUE)
colnames(df) <- c("UMAP1", "UMAP2")
df$Group <- project$TcellSub


# Create UMAP plot with labels for each cluster
umap_plot <- ggplot(df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point() +
  #geom_text(data = summary_df, aes(x = UMAP1, y = UMAP2, label = Group), 
  #           color = "black", size = 8, fontface = "bold", show.legend = FALSE) +
  scale_color_manual(values=colpal.bgp2)+
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme_classic() +
  theme(legend.position = "bottom")
ggsave(plot=umap_plot,"/icbb/projects/igunduz/Tcell_subset/Plots/umap2.pdf", width = 10, height = 10)


# save the project
saveArchRProject(project, outputDirectory = outputDir, load = FALSE)
#table(predictedGroup_tcell_updated)
