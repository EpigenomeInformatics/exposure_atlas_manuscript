
prettyOrderMat <- function(mat, scale=TRUE, cutOff=1, lmat=NULL, clusterCols=TRUE){
  # Reorder mat in a prettier way for plotting
  # Adapted from Jeff's ArchR .binarySort
  ###################################
  # mat = matrix (like) object to sort
  # scale = should mat be scaled before building logical mat
  # cutOff = cutoff for lmat
  # lmat = logical matrix for ordering rows (binary sorting)
  # clusterCols = should columns be clustered?
  mat <- as.matrix(mat)

  if(is.null(lmat)){
    # Compute row Z-scores
    if(scale){
      lmat <- sweep(mat - rowMeans(mat), 1, matrixStats::rowSds(mat), `/`)
    }else{
      lmat <- mat
    }
    # Logical matrix of values passing cutoff 
    lmat <- lmat >= cutOff
  }

  # Transpose:
  mat <- t(mat)
  lmat <- t(lmat)

  # Identify column ordering:
  if(clusterCols){
    hc <- hclust(dist(mat))
    colIdx <- hc$order
    mat <- t(mat[colIdx,])
    lmat <- t(lmat[colIdx,])
  }else{
    mat <- t(mat)
    lmat <- t(lmat)
    hc <- NULL
  }

  # Identify row ordering:
  rowIdx <- do.call("order", c(as.data.frame(lmat)[seq_len(ncol(lmat))], list(decreasing = TRUE)))
  mat <- mat[rowIdx,]

  return(list(mat=mat, hclust=hc))
}


# Heatmap wrapper:
BORHeatmap <- function(
  mat, # Data to plot (matrix or dataframe)
  limits = NULL, # Enforced limits for colormap (2 dimensional array)
  clusterCols = TRUE, # Should columns be clustered
  clusterRows = TRUE, # Should rows be clustered
  labelCols = FALSE, # Should columns be labeled
  labelRows = FALSE, # Should rows be labeled
  dataColors = NULL, # Colormap for plotting data
  dataColorMidPoint = NULL, # The data value to be the middle of the color map
  customRowLabel = NULL,
  customRowLabelIDs = NULL,
  customColLabel = NULL,
  customColLabelIDs = NULL,
  customLabelWidth = 0.15,
  useRaster = TRUE, # Should heatmap be rasterized
  rasterDevice = "CairoPNG",
  rasterQuality = 5, # Raster quality. Higher is {better?}
  fontSize = 6, # Font size for labels
  showColDendrogram = FALSE, # Should the column dendrogram be shown
  showRowDendrogram = FALSE, # Should the row dendrogram be shown
  borderColor = NA, # Color for lines between cells
  mapname = " ", # 'Name' to give heatmap
  legendTitle = " ", # Name of legend
  ...
){
  
  #Packages
  suppressPackageStartupMessages(require(ComplexHeatmap))
  suppressPackageStartupMessages(require(circlize))
  
  # Make sure mat is actually a matrix
  if(!is.matrix(mat)){
    message("'mat' needs to be a matrix. Converting...")
    mat <- as.matrix(mat)
  }
  
  # Prepare color function
  if(!is.null(limits)){
    ll <- limits[1]
    ul <- limits[2]
  }else{
    ll <- min(mat, na.rm=TRUE)
    ul <- max(mat, na.rm=TRUE)
  }
  # If no colormap provided, use solarExtra
  if(is.null(dataColors)){
    dataColors <- c("1"='#3361A5', "2"='#248AF3', "3"='#14B3FF', 
                    "4"='#88CEEF', "5"='#C1D5DC', "6"='#EAD397', 
                    "7"='#FDB31A', "8"= '#E42A2A', "9"='#A31D1D')
  }
  dataColFun <- makeColFun(ll, ul, dataColors, midpoint = dataColorMidPoint)
  
  message("Preparing Heatmap...")
  hm <- Heatmap(
    # Main components:
    matrix = mat,
    name = mapname,
    col = dataColFun,
    
    # Legend options:
    heatmap_legend_param = list(
      color_bar = "continuous",
      legend_direction = "vertical",
      legend_width = unit(1, "cm"),
      title = legendTitle
    ),
    rect_gp = gpar(col = borderColor), 
    
    # Column options:
    show_column_names = labelCols,
    cluster_columns = clusterCols,
    show_column_dend = showColDendrogram,
    clustering_method_columns = "ward.D2",
    #column_names_gp = gpar(fontsize = fontSize), 
    
    # Row options:
    show_row_names = labelRows,
    cluster_rows = clusterRows,
    show_row_dend = showRowDendrogram,
    clustering_method_rows = "ward.D2",
    #row_names_gp = gpar(fontsize = fontSize), 
    
    # Raster info:
    use_raster = useRaster,
    raster_device = rasterDevice,
    raster_quality = rasterQuality,

    # Other
    ...
  )
  
  # Add row labels if provided:
  if(!is.null(customRowLabel)){
    if(is.null(customRowLabelIDs)){
      customRowLabelIDs <- rownames(mat)[customRowLabel]
    }
    hm <- hm + rowAnnotation(
      link = anno_mark(at = customRowLabel, labels = customRowLabelIDs, labels_gp = gpar(fontsize = fontSize)),
      width = unit(customLabelWidth, "cm") + max_text_width(customRowLabelIDs)
    )
  }
  
  return(hm)
}

# This is used primarily for making colormaps for ComplexHeatmap
makeColFun <- function(start, end, cmap, midpoint = NULL){
  # Make a color ramp function from provided start and end breaks,
  # and optionally a midpoint
  cmapLen <- length(cmap)
  if(!is.null(midpoint)){
    interpolate <- function(c1, c2, colorspace = "Lab"){
      rgb2hex(colorRamp(c(c1, c2), space = colorspace)(0.5))
    }
    if(length(cmap) %% 2 == 0){
      # Interpolate middle colors if necessary to get midpoint
      preMidIdx <- floor(cmapLen / 2)
      midCol <- interpolate(cmap[preMidIdx], cmap[preMidIdx + 1])
      cmap <- c(cmap[1:preMidIdx], midCol, cmap[(preMidIdx + 1):cmapLen])
      cmapLen <- length(cmap)
    }
    midIdx <- ceiling(cmapLen / 2)
    breaks <- c(seq(start, midpoint, length.out = midIdx), seq(midpoint, end, length.out = midIdx)[2:midIdx])
  } else {
    breaks <- seq(start, end, length.out = cmapLen)
  }
  colorRamp2(breaks, cmap)
}


# load libraries
suppressPackageStartupMessages({
  library(ArchR)
  library(dplyr)
  library(mclust)
  library(ComplexHeatmap)
  library(ggplot2)
  library(rescueR)
})
project <- ArchR::loadArchRProject(outputDir, showLogo = FALSE)
# construct confusion matrix
df <- as.data.frame(as.matrix(confusionMatrix(
  project$Clusters_0.8,
  project$predictedGroupAtlas
)))
cM <- computeJaccardIndex(df, heatmap = FALSE)
cM <- prettyOrderMat(t(cM),clusterCols=TRUE)$mat %>% t()

whitePurple = c('#f7fcfd','#e0ecf4','#bfd3e6','#9ebcda','#8c96c6',
                  '#8c6bb1','#88419d','#810f7c','#4d004b')

ht_opt$simple_anno_size <- unit(0.25, "cm")

#ra <- HeatmapAnnotation(atac_cluster=rownames(cM),col=list(atac_cluster=atac_label_cmap), which="row", show_legend=c("atac_cluster"=FALSE))
#ta <- HeatmapAnnotation(rna_cluster=colnames(cM),col=list(rna_cluster=rna_label_cmap), show_legend=c("rna_cluster"=FALSE))
hm <- BORHeatmap(
  cM, 
  #limits=c(0,1), 
  dataColorMidPoint = 0.4,
  #clusterCols=TRUE, clusterRows=TRUE,
  labelCols=TRUE, labelRows=TRUE,
  dataColors = whitePurple,
  #left_annotation = ra,
  #top_annotation = ta,
  showColDendrogram = F, # Should the column dendrogram be shown
  showRowDendrogram = F,
  row_names_side = "left",
  width = ncol(cM)*unit(0.5, "cm"),
  height = nrow(cM)*unit(0.5, "cm"),
  border_gp=gpar(col="black") # Add a black border to entire heatmap
  )
pdf(paste0("/icbb/projects/igunduz/jaccard_heatmap.pdf"), width=6, height=6)
draw(hm)
dev.off()
