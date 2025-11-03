# GEO Export Script for ArchR scATAC-seq Data
# Author: [Your Name]
# Date: [Today]
# Description:
#   Exports PeakMatrix and GeneScoreMatrix in GEO-acceptable formats
#   using sparse Matrix Market format (.mtx.gz) + features.tsv.gz + barcodes.tsv.gz.
suppressPackageStartupMessages({
library(ArchR)
library(Matrix)
library(R.utils)
library(rtracklayer)
})

# ====== Load ArchR Project ======
project <- ArchR::loadArchRProject(
  "/icbb/projects/igunduz/archr_projects/icbb/projects/igunduz/archr_project_011023/",
  showLogo = FALSE
)
addArchRThreads(threads = 30)

# ====== Base output path ======
path <- "/icbb/projects/igunduz/ECHO_GEO_submission"
dir.create(path, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Get the fragments from the project
# -------------------------------------------------------------------------
frags <- getFragmentsFromProject(project)
out_dir <- "/icbb/projects/igunduz/ECHO_GEO_submission/fragments"
dir.create(path, recursive = TRUE, showWarnings = FALSE)

extract_barcode <- function(mc) {
  # prefer explicit barcode fields; else parse from RG like "sample#BARCODE-1"
  if ("barcode" %in% colnames(mc))       return(as.character(mc$barcode))
  if ("cell" %in% colnames(mc))          return(as.character(mc$cell))
  if ("CB" %in% colnames(mc))            return(as.character(mc$CB))
  if ("RG" %in% colnames(mc))            return(sub(".*#", "", as.character(mc$RG)))
  return(NULL)
}

for (i in seq_along(frags)) {
  gr <- frags[[i]]
  nm <- names(frags)[i]
  if (is.null(nm) || nm == "") nm <- sprintf("sample_%02d", i)

  # pull barcodes; fallback to sample name if truly missing
  bc <- extract_barcode(mcols(gr))
  if (is.null(bc)) bc <- rep(nm, length(gr))

  # make sure starts are 0-based for BED
  start(gr) <- start(gr) - 1L
  strand(gr) <- "*"

  # keep only the two BED fields we care about: name (=barcode) and score (=count)
  mcols(gr) <- DataFrame(name = as.character(bc), score = as.integer(1))

  # write gzipped BED; rtracklayer compresses if suffix is .gz
  out <- file.path(out_dir, paste0(nm, ".bed.gz"))
  export.bed(gr, con = out, format = "bed")
}

message("Exported ", length(frags), " BED files to: ", out_dir)

# -------------------------------------------------------------------------
# 1. Export PeakMatrix
# -------------------------------------------------------------------------

peak_dir <- file.path(path, "GEO_PeakMatrix")
if(!dir.exists(peak_dir)){
  dir.create(peak_dir, showWarnings = FALSE)
peakMat <- getMatrixFromProject(project, useMatrix = "PeakMatrix")

# Sparse matrix export
writeMM(assay(peakMat), file.path(peak_dir, "matrix.mtx"))
write.table(
  rownames(peakMat),
  file = file.path(peak_dir, "peaks.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)
write.table(
  colnames(peakMat),
  file = file.path(peak_dir, "barcodes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

# Compress
gzip(file.path(peak_dir, "matrix.mtx"), overwrite = TRUE)
gzip(file.path(peak_dir, "peaks.tsv"), overwrite = TRUE)
gzip(file.path(peak_dir, "barcodes.tsv"), overwrite = TRUE)
}
# -------------------------------------------------------------------------
# 2. Export GeneScoreMatrix (Gene Activity)
# -------------------------------------------------------------------------

ga_dir <- file.path(path, "GEO_GeneActivity")
#if(!dir.exists(ga_dir)){
dir.create(ga_dir, showWarnings = FALSE)
geneAct <- getGroupSE(project, useMatrix = "GeneScoreMatrix",groupBy= "Sample")
matrix_data <- as(assay(geneAct), "dgCMatrix") # Convert to sparse matrix

writeMM(matrix_data, file.path(ga_dir, "matrix.mtx"))
write.table(
  rowData(geneAct),
  file = file.path(ga_dir, "genes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)
write.table(
  colnames(geneAct),
  file = file.path(ga_dir, "barcodes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

gzip(file.path(ga_dir, "matrix.mtx"), overwrite = TRUE)
gzip(file.path(ga_dir, "genes.tsv"), overwrite = TRUE)
gzip(file.path(ga_dir, "barcodes.tsv"), overwrite = TRUE)

#}

# -------------------------------------------------------------------------
# 3. Export chromVAR deviations
# -------------------------------------------------------------------------

cv_dir <- file.path(path, "GEO_chromVAR")
if(!dir.exists(cv_dir)){
  dir.create(cv_dir, showWarnings = FALSE)
chromVarDev <- getMatrixFromProject(project, useMatrix = "MotifMatrix")

# Save matrix in sparse format
writeMM(assay(chromVarDev), file.path(cv_dir, "matrix.mtx"))

# Save motifs/features (rownames are motif names)
write.table(
  rownames(chromVarDev),
  file = file.path(cv_dir, "motifs.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

# Save barcodes/cell IDs
write.table(
  colnames(chromVarDev),
  file = file.path(cv_dir, "barcodes.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
)

# Compress all
gzip(file.path(cv_dir, "matrix.mtx"), overwrite = TRUE)
gzip(file.path(cv_dir, "motifs.tsv"), overwrite = TRUE)
gzip(file.path(cv_dir, "barcodes.tsv"), overwrite = TRUE)

}
message("✅ GEO export complete! Files are saved in: ", path)
