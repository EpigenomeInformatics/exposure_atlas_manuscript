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
out_dir <- "/icbb/projects/igunduz/ECHO_GEO_submission/fragments_141125"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (i in seq_along(frags)) {
  gr <- frags[[i]]
  nm <- names(frags)[i]
  if (is.null(nm) || nm == "") nm <- sprintf("sample_%02d", i)

  # Remove ".tsv.gz" from the name if it exists
  nm <- sub("\\.tsv\\.gz$", "", nm)

  # make sure starts are 0-based for BED
  start(gr) <- start(gr) - 1L
  strand(gr) <- "*"

  # keep only RG as the 'name' field in BED
  mcols(gr) <- DataFrame(name = as.character(mcols(gr)$RG))

  # write gzipped BED; rtracklayer compresses if suffix is .gz
  out <- file.path(out_dir, paste0(nm, ".bed.gz"))
  export.bed(gr, con = out, format = "bed")
  message("Exported ", out, " as a BED file")
}

message("Exported ", length(frags), " BED files to: ", out_dir)
