# Load necessary libraries
suppressPackageStartupMessages(library(JASPAR2018))
suppressPackageStartupMessages(library(TFBSTools))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(motifmatchr))
suppressPackageStartupMessages(library(ChrAccR))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(BSgenome))
suppressPackageStartupMessages(library(muLogR))

# Define the directory to save motif logos
motif_dir <- "/icbb/projects/igunduz/finalize_echo_050824/figures/motif_logos/"
if (!dir.exists(motif_dir)) {
  dir.create(motif_dir)
}
#source("/icbb/projects/igunduz/irem_github/ChrAccR/R/utils_motifs.R")
hmSeqLogo <- function(pwm, x=unit(0.5, "npc"), y=unit(0.5, "npc"), width=1, height=1, ic.scale=TRUE){
	if (!requireNamespace("grid")) logger.error(c("Could not load dependency: grid"))
	# # convert units to numbers
	# unitType <- attr(x, "unit")
	unitType <- "npc"
	x <- as.numeric(x)
	y <- as.numeric(y)
	width <- as.numeric(width)
	height <- as.numeric(height)


	# convert the PWM to matrix
	if (class(pwm) == "pwm") {
		pwm <- pwm@pwm
	} else if (class(pwm) == "PWMatrix") {
		pwm <- ChrAccR:::PWMatrixToProbMatrix(pwm)
	} else if (class(pwm) == "data.frame") {
		pwm <- as.matrix(pwm)
	} else if (class(pwm) != "matrix"){
		stop("pwm must be of class matrix or data.frame")
	}
	if (any(abs(1 - apply(pwm,2,sum)) > 0.01)) stop("Columns of PWM must add up to 1.0")

	chars <- c("A","C","G","T")
	letters <- list(x=NULL,y=NULL,id=NULL,fill=NULL)
	npos <- ncol(pwm)

	if (ic.scale) {
		facs <- seqLogo:::pwm2ic(pwm)
		facs <- facs/max(facs) # scale columns to max information content
	} else {
		facs <- rep(1, npos)
	}

	wt <- width / npos
	x.pos <- x - width/2
	for (j in 1:npos) {
		column <- pwm[,j]
		hts <- 0.99*column*facs[j]*height
		letterOrder <- order(hts)

		y.pos <- y-height/2
		for (i in 1:length(chars)) {
			letter <- chars[letterOrder[i]]
			ht <- hts[letterOrder[i]]
			if (ht>0) letters <- seqLogo:::addLetter(letters, letter, x.pos, y.pos, ht, wt, fill=c("A"="green","C"="red","G"="blue","T"="orange"))
			y.pos <- y.pos + ht #+ 0.01
		}
		x.pos <- x.pos + wt
	}
	# print(str(letters))
	grid::grid.polygon(x=unit(letters$x, unitType), y=unit(letters$y, unitType), id=letters$id, gp=grid::gpar(fill=letters$fill,col="transparent"))
}

# Function to extract and save motif logos for specific transcription factors
save_motif_logos <- function(tfNames, motifDb = "jaspar2020", motif_dir) {
  
  # Fetch the motif data from the specified database
  motifObj <- prepareMotifmatchr("hg38", motifDb)$motifs
  
  for (tf in tfNames) {
    # Use a regex pattern to match the exact transcription factor name
    full_mn <- grep(tf, names(motifObj), value = TRUE, ignore.case = TRUE)
    
    if (length(full_mn) > 0) {
      # Open a PDF device
      pdf_file <- paste0(motif_dir, tf, "_motif_logo.pdf")
      pdf(pdf_file, width = 6, height = 2)
      
      for (mn in full_mn) {
        pwm <- motifObj[[mn]]  # Extract the PWM for the current match

        # Generate the motif logo using hmSeqLogo and save it
        grid.newpage()
        hmSeqLogo(pwm)
        
        # Add a title with the original motif name
        grid.text(mn, y = unit(1, "npc") - unit(1, "lines"), just = "center", gp = gpar(fontsize = 10))
      }

      # Close the PDF device
      dev.off()
      
      cat("Saved motif logos for", tf, "as", pdf_file, "\n")
    } else {
      cat("Motif for transcription factor", tf, "not found in the database.\n")
    }
  }
}

# Example usage
tfNames <- c("NRF1")#"RELA", "BATF", "RUNX3")
#tfNames <- c("CREB1", "CEBPA", "ATF4", "CEBPD", "REL", "BATF::JUN", "SPIC", "IRF2") 
save_motif_logos(tfNames, motifDb = "jaspar2018", motif_dir = motif_dir)

