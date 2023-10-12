#!/usr/bin/env Rscript

#####################################################################
# 01_02_create_sampleannot.R
# created on 2022-08-24 by Irem Gunduz
# Sample annotation for snmC-seq2 data, including sample IDs, file paths, and conditions
#####################################################################

library(dplyr)

filePathCol <- "allC_FilePathfull"
sampleIdCol <- "Common_Minimal_Informative_ID"
indexed <- F
if (indexed) {
  extend <- ".tbi.tsv.gz"
} else {
  extend <- ".tsv.gz"
}
sampleAnnot <- data.table::fread("/icbb/projects/igunduz/FinalFinal_Rvsd_Consol_All_allCFileNames_allCMetadata_20221026.tsv")
sampleAnnot$Common_Minimal_Informative_ID <-	gsub("(.+_[0-9]+)_2$", "\\1",
					gsub("_Rep_?[0-9]", "", sampleAnnot$Common_Minimal_Informative_ID)
			)
sidsu <- sort(sampleAnnot$Common_Minimal_Informative_ID )

#add addittional info
s2 <- data.table::fread("/icbb/projects/igunduz/DARPA/samples_all_20221026.tsv")
sids_annot <- sort(unique(s2[["CommonMinID"]]))

colnames(s2) <-c("Salk_ID",colnames(s2)[2:19])
idx <- grepl("^Dx_Ctrl", s2[["Salk_ID"]])
missingSampleMap <- s2[idx == TRUE,"CommonMinID"]$CommonMinID
names(missingSampleMap) <- s2[idx==TRUE,"Salk_ID"]$Salk_ID
idx <- sidsu %in% names(missingSampleMap)
sidsu[idx == TRUE] <- missingSampleMap[sidsu[idx == TRUE]]
sampleAnnot$Common_Minimal_Informative_ID  <- sidsu
#idx <- sidsu[sidsu %in% names(missingSampleMap)== TRUE]
#sampleAnnot[idx] <- missingSampleMap[sidsu[idx]]
#sampleAnnot[(sampleAnnot$Common_Minimal_Informative_ID %in% idx) == TRUE, ] <- missingSampleMap[idx]

s2 <- s2 %>%
    dplyr::mutate(condition = paste0(exposure_type,"_",exposure_group)) %>%
    dplyr::rename(Common_Minimal_Informative_ID = CommonMinID) %>%
    #dplyr::filter(exposure_type != c("BA","MRSA-MSSA")) %>%
    dplyr::select(condition,age,Common_Minimal_Informative_ID) 

sampleAnnot <-merge(sampleAnnot,s2,"Common_Minimal_Informative_ID")
sampleAnnot$allC_FilePath <- gsub("allc_CGN","allcFiles",sampleAnnot$allC_FilePath)
#sampleAnnot$InfoID <- sub("_Rep2","",sampleAnnot$Common_Minimal_Informative_ID)
#sampleAnnot$InfoID <- sub("_Rep1","",sampleAnnot$InfoID)
#sampleAnnot$InfoID <- sub("_Rep3","",sampleAnnot$InfoID)

sampleAnnot$allC_FilePathfull <- paste0("/icbb/projects/igunduz/DARPA/", sampleAnnot$allC_FilePath)
sampleAnnot <- sampleAnnot[!grepl("MRSA_", sampleAnnot$Salk_ID)]
sampleAnnot <- sampleAnnot[!grepl("BA", sampleAnnot$Salk_ID)]
sampleAnnot <- sampleAnnot[!grepl(
  "/icbb/projects/igunduz/DARPA/allcFiles/Ctrl/Ctrl_3_2-Ctrl_4_2-A8-AD006.CGN-Merge.allc.tsv.gz",
  sampleAnnot$allC_FilePathfull
)]

inputFilenames <- as.character(sampleAnnot[[filePathCol]])
sampleIds <- as.character(sampleAnnot[[sampleIdCol]])
if (indexed) {
  extend <- ".tbi"
} else {
  extend <- ".tsv"
}
# make sure we have the final paths
inputFilenames <- sampleAnnot[grepl(extend, inputFilenames), ] %>%
  dplyr::select(dplyr::all_of(filePathCol))
inputFilenames <- as.character(sampleAnnot[[filePathCol]])
names(inputFilenames) <- sampleIds

missingSamples <- sampleIds[!file.exists(inputFilenames)]
write.csv(inputFilenames[missingSamples],"/icbb/projects/igunduz/DARPA/meth_missing_samples.csv")
sampleAnnot <- sampleAnnot[!sampleAnnot[[sampleIdCol]]%in% missingSamples]

inputFilenames <- as.character(sampleAnnot[[filePathCol]])
sampleIds <- as.character(sampleAnnot[[sampleIdCol]])
inputFilenames <- sampleAnnot[grepl(extend, inputFilenames), ] %>%
  dplyr::select(dplyr::all_of(filePathCol))
inputFilenames <- as.character(sampleAnnot[[filePathCol]])
names(inputFilenames) <- sampleIds

sampleAnnot <- data.table::fread("/icbb/projects/igunduz/DARPA/allc_sample_annot_final.csv")%>%
  dplyr::select(!V1)%>%
  dplyr::filter(!cell_type %in% c("Th-Eff","Tc-Eff"))
sampleAnnot$condition <- ifelse(sampleAnnot$condition == "HIV_Cro","HIV_chronic",sampleAnnot$condition)
sampleAnnot$condition <- ifelse(sampleAnnot$condition == "HIV_Pre","HIV_pre",sampleAnnot$condition)
sampleAnnot$condition <- ifelse(sampleAnnot$condition == "HIV_Acu","HIV_acute",sampleAnnot$condition)

write.csv(sampleAnnot,"/icbb/projects/igunduz/DARPA/allc_sample_annot_final.csv")
#####################################################################
