library('derfinder')
library('derfinderData')
library('GenomicRanges')
library('EDASeq')
library("DESeq2")
library("SummarizedExperiment")
library("Biobase")
library("DelayedArray")
library("matrixStats")
library("BiocParallel")
library("bumphunter")
library("GenomeInfoDbData")
library("AnnotationDbi")
library("GenomicFeatures")
library("RSQLite")
library("rtracklayer")
library("dplyr")
library("purrr")
library("ggplot2")
library("derfinderPlot")
library("tidyr")
library("GenomicAlignments")
library("RMySQL")
library("DBI")
library("devtools")
library("AnnotationHub")
library("biomaRt")
library('xml2')
library('curl')
library("regionReport")
library("biovizBase")
library("tidyverse")
library("knitr")
library("limma")
library("ensembldb")
library("AnnotationHub")
library("GOstats")
library("AnnotationForge")
library("org.Oaries.eg.db")
library("ggbio")
library("RMariaDB")
library("rmarkdown")
library("pryr")
library("furrr")
library("rstudioapi")

plan("multisession")
options(future.globals.maxSize= 8912896000)
options(future.globals.onMissing = "ignore")

##### Creating all RegionMats for each MCC
fullcov_path <-file.path("/")
fullcov <- load(fullcov_path)

### Creating the RegionsMat's folder

if (!"RegionMats" %in% list.dirs(SavePath, full.names = FALSE)) {
  dir.create(path = paste0(SavePath, "RegionMats"))
  RegionMatsPath <<- file.path(paste0(SavePath, "RegionMats/"))
} else {
  RegionMatsPath <<- file.path(paste0(SavePath, "RegionMats/"))
}

### Creating the RegionMats

for (c in seq_along(MCC)){   ##### Create seperate dirs for each RegionMat_MCC to store RegionMat_MCC.RData and MRG derived regions.rds 
  if (!any(str_detect(pattern = paste0(paste("RegionMats", MCC[c], sep = "_"),".RData"),
    string = list.files(RegionMatsPath, full.names = FALSE, recursive=TRUE)))) {

    file.create(file.path(paste0(paste0(paste0(RegionMatsPath,
     paste("RegionMats", MCC[c], sep = "_"),"/"), paste("RegionMats", MCC[c], sep = "_")), ".RData")))

    New_RegionMat <- regionMatrix(fullCov, cutoff = MCC[c], runfilter = FALSE ,L = c(rep(150,3), rep(100,2),rep(150, 4)),
                        chrsStyle = "Ensembl", species = "ovis_aries", currentStyle = "Ensembl")
    
    save(New_RegionMat,
         file = file.path(paste0(paste0(paste0(RegionMatsPath,
     paste("RegionMats", MCC[c], sep = "_"),"/"), paste("RegionMats", MCC[c], sep = "_")), ".RData")))
    rm(New_RegionMat)
    gc()
    print(paste0("RegionMats data of cutoff value: ", MCC[c], " and mean filter policy has been saved to ",
      file.path(paste0(paste0(RegionMatsPath, paste("RegionMats", MCC[c], sep = "_")), ".RData"))))
  }
}
rm(c)
