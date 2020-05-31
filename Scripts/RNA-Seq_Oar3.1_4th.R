####################  Following intact instruction of DER Finder  ####################

##########  DER Finder instruction  ##########

##### Step1- Loading vital packages #####


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
#library("TxDb.Hsapiens.UCSC.hg19.knownGene")
#library("org.Hs.eg.db")


##### Working with non-human data #####


options(species = "Ovis_aries")

options(chrsStyle = "Ensembl")

extendedMapSeqlevels(seqnames = c(1:26, "X"), style = "Ensembl", currentStyle = "Ensembl")                               

OarSeqinfo <- Seqinfo(seqnames = c(1:26, "X"), seqlengths = c( 275612895, 248993846, 224283230, 119255633, 107901688, 117031472, 
                                                               100079507, 90695168, 94726778, 86447213, 62248096, 79100223, 83079144, 62722625, 80923592,
                                                               71719816, 72286588, 68604602, 60464314, 51176841, 50073674, 50832532, 62330649, 42034648,
                                                               45367442, 44077779, 135437088),isCircular = rep(FALSE, 27), genome = "Oar3.1")

##### Save path of big dataframes

SavedData <- file.path("/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_4th/SavedData/")


##### Step2- Loading the Data

## Step2-1 Phenotype data

# Get pheno table
Fecuntrt <- c( Fecun = c(rep("LF",3), "HF",rep("LF",2), rep("HF",2),"LF"))
Breedtrt <- c(FecunBreed = c(rep("San",3), rep("Shall",6)))
pheno <- as.data.frame( cbind(Breed = Breedtrt, Fecun = Fecuntrt))#[-c(4,5,9),]
pheno <- map_dfc(colnames(pheno), ~ as.factor(pheno2[,.]))
colnames(pheno) <- c("Breedtrt","Fecuntrt")
rownames(pheno) <- NULL
pheno <- pheno[-c(4,5,9),]


##  Step2-2 Loading Data

bamfilespath <- list.files(path = "/home/MPourjam/RNA-Seq/DATA/Analyze/Out/Sheep/", pattern = ".Aligned.out.bam", all.files = TRUE, full.names = TRUE, recursive = TRUE )
readsnames <- gsub(".*/(.*)_Aligned.out.bam", "\\1", bamfilespath)
names(bamfilespath) <- readsnames
Coordbamfilespath <- list.files(path = "/home/MPourjam/RNA-Seq/DATA/Analyze/Out/Sheep", pattern = "sortedByCoord.out.bam$", full.names = TRUE, recursive = TRUE)
Coordbamfilespathindex <- paste0(Coordbamfilespath, ".bai")
Bamfiles <- map(seq_along(Coordbamfilespath)[-c(4,5,9)], function(i) {BamFile(file = Coordbamfilespath[i], index = Coordbamfilespathindex[i])})
bamfileslist <-  BamFileList(Bamfiles)
names(bamfileslist) <- readsnames[-c(4,5,9)]

# Creating Metadata columns and TxDb Object
met <- data.frame( IllSeqSheep = c(paste("LF_San_", rep(1:3,1), sep = "" ), c("LF_Shall_2","HF_Shall_2","HF_Shall_3")))
metadata(bamfileslist) <- met

Mart <- useMart(biomart = "ensembl_mart_94", path = "/biomart/martservice", port = 80 ,archive = FALSE)

TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1 <- makeTxDbFromBiomart(taxonomyId = 9940, port = 80, biomart = "ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org", dataset = "oaries_gene_ensembl", filter = NULL, circ_seqs = "MT", id_prefix="ensembl_",miRBaseBuild=NA)
#load(paste0(SavedData, "TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1.RData", sep = ""))
#TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1 <- makeTxDbFromGFF(file = file.path("/home/MPourjam/RNA-Seq/DATA/Reference/Ovis Areis/Ensembl/Ovis_aries.Oar_v3.1.94.gff3")
#                                                                  ,format = "gff3"
#                                                                  ,organism = "Ovis aries"
#                                                                  ,taxonomyId = 9940
#                                                                 #,chrominfo = OarSeqinfo
#                                                                  )
                
keepSeqlevels(TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1,c(1:26, "X"))
save(TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1, file = paste0(SavedData, ls()[which(str_detect(ls(), "^.*TxDb.*aries.*$"))], ".RData", sep = ""))
genomicState <- makeGenomicState(txdb =  TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1, chrs = names(seqlengths(bamfileslist)[1:27]), style = getOption("chrsStyle", "Ensembl"), currentStyle = "Ensembl")
seqinfo(genomicState) <- OarSeqinfo
saveRDS(genomicState, file = paste(SavedData, "genomicState.RDS",sep = ""))


# TxDb package building
#makeTxDbPackageFromBiomart(version = "3.5", maintainer = "<mohsenpm50@gmail.com>", author = "Mohsen Pourjam", destDir = "/home/MPourjam/RNA-Seq/RProject",
#                           biomart = "ENSEMBL_MART_ENSEMBL", host = "www.ensembl.org", dataset = "oaries_gene_ensembl", filter = NULL, circ_seqs = "MT", id_prefix="ensembl_",miRBaseBuild=NA)
#install.packages("/home/MPourjam/RNA-Seq/RProject/TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1", repos = NULL, type = "source")
library(TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1)

# Totalmapped
TotalMapped <- as.vector(map_dbl(Coordbamfilespath[-c(4,5,9)], getTotalMapped))
names(TotalMapped) <- names(bamfileslist)

# fullCov (unfiltered)
if (any(grepl("^fullCov\\.RData$", dir(SavedData)))) {
load(paste0(SavedData, list.files(SavedData, "^fullCov\\.RData$"), sep = ""))
} else {
fullCov <- fullCoverage(bamfileslist, chrs = c(rep(1:26), "X"), totalMapped = TotalMapped, targetSize = 6e+07)
save(fullCov, file = paste0(SavedData, ls()[which(str_detect(ls(), "^fullCov$"))], ".RData", sep = ""))
}
file.create("/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_4th/NewfullCovwarnings.txt")
NewfullCovwarnings <- warnings()
write_lines(as.character(NewfullCovwarnings), "/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_4th/NewfullCovwarnings.txt")


##  Step2-3 Filter coverage

# Filter coverage
if (any(grepl("^filteredCov\\.RData$", dir(SavedData)))) {
 load(paste0(SavedData, "filteredCov.RData", sep = ""))
} else {
  filteredCov <- lapply(fullCov, filterData, cutoff = 2)
  save(filteredCov, file = paste0(SavedData, ls()[which(str_detect(ls(),"^filteredCov$"))], ".RData"))
}


##### Step3 Expressed region-level analysis

##  Step3-1 Via regionMatrix()

if (any(grepl("^regionMat\\.RData$", dir(SavedData)))) {
  load(paste0(SavedData, "regionMat.RData" , sep = ""))
} else {
  regionMat <- regionMatrix(fullCov, cutoff = 30 , L = 150, chrsStyle = "Ensembl", species = "ovis_aries", currentStyle = "Ensembl")
  save(regionMat, file = paste0(SavedData, ls()[which(str_detect(ls(), "^regionMat$"))], ".RData", sep = ""))
}



##  Step3-2 Find DERs with DESeq2

# Round matrix
counts <- lapply(names(regionMat), function(i) {round(regionMat[[i]]$coverageMatrix)})

# Round matrix and specify design
dse <- lapply(counts, DESeqDataSetFromMatrix, colData = pheno, design = ~ Breedtrt + Fecuntrt)

# Perform DE analysis
dse <- lapply(dse, DESeq, test = "LRT", fitType = "local", reduced = ~ Fecuntrt)

# Extract results
deseq <- lapply(seq_along(regionMat), function(n) {cbind(regionMat[[n]]["regions"], results(dse[[n]]))})
names(deseq) <- names(regionMat)


##  Step3-3 Find DERs with limma

# We skip due to time saving


##### Step4- Single base-level F-statistics analysis

##  Step4-1 Models

#Get some idea of the library sizes
sampleDepths <- sampleDepth(collapseFullCoverage(fullCov),0.75) # Note that you would normally use the unfiltered data from all the chromosomes in this step and not just one.

# Define models
models <- makeModels(sampleDepths, testvars = pheno$Breedtrt
                     ) # Use the same models for all your chromosomes unless you have a specific reason to use chromosome-specific models. 

##  Step4-2 Find candidate DERs

# Create a analysis directory
setwd(file.path("/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_4th"))
dir.create('analysisResults3')
originalWd <- getwd()
setwd(file.path(originalWd,'analysisResults3'))

# Perform differential expression analysis
if (any(grepl("^results\\.RData$", dir(SavedData)))){
  load(paste0(SavedData, "results.RData"))
} else {
  system.time(results <- map(seq_along(filteredCov), ~ analyzeChr(chr = names(filteredCov[.]), coverageInfo = filteredCov[[.]], models,
                                                                  groupInfo = pheno$Breedtrt, writeOutput = TRUE, cutoffFstat = 5e-02, cutoffType = "theoretical",
                                                                  nPermute = 50, seeds = 20190120 + seq_len(50), returnOutput = TRUE, txdb = TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1,
                                                                  runAnnotation = TRUE, genomicState = genomicState$fullGenome, annotationPackage = NULL)))
  names(results) <- names(regionMat)
  analyzeChrWarnings <- warnings()
  save(results, file = paste0(SavedData, ls()[which(str_detect(ls(), "^results$"))], ".RData", sep = ""))
}
#######################################################################################################################################################################
# Note that for this type of analysis you might want to try a few coverage cutoffs and/or F-statistic cutoffs. One quick way to evaluate the results is to compare the width of the regions.
summary(width(results$'2'$regions$regions))
results$`2`$regions$regions[width(results$'2'$regions$regions) > 500 & as.logical(results$`2`$regions$regions$significant),]

# Width of candidate DERs
sig <- map(names(results), ~ as.logical(results[[.]]$regions$regions$significantQval))
names(sig) <- names(results)
suppressWarnings(GenomeSigQvalRegions <- makeGRangesFromDataFrame(bind_rows(map(names(results), ~ data.frame(results[[.]]$regions$regions[ sig[[.]], ]))), keep.extra.columns = TRUE))

# Use this information to make a plot
#timed <- diff(results$timeinfo)
#timed.df <- data.frame(Seconds = as.numeric(timed), Step = factor(names(timed),
#                                                                  levels = rev(names(timed))))
#ggplot(timed.df, aes(y = Step, x = Seconds)) + geom_point()


##  Step4-3 Merge results

setwd(originalWd)

# Merge results from several chromosomes.
mergeResults(chrs=c(1:26,"X"), prefix="analysisResults3",
             genomicState = genomicState$fullGenome, 
             optionsStats = results$`1`$optionsStats)

# Files created by mergeResults()
dir('analysisResults3', pattern = '.Rdata')

# fullRegions ---> The main result from mergeResults() is in fullRegions
load(file.path('analysisResults3', 'fullRegions.Rdata'))
seqinfo(fullRegions) <- OarSeqinfo

# fullAnotatedRegions
load(file.path('analysisResults3', 'fullAnnotatedRegions.Rdata'))


##### Step5 Visually explore results

# Find overlaps between regions and summarized genomic annotation
annoRegs <- suppressWarnings(annotateRegions(fullRegions, genomicState$fullGenome))

# Get the region coverage
regionCov <- getRegionCoverage(fullCov, fullRegions)
txdb <- TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1

# Now prepare some data.frames
#setwd("/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_4th/NCBIfiles/")
#makeOrgPackageFromNCBI(version = "0.1",
#                       author = "Mohsen Pourjam <mohsenpm50@gmail.com>",
#                       maintainer = "Mohsen Pourjam <mohsenpm50@gmail.com>",
#                       outputDir = ".",
#                       tax_id = "9940",
#                       genus = "Ovis",
#                       species = "aries",rebuildCache = FALSE,
#                       NCBIFilesDir = "/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_4th/NCBIfiles/",
#                       useDeprecatedStyle = FALSE
#                       )

# Overview of the candidate DERs in the genome
#undebug(derfinderPlot::plotOverview)
#plotOverview(regions = Chr2DropedSeqlevel, annotation = results$`2`$annotation, type = "qval")
map(names(results), ~ plotOverview(regions = keepSeqlevels(fullRegions, as.character(.), pruning.mode = "tidy"), annotation = results[[.]]$annotation,
            type = 'fwer'))

# Base-levle coverage plots for the first 10 regions
map(names(results), ~ plotRegionCoverage(regions = fullRegions, regionCoverage = regionCov, 
                   groupInfo = pheno$Breedtrt, nearestAnnotation = results[[.]]$annotation, 
                   annotatedRegions = annoRegs, whichRegions=1:10, txdb = txdb, scalefac = 1, 
                   ask = FALSE))
undebug(plotCluster)
# Cluster plot for the first region
plotCluster(idx = 1, regions = results$'2'$regions$regions, annotation = results$'2'$annotation,
            coverageInfo = fullCov$X, txdb = txdb, groupInfo = pheno$Breedtrt,
            titleUse = 'qval',p.ideogram = NULL)
debug(derfinderReport)
# Derfinder Report
derfinderReport(prefix = "analysisResults3",outdir = "Report",output = "DerfinderReport",project = "SheepRNA-Seq",browse = TRUE,nBestRegions = 30
                , makeBestClusters = FALSE
                ,fullCov = fullCov,hg19 = FALSE
                ,p.ideos = NULL)
               

































































































































