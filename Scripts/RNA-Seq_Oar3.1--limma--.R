# TEMPORARY
####################  Following intact instruction of DER Finder  ####################

##########  DER Finder instruction - Testing using various cutoff and TestVariable arguments ##########

##### Step1-1 Loading vital packages #####


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

ggplot(as.data.frame(fullFstats$`1`), aes(x = c(1:nrow(as.data.frame(fullFstats$`1`))), y = value)) + geom_point()
# Working with non-human data 

options(species = "Ovis_aries")

options(chrsStyle = "Ensembl")

extendedMapSeqlevels(seqnames = c(1:26, "X"), style = "Ensembl", currentStyle = "Ensembl")                               

OarSeqinfo <- Seqinfo(seqnames = c(1:26, "X"), seqlengths = c( 275612895, 248993846, 224283230, 119255633, 107901688, 117031472, 
                                                               100079507, 90695168, 94726778, 86447213, 62248096, 79100223, 83079144, 62722625, 80923592,
                                                               71719816, 72286588, 68604602, 60464314, 51176841, 50073674, 50832532, 62330649, 42034648,
                                                               45367442, 44077779, 135437088),isCircular = rep(FALSE, 27), genome = "Oar3.1")
SavedData <- file.path("/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_V/SavedData/")


## Step1-2 Setting critical variables 

# Get pheno table
Fecuntrt <- c(c(rep("LF",3), "HF",rep("LF",2), rep("HF",2),"LF"))
Breedtrt <- c(c(rep("San",3), rep("Shall",6)))
ThreeGroup <- c(rep("San", 3), rep("Shall_OUT",2), rep("Shall",3), "Shall_OUT")
ReverseBreed <- c(rep("Shall", 5), rep("San",3), rep("Shall",1))
pheno <- as.data.frame( cbind(Breed = Breedtrt, Fecun = Fecuntrt, ThreeGroup = ThreeGroup, Reverse = ReverseBreed))
pheno <- map_dfc(colnames(pheno), ~ as.factor(pheno[,.]))
colnames(pheno) <- c("Breedtrt","Fecuntrt","ThreeGroup", "ReverseBreed")
rownames(pheno) <- NULL

# setting result determiner variables
tstVar <- colnames(pheno)
adjVar <- c(colnames(pheno)[-3], 'NULL')
CF_Type <- c("mean","one")
CF_Value <- seq(10,12 , by = 0.1)
CriticalVar <- list(tstVar = tstVar, adjVar =adjVar, CF_Type = CF_Type, CF_Value = CF_Value)

#pmap(CriticalVar, as.character(CF_Value), ~file.path(paste("/home/MPourjam/RNA-Seq/RProject","RNA-Seq_Oar3.1_5th"
#                                                               ,paste0(as.character(..1),"_",as.character(..2),"_"
#                                                                       ,as.character(..3),as.character(..4)), sep = "/"),recursive = TRUE))


dirNamePermute <- function(list) {
  stopifnot(is.list(list))
  stopifnot(all(map_lgl(seq_len(length(list)), ~ !is.list(list[[.]]))))
  # i <- seq(0,length(list)-1,+1)[-1] 
  # x <- map(seq(1,length(list)-1,1), ~.)
  if ((length(list) %% 2) == 0) {
    i <- length(list) / 2
    x <- map(seq_len(i), ~ .)
    for(n in seq_len(i)) {
      x[[n]] <- map(as.character(list[[((n*2)-1)]]), ~paste0(.,"_",as.character(list[[(n*2)]])))
    }
    y <- map(x, ~as.character(flatten(.)))
    yTwoByTwo <- map(seq_len(i), ~as.table(str_split_fixed(y[[.]], "_",2)))
    s <- map(yTwoByTwo, ~ .[,1] != .[,2])
    for ( x in seq_len(i)){
      y[[x]] <- y[[x]][s[[x]]]
    }
    FirstTwo_LastTwo <- map_dfc(seq_along(y[[1]]), ~paste0(y[[1]][.],"_",y[[2]])) 
    FinalComb <- str_split_fixed(unlist(flatten(FirstTwo_LastTwo)), pattern = "_", n = length(list))
    base::colnames(FinalComb) <- names(list)
    print(list(EachVarComb_df = FinalComb, JoinedComb = data.frame(FinalComb = unlist(flatten(FirstTwo_LastTwo)))))
    #print(data.frame(FinalComb = unlist(flatten(FirstTwo_LastTwo))))
    
  } else {
    i <- (length(list)-1) / 2
    x <- map(seq_len(i) , ~ .)
    for(n in seq_len(i)) {
      x[[n]] <- map(as.character(list[[((n*2)-1)]]), ~paste0(.,"_",as.character(list[[(n*2)]])))
    }
    y <- map(x, ~as.character(flatten(.)))
    yTwoByTwo <- map(seq_len(i), ~as.table(str_split_fixed(y[[.]], "_",2)))
    s <- map(yTwoByTwo, ~ .[,1] != .[,2])
    for ( x in seq_len(i)){
      y[[x]] <- y[[x]][s[[x]]]
    }
    FirstTwo_LastTwo <- map_dfc(seq_along(y[[1]]), ~paste0(y[[1]][.],"_",y[[2]]))
    FinalComb <- str_split_fixed(unlist(flatten(map_dfc(flatten(FirstTwo_LastTwo),~paste0(as.character(.), "_", as.character(list[[length(list)]]))))),pattern = "_", n = length(list))
    base::colnames(FinalComb) <- names(list)
    print(list(EachVarComb_df = FinalComb, JoinedComb = data.frame(AllComb = unlist(flatten(map_dfc(flatten(FirstTwo_LastTwo),~paste0(as.character(.), "_", as.character(list[[length(list)]]))))))))
    #print(data.frame(AllComb = unlist(flatten(map_dfc(flatten(FirstTwo_LastTwo),~paste0(as.character(.), "_", as.character(list[[length(list)]])))))))
  }
  
}

dirNamesPer <- dirNamePermute(CriticalVar)
dirNamesPer <- as.character(dirNamesPer[[2]][,1])


##### Step2 Loading Data

## Step2-1 Loading and creating Raw Data

bamfilespath <- list.files(path = "/home/MPourjam/RNA-Seq/DATA/Analyze/Out/Sheep/"
                           , pattern = ".Aligned.out.bam", all.files = TRUE, full.names = TRUE
                           , recursive = TRUE )

readsnames <- gsub(".*/(.*)_Aligned.out.bam", "\\1", bamfilespath)
names(bamfilespath) <- readsnames

Coordbamfilespath <- list.files(path = "/home/MPourjam/RNA-Seq/DATA/Analyze/Out/Sheep"
                                , pattern = "sortedByCoord.out.bam$", full.names = TRUE
                                , recursive = TRUE)

Coordbamfilespathindex <- paste0(Coordbamfilespath, ".bai")
Bamfiles <- map(seq_along(Coordbamfilespath), function(i) {BamFile(file = Coordbamfilespath[i], index = Coordbamfilespathindex[i])})
bamfileslist <-  BamFileList(Bamfiles)
names(bamfileslist) <- readsnames
pheno$Samples <- names(bamfileslist)

# Creating Metadata columns and TxDb Object
met <- data.frame( IllSeqSheep = c(paste("LF_San_", rep(1:3,1), sep = "" )
                                   , c("HF_Shall_1", "LF_Shall_1","LF_Shall_2","HF_Shall_2","HF_Shall_3", "LF_Shall_3")))
metadata(bamfileslist) <- met

Mart <- useMart(biomart = "ensembl_mart_96", path = "/biomart/martservice", port = 80 ,archive = FALSE)
ensembl = useMart("ensembl", dataset = "oaries_gene_ensembl")

txdb <- TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1 <- makeTxDbPackageFromBiomart(taxonomyId = 9940, port = 80, biomart = "ENSEMBL_MART_ENSEMBL"
                                                                      , host = "www.ensembl.org", dataset = "oaries_gene_ensembl"
                                                                      , filter = NULL, circ_seqs = "MT", id_prefix="ensembl_", miRBaseBuild = NA
                                                                      , version = "1.0"
                                                                      , maintainer = "Mohsen Pourjam <Pourjam.cs@hotmail.com>"
                                                                      , author = "Mohsen Pourjam")
library(package= "TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1",lib.loc = NULL)
keepSeqlevels(TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1::TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1,c(1:26, "X"))
#plotkeepSeqlevels(TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1,c(1:26, "X"))
save(TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1, file = paste0(SavedData, ls()[which(str_detect(ls(), "^.*TxDb.*aries.*$"))], ".RData", sep = ""))
library(TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1)

genomicState <- makeGenomicState(txdb =  TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1::TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1
                                 , chrs = names(seqlengths(bamfileslist)[1:27])
                                 , style = getOption("chrsStyle", "Ensembl"), currentStyle = "Ensembl")

seqinfo(genomicState) <- OarSeqinfo
save(genomicState, file = paste(SavedData, "genomicState.RData",sep = "")) 

# Totalmapped
TotalMapped <- as.vector(map_dbl(Coordbamfilespath, getTotalMapped))
names(TotalMapped) <- names(bamfileslist)


##  Step2-2 Full Coverage DataFrame

# fullCov (unfiltered)
if (any(grepl("^fullCov\\.RData$", dir(SavedData)))) {
  load(paste0(SavedData, list.files(SavedData, "^fullCov\\.RData$"), sep = ""))
} else {
  fullCov <- fullCoverage(bamfileslist, chrs = c(rep(1:26), "X"), totalMapped = TotalMapped, targetSize = 4e+07, returnMean = TRUE)
  save(fullCov, file = paste0(SavedData, ls()[which(str_detect(ls(), "^fullCov$"))], ".RData", sep = ""))
}
file.create(paste0(SavedData,"/", "NewfullCovwarnings.txt"))
NewfullCovwarnings <- warnings()
write_lines(as.character(NewfullCovwarnings), path = paste0(SavedData,"/", "NewfullCovwarnings.txt"))

filteredCov <-  lapply(fullCov, filterData, cutoff = 2, filter = "mean" , TotalMapped = rep(1, length(bamfileslist)), targetSize = 1)
##### Step3 Expressed region-level analysis

##  Step3-1 Via regionMatrix()

CF_Value_Type <- unique(str_split_fixed(dirNamesPer, "_",4)[,3:4])
CF_Value_Type <- CF_Value_Type[CF_Value_Type[,1] == "mean",] # Make it in case you only need "mean" filter policy # Critical as deciding wether to use 'RegionMat' or ''

if (!"RegionMats" %in% list.dirs(SavedData, full.names = FALSE)) {
  dir.create(path = paste0(SavedData, "RegionMats"))
  RegionMatsPath <<- file.path(paste0(SavedData, "RegionMats"))
} else {
  RegionMatsPath <<- file.path(paste0(SavedData, "RegionMats"))
}

for (c in seq_len(nrow(CF_Value_Type))){
  if (!any(str_detect(pattern = paste0(paste("RegionMats", CF_Value_Type[c,1], CF_Value_Type[c,2], sep = "_"),".RData"),string = list.files(RegionMatsPath, full.names = FALSE)))) {
    assign(paste("RegionMats", CF_Value_Type[c,1], CF_Value_Type[c,2], sep = "_"), 
           regionMatrix(fullCov, cutoff = as.numeric(CF_Value_Type[c,2]), runfilter = FALSE ,L = c(rep(150,3), rep(100,2),rep(150, 4)),
                        chrsStyle = "Ensembl", species = "ovis_aries", currentStyle = "Ensembl"))
    file.create(file.path(paste0(paste(RegionMatsPath, paste("RegionMats", CF_Value_Type[c,1], CF_Value_Type[c,2], sep = "_"),sep = "/"), ".RData")))
    save(list = ls()[which(str_detect(ls(), paste("RegionMats", CF_Value_Type[c,1], CF_Value_Type[c,2], sep = "_")))],
         file = file.path(paste0(paste(RegionMatsPath, paste("RegionMats", CF_Value_Type[c,1], CF_Value_Type[c,2], sep = "_"),sep = "/"), ".RData")))
    rm(list = ls()[which(str_detect(ls(), paste("RegionMats", CF_Value_Type[c,1], CF_Value_Type[c,2], sep = "_")))]) 
    print(paste0("RegionMats data of cutoff value ", CF_Value_Type[c,2], " and (", CF_Value_Type[c,1] , ") filter policy has been saved."))
  }
}
rm(c)

RegionMatFilesPath <- list.files(path = "/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_V/SavedData/RegionMats", full.names = TRUE)


##  Step3-2 Find DERs with DESeq2

#Round matrix
#counts <- lapply(names(regionMat), function(i) {round(regionMat[[i]]$coverageMatrix)})

# Round matrix and specify design
#dse <- lapply(counts, DESeqDataSetFromMatrix, colData = pheno, design = ~ Breedtrt + Fecuntrt)

# Perform DE analysis
#dse <- lapply(dse, DESeq, test = "LRT", fitType = "local", reduced = ~ Fecuntrt)
#names(dse) <- names(regionMat)

# Extract results
#deseq <- lapply(seq_along(regionMat), function(n) {cbind(regionMat[[n]]["regions"], results(dse[[n]]))})
#names(deseq) <- names(regionMat)


##  Step3-3 Find DERs with limma

# We skip due to save time


##### Step4- Single base-level F-statistics analysis

##  Step4-1 Models

#Get some idea of the library sizes
sampleDepths <- sampleDepth(collapseFullCoverage(fullCov),0.75) 

#Define models
tst_adj_Comb <- expand.grid(tst = CriticalVar[[1]], adj = CriticalVar[[2]], stringsAsFactors = FALSE)
tst_adj_Comb <- tst_adj_Comb[tst_adj_Comb$tst != tst_adj_Comb$adj & !(tst_adj_Comb$tst == "ThreeGroup" & tst_adj_Comb$adj != 'NULL'),]
models <- map(seq_len(nrow(tst_adj_Comb)), ~ makeModels(sampleDepths = sampleDepths, testvars = pheno[[tst_adj_Comb[.x,1]]], adjustvars = pheno[[tst_adj_Comb[.x,2]]]))
    
names(models) <- apply(tst_adj_Comb, 1, str_flatten, collapse = "_")

##  Step4-2 Find candidate DERs

# Create an analysis directory
setwd(file.path("/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_V"))
for (p in as.list(file.path(list.dirs(path = "/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_V/Results",full.names = TRUE,recursive = FALSE )))) {
  if (!dir.exists(paste(p[[1]], 'analysisResults', sep = "/" ))) {
    dir.create(path = paste(p[[1]], 'analysisResults', sep = "/" ),recursive = TRUE)
  }
}
rm(p)
originalWd <- getwd()

##  Step4-3 Finding the best region Cuttoff value

# finding non-Overlapping exons
codingexons <- genomicState$codingGenome[mcols(genomicState$codingGenome)[,1] == "exon", ]
codingexonsByChr <- map(levels(seqnames(codingexons)), ~ ranges(codingexons[seqnames(codingexons) == .x,]))
names(codingexonsByChr) <-  seqnames(OarSeqinfo)
OverlappedPairs <- map(codingexonsByChr, ~ findOverlapPairs(.x,.x))
OverlappedExonsNames <- map(seq_len(27),~ as.data.frame(table(as.data.frame(OverlappedPairs[[..1]])[,4]))[as.data.frame(table(as.data.frame(OverlappedPairs[[..1]])[,4]))[,2] > 1,1])
NonOverlappedPairs <- map(seq_len(length(codingexonsByChr)), ~ codingexonsByChr[[..1]][! as.data.frame(codingexonsByChr[[..1]])[,4] %in% OverlappedExonsNames[[..1]]] )
names(NonOverlappedPairs) <-  seqnames(OarSeqinfo)

# Setting Delta for ER's alignment on each RegionMat(MCC)  "RegionMats_[[:alpha:]]*_[[:digit:]]*"

if (!file.exists(file.path("/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_V/SavedData/Deltas.RData"))) {
  Deltas <<-  map(str_extract_all(RegionMatFilesPath,"RegionMats_[[:alpha:]]*_.*"), ~ .x)
for (g in seq_len(length(RegionMatFilesPath))) {
  load(file = RegionMatFilesPath[g])
  regionMatVar <- get(ls()[which(str_detect(ls(), "RegionMats_[[:alpha:]]*_.*" ))])
   Deltas[[g]] <- map(seq_len(length(NonOverlappedPairs)), ~as_tibble(findOverlapPairs(ranges(NonOverlappedPairs[[..1]]), ranges(regionMatVar[[..1]][['regions']]))))
  rm(list = ls()[which(str_detect(ls(), "RegionMats_[[:alpha:]]*_.*" ))])
}
rm(g)
names(Deltas) <- str_replace(string =  list.files(path = "/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_V/SavedData/RegionMats", full.names = FALSE),
                             pattern = "[[:alpha:]]*[[:punct:]]+[[:alpha:]]*[[:punct:]]+(.*)[[:punct:]]+[[:alpha:]]*", replacement = "\\1")
Deltas <- modify(Deltas, ~ `names<-`(.x, names(NonOverlappedPairs)))
Deltas <- map(Deltas, ~ bind_rows(.x, .id = "Chr"))
Deltas <- map(Deltas, ~ mutate(.x, DeltaVal = abs(first.start - second.start) + abs(first.end - second.end)))
Deltas <- bind_rows(Deltas, .id = "MCC")
Deltas[["MCC"]] <- as.factor(Deltas[["MCC"]])
Deltas[["Chr"]] <- as.factor(Deltas[["Chr"]])
Deltas <- arrange(Deltas, c(MCC))
save(Deltas,file = file.path("/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_V/SavedData/Deltas.RData"))
} else{
  load(file = file.path("/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_V/SavedData/Deltas.RData"))
}


# Finding the Ideal MCC Value

Deltas <- arrange(Deltas, MCC)
DeltaValMedian_Zero <- Deltas %>% group_by(MCC) %>% summarise(DeltaValZero = sum(DeltaVal == 0) * 100/ n(), DeltasValMedian = median(DeltaVal)) %>% ungroup() %>%
                            mutate(IdealPoint  = DeltaValZero == max(DeltaValZero) & DeltasValMedian == min(DeltasValMedian)) %>% gather(key = "Variables", value = "Value", c(-MCC, -IdealPoint))

IdealPointPlot <- ggplot(DeltaValMedian_Zero, aes(x = MCC, y = Value)) + geom_point(aes(shape = IdealPoint, col = IdealPoint),size = 3.5) +
                                                                          ylim(c(53,54))
                                                                          #facet_wrap(~ DeltaValMedian_Zero$Variables, nrow = 2) + geom_line(aes(group =1))
  
IdealMCC <-  DeltaValMedian_Zero$MCC[which(DeltaValMedian_Zero$IdealPoint)[1]]


## Step4-4 Creating filteredCov using IdealMCC

IdealMCCfileName <- str_replace(list.files(path = "/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_V/SavedData/RegionMats", 
                                           full.names = FALSE)[which(str_detect(list.files(path = "/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_V/SavedData/RegionMats",full.names = FALSE),
                                                                                as.character(IdealMCC)))], pattern = "([[:alpha:]]*[[:punct:]]+[[:alpha:]]*[[:punct:]]+[[:digit:]]+[[:punct:]]{1}[[:digit:]]+)[[:punct:]]{1}[[:alpha:]]*", 
                                replacement = "\\1")

#file.create(path = paste(originalWd,"SavedData", "filteredCovs", paste0("filteredCov_",IdealMCCfileName,".RData"), sep = "/"), recursive = TRUE)
#filteredCov <-  lapply(fullCov, filterData, cutoff = as.numeric(IdealMCC), filter = "one" , TotalMapped = rep(1, length(bamfileslist)), targetSize = 1)
#save(filteredCov,file = file.path(paste(originalWd,"SavedData", "filteredCovs", paste0("filteredCov_",IdealMCCfileName,".RData"), sep = "/")))
load(file = paste0(SavedData, "RegionMats/", IdealMCCfileName, ".RData")) 
regionMat <-  get(ls()[which(str_detect(string = ls(),pattern = IdealMCCfileName) )])
rm(list = ls()[which(str_detect(string = ls(),pattern = IdealMCCfileName) )])


##  Step5-1 Save path of big dataframes

ResultsPath <- file.path("/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_V/Results")
if ( any(!names(models) %in% list.dirs(ResultsPath, recursive = TRUE, full.names = FALSE))) {
  map(paste0(ResultsPath, "/", names(models)[which(!names(models) %in% list.dirs(ResultsPath, recursive = TRUE, full.names = FALSE))]), ~dir.create(path = .x, recursive = TRUE))
}
setwd(ResultsPath)

## Implementing ER approach as suggested by  Leonardo Collado-Torres

##  Step3-2 Find DERs with DESeq2

#Round matrix
#counts <- lapply(names(regionMat), function(i) {round(regionMat[[i]]$coverageMatrix)})

# Round matrix and specify design
#dse <- lapply(counts, DESeqDataSetFromMatrix, colData = pheno, design = ~ Breedtrt + Fecuntrt)

# Perform DE analysis
#dse <- lapply(dse, DESeq, test = "LRT", fitType = "local", reduced = ~ Fecuntrt)
#names(dse) <- names(regionMat)

# Extract results
#deseq <- lapply(seq_along(regionMat), function(n) {cbind(regionMat[[n]]["regions"], results(dse[[n]]))})
#names(deseq) <- names(regionMat)


##  Step3-3 Find DERs with limma

mod <- model.matrix(~ pheno$Breedtrt + pheno$Fecuntrt)
mod0 <- model.matrix(~ pheno$Fecuntrt)

transformedCov <- map(regionMat, ~ log2(.x$coverageMatrix + 32))

fit <- map(transformedCov, ~lmFit(.x, mod))
fit0 <- map(transformedCov, ~lmFit(.x, mod0))

getF <- function(fit, fit0, theData){
  rss1 = rowSums((fitted(fit)-theData)^2)
  df1 = ncol(fit$coefficients)
  rss0 = rowSums((fitted(fit0)-theData)^2)
  df0 = ncol(fit0$coefficients)
  fstat = ((rss0-rss1)/(df1-df0))/(rss1/(ncol(theData)-df1))
  f_pval = pf(fstat, df1-df0, ncol(theData)-df1,lower.tail=FALSE)
  fout = cbind(fstat,df1-1,ncol(theData)-df1,f_pval)
  colnames(fout)[2:3] = c("df1","df0")
  fout = data.frame(fout)
  return(fout)
}

getF_arg_list <- list(fit, fit0, transformedCov)
names(getF_arg_list) <- c("fit","fit0","transformedCov")

ff <-  pmap(getF_arg_list, ~getF(fit = ..1, fit0 = ..2,theData = ..3 ))

## Get the p-value and assign it to the regions

limma <- map(regionMat,~.x$regions)
limma <- map2(.x = limma, .y = ff, ~ makeGRangesFromDataFrame(cbind.data.frame(.x, fstat =.y[["fstat"]]), keep.extra.columns = TRUE))
limma <- map2(.x = limma, .y = ff, ~ makeGRangesFromDataFrame(cbind.data.frame(.x, f_pval =.y[["f_pval"]]), keep.extra.columns = TRUE))
limma <- map2(.x = limma, .y = ff, ~ makeGRangesFromDataFrame(cbind.data.frame(.x, padj = p.adjust(.y[["f_pval"]], method = "fdr")), keep.extra.columns = TRUE))

DERs_limma <- bind_rows(map(limma, ~data.frame(indx = which(.x$padj < 0.1))),.id = "Chr")


# Acquiring fullRegion for each model design  
for (b in seq_along(names(models))) {
  if (!any(grepl("^fullRegions.Rdata$", dir(paste0(ResultsPath, "/", names(models)[b] ))))) {
    setwd(paste0(ResultsPath,"/", as.character(names(models)[b])))
    system.time( assign(paste0("Results_", as.character(names(models)[b]), "_", IdealMCCfileName) , map(seq_along(filteredCov), ~ analyzeChr(chr = names(filteredCov[.]),
                                                                                                                                             coverageInfo = filteredCov[[.]],
                                                                                                                                             models = models[[as.character(names(models)[b])]],
                                                                    groupInfo = pheno$Breedtrt, writeOutput = TRUE, cutoffFstat = 5e-02, cutoffType = "theoretical",
                                                                    nPermute = 50, seeds = 20190226 + seq_len(50), returnOutput = TRUE, txdb = TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1,
                                                                    runAnnotation = TRUE, genomicState = genomicState$fullGenome, annotationPackage = NULL, mc.cores = 7))))
    assign(paste0("Results_", as.character(names(models)[b]), "_", IdealMCCfileName),`names<-`( get(paste0("Results_", as.character(names(models)[b]), "_", IdealMCCfileName)), names(regionMat) ) )
   #names(results) <- names(regionMat)
    analyzeChrWarnings <- warnings()
    save( list = ls()[which(str_detect(ls(), paste0("Results_", as.character(names(models)[b]), "_", IdealMCCfileName)))] ,
         file = paste0(paste0(ResultsPath, "/", as.character(names(models)[b])), "/", ls()[which(str_detect(ls(), paste0("Results_", as.character(names(models)[b]), "_", IdealMCCfileName)))], ".RData"))
    mergeResults(chrs=c(1:26,"X"), prefix=".",
                 genomicState = genomicState$fullGenome, 
                 optionsStats = get(paste0("Results_", as.character(names(models)[b]), "_", IdealMCCfileName))[[1]]$optionsStats)
    rm(list = ls()[which(str_detect(ls(), paste0("Results_", as.character(names(models)[b]), "_", IdealMCCfileName)))])
    
  } else {
    print(paste0("fullRegions data file do exist in ",paste0(ResultsPath, "/", as.character(names(models)[b])), "."))
  }
}
rm(b)


##  5-2 Visually exploring the results


for (w in list.dirs("~/RNA-Seq/RProject/RNA-Seq_Oar3.1_V/Results", recursive = FALSE)) {
  setwd(w)
  load("./fullRegions.Rdata")
  load("./fullNullSummary.Rdata")
  load("./fullAnnotatedRegions.Rdata")
  load("./optionsMerge.Rdata")
  derfinderReport(prefix='.', browse=FALSE,
                  nBestRegions=25, makeBestClusters=TRUE, outdir='analysisResults',
                  fullCov= fullCov, optionsStats = optionsMerge$optionsStats, hg19 = FALSE,txdb = TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1, p.ideos = FALSE )
 
}
rm(w)



########################################################################################################
load("./fullRegions.Rdata")
load("./fullNullSummary.Rdata")
load("./fullAnnotatedRegions.Rdata")
load("./Results_Fecuntrt_NULL_mean_11.RData")
annoRegs <- annotateRegions(fullRegions, genomicState$fullGenome)
identical(annoRegs, fullAnnotatedRegions)

regionCov <- getRegionCoverage(fullCov, fullRegions)
load(file.path('Results', 'fullRegions.Rdata'))
seqinfo(fullRegions) <- OarSeqinfo
 
#############################  Take an overall and each Overview plots   fullRegions[as.logical(seqnames(fullRegions) == 1)] 
plotOverview(regions = keepSeqlevels(fullRegions, value = 1, pruning.mode = "tidy")  , annotation = Results_Breedtrt_Fecuntrt_mean_11$`1`$annotation,
             type = 'fwer' )
plotOverview(regions = fullRegions)
fullAnnotationDF <- map_dfr(seq_along(Results_Breedtrt_Fecuntrt_mean_11) , ~ Results_Breedtrt_Fecuntrt_mean_11[[.x]]$annotation, .id = "seqnames")
fullAnnotationDF$seqnames <- Rle(fullAnnotationDF$seqnames)

plotRegionCoverage(regions = fullRegions, regionCoverage = regionCov, 
                   groupInfo = pheno$ThreeGroup, nearestAnnotation = fullAnnotationDF, 
                   annotatedRegions = annoRegs, whichRegions=1:50, txdb = TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1, scalefac = 1, 
                   ask = FALSE)


###################################################################################################################
setwd("~/RNA-Seq/RProject/RNA-Seq_Oar3.1_V/Results/Breedtrt_NULL_2/")

#### I think the coverageinfo must include the filteredCov been filtered for IdealMCC value not cutoff value of 2
filteredCov_II <-  lapply(fullCov, filterData, cutoff = 1, filter = "one" , TotalMapped = rep(1, length(bamfileslist)), targetSize = 1)

Breedtrt_NULL_2nd <- map(seq_along(filteredCov_II), ~ analyzeChr(chr = names(filteredCov_II[.]),
                                            coverageInfo = filteredCov_II[[.]],
                                            models = models$Breedtrt_NULL,
                                            groupInfo = pheno$Breedtrt, writeOutput = TRUE, cutoffFstat = 0.99, cutoffType = "manual",
                                            nPermute = 10, seeds = 20190614 + seq_len(10), returnOutput = TRUE, txdb = TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1,
                                            runAnnotation = TRUE, genomicState = genomicState$fullGenome, annotationPackage = NULL, mc.cores = 7))

mergeResults(chrs=c(1:26,"X"), prefix=".",
             genomicState = genomicState$fullGenome,
             optionsStats = Breedtrt_NULL_2nd[[1]]$optionsStats )
load("./fullRegions.Rdata")
load("./fullNullSummary.Rdata")
load("./fullAnnotatedRegions.Rdata")
load("./optionsMerge.Rdata")

derfinderReport(prefix='.', browse=FALSE,
                nBestRegions=50, makeBestClusters=TRUE, outdir='analysisResults',
                fullCov= fullCov, optionsStats = optionsMerge$optionsStats, hg19 = FALSE,txdb = TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1, p.ideos = FALSE )














##  Step5-2 Merge results

setwd(originalWd)

# Merge results from several chromosomes.
mergeResults(chrs=c(1:26,"X"), prefix="Results",
             genomicState = genomicState$fullGenome, 
             optionsStats = results$`1`$optionsStats)

# Files created by mergeResults()
dir('Results', pattern = '.Rdata')

# fullRegions ---> The main result from mergeResults() is in fullRegions
load(file.path('Results', 'fullRegions.Rdata'))
seqinfo(fullRegions) <- OarSeqinfo

# fullAnotatedRegions
load(file.path('Results', 'fullAnnotatedRegions.Rdata'))


##### Step5 Visually explore results

# Find overlaps between regions and summarized genomic annotation
annoRegs <- suppressWarnings(annotateRegions(fullRegions, genomicState$fullGenome))

# Get the region coverage
regionCov <- getRegionCoverage(fullCov, fullRegions)
txdb <- TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1





















