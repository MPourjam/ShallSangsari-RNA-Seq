### Firstly
# finding gaps
###  Secondly
# finding nearest regions to gaps with median size on both sides

############# Finding Ideal MRG #####################

IdealMCC <- 10.4
load(paste0("/media/animalscience/nikookalam/MPourjam/RProject/RNA-Seq_Oar3.1_V/SavedData/RegionMats/", "RegionMats_mean_", as.character(IdealMCC),".RData"))
IdealMCC_RegionMat <- get(ls()[which(str_detect(ls(), paste0("RegionMats_mean_", as.character(IdealMCC)) ))])
rm(list = c(paste0("RegionMats_mean_", as.character(IdealMCC))))


Gaps_Grange <- makeGRangesFromDataFrame(map_dfr(IdealMCC_RegionMat,~ as.data.frame(gaps(.x$regions)[-c(1,2)])[as.data.frame(gaps(.x$regions)[-c(1,2)])$start > 1,]), # To exclude initial gaps
                                        keep.extra.columns = TRUE) ### The first two records for gaps seems to conatain length of chromosome
sort.GenomicRanges(Gaps_Grange)
Gaps_width <- tibble(width= width(Gaps_Grange))

#### Validity Tests ####
median(Gaps_width$width)
mean(Gaps_width$width)

ggplot(Gaps_width, aes(x = log2(width))) + geom_density() #+ xlim(c(0,5))
ggplot(Gaps_width, aes(x = ))
ggplot(subset.data.frame(Gaps_width, subset = width < 500), aes(x = width)) + geom_density()

#### Validity Tests ####
########################

filter_limits <- data.frame(maxl = seq(10,100, 10))
filter_limits <- unite(filter_limits, col = "range", remove = FALSE)

Gaps_filtered_length <- map(seq_len(nrow(filter_limits)), ~ Gaps_Grange[which(width(Gaps_Grange) <= filter_limits[.x,2]),])
names(Gaps_filtered_length) <- filter_limits$range


## reformating the Regionamt regions would be much more useful  ### CoverageToExon and getRegionCoverage would be a good option 

IdealMCC_RegionMat_tibble <- map_dfr(seq_along(IdealMCC_RegionMat), ~ as.tibble(IdealMCC_RegionMat[[.x]]$regions ))

IdealMCC_RegionMat_Unified <- makeGRangesFromDataFrame(IdealMCC_RegionMat_tibble, keep.extra.columns = TRUE)

nearest_regions_gap_filtered <- map(Gaps_filtered_length, ~ tibble(Leader_Region = precede(.x, IdealMCC_RegionMat_Unified)))

is.sequential <- function(vector, diff = 1){
  gap <- c()
  i <- seq_len(length(vector))
  for (iter in i[-length(i)]){
    if (vector[iter+ 1] == vector[iter] + diff){
      # it will go to upper iteration
    } else{
    gap <- append(gap, iter)
    }
  }
 append(gap, i[length(i)])
}

#map(seq_along(nearest_regions_gap_filtered), ~ tibble(LeaderIndex = is.sequential(as.vector(nearest_regions_gap_filtered[[.x]][[1]]))))
CandidateRegionsIndices <- map(map(seq_along(nearest_regions_gap_filtered),
                                   ~ tibble(LeaderIndex = is.sequential(as.vector(nearest_regions_gap_filtered[[.x]][[1]])))),
                               ~ mutate(.x, FollowerIndex = lag(c(.x[["LeaderIndex"]] + 1), default = 1) ))
# Retreiving the regions from the indices 
S_and_E_Regions <- map(seq_along(CandidateRegionsIndices)
    ,~ tibble(Leader = nearest_regions_gap_filtered[[.x]][[1]][CandidateRegionsIndices[[.x]][[1]]]
              ,Follower = as.integer(nearest_regions_gap_filtered[[.x]][[1]][CandidateRegionsIndices[[.x]][[2]]]-1) ))


### Creating candidate unifeid regions IRanges and Granges

Raw_Gap_filled_IRanges <- map(S_and_E_Regions, ~ bind_cols(IdealMCC_RegionMat_tibble[.x[[2]],c("seqnames","start","strand")],
                                                IdealMCC_RegionMat_tibble[.x[[1]], c("end")]) %>% mutate(width = (end - start)+1))

gap_filled_regions <- map(S_and_E_Regions, ~ as.integer(levels(as.factor(unlist(map2(.x[[2]], .x[[1]], ~ seq(from = ..1, to = ..2) ))))))  ### ??? why there are some regions repeated more than once ????? negative subsetting (i[-1]) doesn't mean t
map(gap_filled_regions, ~ IdealMCC_RegionMat_tibble[-c(.x),c("seqnames","start","end","width","strand")])

### Calculating number of not included in unified regions
nonUnified_Regions_No <- map_int(gap_filled_regions, ~ length(.x)  - NROW(IdealMCC_RegionMat_Unified))

compiled_RegionMat <- map(seq_along(Raw_Gap_filled_IRanges),
                          ~ makeGRangesFromDataFrame(bind_rows(Raw_Gap_filled_IRanges[[.x]],
                                                               IdealMCC_RegionMat_tibble[-c(gap_filled_regions[[.x]]), c("seqnames", "start","strand","end","width")])  ))


for (c in seq_along(compiled_RegionMat)){
  seqinfo(compiled_RegionMat[[c]]) <- OarSeqinfo
}
rm(c)


############# Finding Ideal MRG #####################
#####################################################




############# Memory Management Tricks ##############

.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.size <- napply(names, function(x) {object.size(x)/1024^3})
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.dim)
  names(out) <- c("Type", "Size", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}
lsos <- function(..., n=10) {
  .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
}

############# Memory Management Tricks ##############
#####################################################



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
library("pryr")
library("furrr")


plan("multisession")
options(future.globals.maxSize= 8912896000)
options(future.globals.onMissing = "ignore")


# Get pheno table
Fecuntrt <- c(c(rep("LF",3), "HF",rep("LF",2), rep("HF",2),"LF"))
Breedtrt <- c(c(rep("San",3), rep("Shall",6)))
ThreeGroup <- c(rep("San", 3), rep("Shall_OUT",2), rep("Shall",3), "Shall_OUT")
ReadLength <- c(rep(150,3), rep(100,2),rep(150, 4))
pheno <- as.data.frame( cbind(Breed = Breedtrt, Fecun = Fecuntrt, ThreeGroup = ThreeGroup, ReadLength = ReadLength))
pheno <- map_dfc(colnames(pheno), ~ as.factor(pheno[,.]))
colnames(pheno) <- c("Breedtrt","Fecuntrt","ThreeGroup","ReadLength")
rownames(pheno) <- NULL

##  Step3-1 Via regionMatrix()

Breeds <- c("Shall","San")
CriticalVarsList <-  list(MRG =  c(seq(10,100, 10)), MCC = c(seq(1,20,0.5)))
MCC_MRG_Comb <- expand_grid(MCC = CriticalVarsList$MCC, MRG = CriticalVarsList$MRG)
SavedData <- file.path("/home/animalscience/MPourjam/")

### Creating RegionMats Directory
if (!"RegionMats" %in% list.dirs(SavedData, full.names = FALSE)) {
  dir.create(path = paste0(SavedData, "RegionMats"))
  RegionMatsPath <<- file.path(paste0(SavedData, "RegionMats/"))
} else {
  RegionMatsPath <<- file.path(paste0(SavedData, "RegionMats/"))
}
###

## Step2-1 Loading and creating Raw Data

options(species = "Ovis_aries")

options(chrsStyle = "Ensembl")

extendedMapSeqlevels(seqnames = c(1:26, "X"), style = "Ensembl", currentStyle = "Ensembl")                               

OarSeqinfo <- Seqinfo(seqnames = c(1:26, "X"), seqlengths = c( 275612895, 248993846, 224283230, 119255633, 107901688, 117031472, 
                                                               100079507, 90695168, 94726778, 86447213, 62248096, 79100223, 83079144, 62722625, 80923592,
                                                               71719816, 72286588, 68604602, 60464314, 51176841, 50073674, 50832532, 62330649, 42034648,
                                                               45367442, 44077779, 135437088),isCircular = rep(FALSE, 27), genome = "Oar3.1")

bamfilespath <- list.files(path = "/media/animalscience/nikookalam/MPourjam/DATA/Analyze/Out/Sheep"
                           , pattern = ".Aligned.out.bam", all.files = TRUE, full.names = TRUE
                           , recursive = TRUE )

readsnames <- gsub(".*/(.*)_Aligned.out.bam", "\\1", bamfilespath)
names(bamfilespath) <- readsnames

Coordbamfilespath <- list.files(path = "/media/animalscience/nikookalam/MPourjam/DATA/Analyze/Out/Sheep"
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

#Mart <- useMart(biomart = "ensembl_mart_97", path = "/biomart/martservice", port = 80 ,archive = FALSE)
#ensembl = useMart("ensembl", dataset = "oaries_gene_ensembl")

#txdb <- TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1 <- makeTxDbPackageFromBiomart(taxonomyId = 9940, port = 80, biomart = "ENSEMBL_MART_ENSEMBL"
#                                                                                     , host = "www.ensembl.org", dataset = "oaries_gene_ensembl"
#                                                                                     , filter = NULL, circ_seqs = "MT", id_prefix="ensembl_", miRBaseBuild = NA
#                                                                                     , version = "1.0"
#                                                                                     , maintainer = "Mohsen Pourjam <Pourjam.cs@hotmail.com>"
#                                                                                     , author = "Mohsen Pourjam")
#install.packages("./SavedData/TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1", repos = NULL, type = "source")
library( "TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1")
keepSeqlevels(TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1::TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1,c(1:26, "X"))
#plotkeepSeqlevels(TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1,c(1:26, "X"))
#save(TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1, file = paste0(SavedData, ls()[which(str_detect(ls(), "^.*TxDb.*aries.*$"))], ".RData", sep = ""))
#library(TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1)



SavedData <- file.path("/media/animalscience/nikookalam/MPourjam/RProject/RNA-Seq_Oar3.1_V/SavedData") ### Changing SavedData path to load fullcov 
if (any(grepl("^fullCov\\.RData$", dir(SavedData)))) {
  load(paste(SavedData, list.files(SavedData, "^fullCov\\.RData$"), sep = "/"))
} else {
  fullCov <- fullCoverage(bamfileslist, chrs = c(rep(1:26), "X"), totalMapped = TotalMapped, targetSize = 4e+07, returnMean = TRUE)
  save(fullCov, file = paste0(SavedData, ls()[which(str_detect(ls(), "^fullCov$"))], ".RData", sep = ""))
}
SavedData <- file.path("/home/animalscience/MPourjam/") ### Changing back to SavedData defined in line 196




genomicState <- makeGenomicState(txdb =  TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1::TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1
                                 , chrs = names(seqlengths(bamfileslist)[1:27])
                                 , style = getOption("chrsStyle", "Ensembl"), currentStyle = "Ensembl")

seqinfo(genomicState) <- OarSeqinfo

codingexons <- genomicState$codingGenome[mcols(genomicState$codingGenome)[,1] == "exon", ]
codingexonsByChr <- map(levels(seqnames(codingexons)), ~ ranges(codingexons[seqnames(codingexons) == .x,]))
names(codingexonsByChr) <-  seqnames(OarSeqinfo)
OverlappedPairs <- map(codingexonsByChr, ~ findOverlapPairs(.x,.x))
OverlappedExonsNames <- map(seq_len(27),~ as.data.frame(table(as.data.frame(OverlappedPairs[[..1]])[,4]))[as.data.frame(table(as.data.frame(OverlappedPairs[[..1]])[,4]))[,2] > 1,1])
NonOverlappedPairs <- map(seq_len(length(codingexonsByChr)), ~ codingexonsByChr[[..1]][! as.data.frame(codingexonsByChr[[..1]])[,4] %in% OverlappedExonsNames[[..1]]] )
names(NonOverlappedPairs) <-  seqnames(OarSeqinfo)


is.sequential <- function(vector, diff = 1){
  gap <- c()
  i <- seq_len(length(vector))
  for (iter in i[-length(i)]){
    if (vector[iter+ 1] == vector[iter] + diff){
      # it will go to upper iteration
    } else{
      gap <- append(gap, iter)
    }
  }
  append(gap, i[length(i)])
}


 #CreateRegionmats <- function(MCC_MRG_Comb = MCC_MRG_Comb, fullCov , SavePath){
  for (c in seq_len(nrow(MCC_MRG_Comb))){
    gc()
    IdealMCC_RegionMat   <<-   regionMatrix(fullCov, cutoff = as.numeric(MCC_MRG_Comb[c,1]), runfilter = FALSE ,L = c(rep(150,3), rep(100,2),rep(150, 4)), ### L argument can be given in mother function
                                           chrsStyle = "Ensembl", species = "ovis_aries", currentStyle = "Ensembl")
    save(IdealMCC_RegionMat, file = file.path(paste0(SavedData,"RegionMat_MCC-", as.character(MCC_MRG_Comb[c,1]),"_MRG-",as.character(MCC_MRG_Comb[c,2]),".RData" )))
    rm(IdealMCC_RegionMat, pos = 1)
    }
}
### Test ###
#CreateRegionmats(MCC_MRG_Comb = MCC_MRG_Comb[1:10,], fullCov = fullCov, SavePath = SavedData)
#####





CalculateDeltas <- function(MCC_MRG_Comb = MCC_MRG_Comb, fullCov , NonOverlappedPairs, SavePath) {
  Deltas <<- as.list(as.data.frame(matrix(1:nrow(MCC_MRG_Comb), ncol = nrow(MCC_MRG_Comb) , dimnames = list(  NULL  , paste0("RegionMats_", as.character(unite(MCC_MRG_Comb,
                                                                                            col = "names", sep = "_")[[1]]))   )  )))
  # The "future" package does not work properly when a global variables are refered so that I am using map instead of future_map
  map(seq_len(nrow(MCC_MRG_Comb)), .f = function(c) {
    #  if (!any(str_detect(pattern = paste0(paste("RegionMats", MCC_MRG_Comb[c,1], MCC_MRG_Comb[c,2], sep = "_"),".RData"),string = list.files(RegionMatsPath, full.names = FALSE)))) {
    ### Creating RegjionMat using MCC 
    #  assign(paste("RegionMats", as.character(MCC_MRG_Comb[c,1]), as.character(MCC_MRG_Comb[c,2]), sep = "_"), 
    gc()
    IdealMCC_RegionMat   <<-   regionMatrix(fullCov, cutoff = as.numeric(MCC_MRG_Comb[c,1]), runfilter = FALSE ,L = c(rep(150,3), rep(100,2),rep(150, 4)), ### L argument can be given in mother function
                                           chrsStyle = "Ensembl", species = "ovis_aries", currentStyle = "Ensembl") #)
    
    ### Modifying RegionMat using MRG 
    
    Gaps_Grange <<- map(IdealMCC_RegionMat,~ makeGRangesFromDataFrame(as.data.frame(gaps(.x$regions)[-c(1,2)])[as.data.frame(gaps(.x$regions)[-c(1,2)])$start > 1,],keep.extra.columns = TRUE ) )# To exclude initial gaps
    ### The first two records for gaps seems to conatain length of chromosome
    map(Gaps_Grange, ~sort.GenomicRanges(.x))
    Gaps_width <<- map(Gaps_Grange, ~tibble(width= width(.x)))
    
    Gaps_filtered_length <<-  map(Gaps_Grange, ~ .x[which(width(.x) <= as.integer(MCC_MRG_Comb[c,2])),])
    
    # reformating the Regionamt regions would be much more useful  ### "getRegionCoverage" would be a good option 
    
    IdealMCC_RegionMat_tibble <<- map(seq_along(IdealMCC_RegionMat), ~ as.tibble(IdealMCC_RegionMat[[.x]]$regions ))
    
    IdealMCC_RegionMat_Unified <<- makeGRangesFromDataFrame(map_dfr(IdealMCC_RegionMat_tibble,~.x), keep.extra.columns = TRUE)
    
    nearest_regions_gap_filtered <<- map(seq_along(Gaps_filtered_length), ~tibble(Leader_Region = precede(Gaps_filtered_length[[.x]], IdealMCC_RegionMat[[.x]][["regions"]])) %>% drop_na())
    
    #map(seq_along(nearest_regions_gap_filtered), ~ tibble(LeaderIndex = is.sequential(as.vector(nearest_regions_gap_filtered[[.x]][[1]]))))
    CandidateRegionsIndices <<- map(nearest_regions_gap_filtered, ~ tibble(LeaderIndex = is.sequential(as.vector(.x[[1]]))) %>% mutate( FollowerIndex = lag(LeaderIndex + 1, default = 1)))
    # Retreiving the regions from the indices 
    S_and_E_Regions <<- map2(nearest_regions_gap_filtered, CandidateRegionsIndices, ~  tibble(Leader = ..1[[1]][..2[[1]]]
                                                                                             ,Follower = as.integer(..1[[1]][..2[[2]]]-1) ))
    
    
    ### Creating candidate unifeid regions IRanges and Granges
    
    Raw_Gap_filled_IRanges <<-  map2(IdealMCC_RegionMat_tibble, S_and_E_Regions, ~ bind_cols(..1[..2[[2]],c("seqnames","start","strand")], ..1[..2[[1]], c("end")]) %>% mutate(width = (end - start)+1))
    
    gap_filled_regions <<- map(S_and_E_Regions, ~as.integer(levels(as.factor(unlist(map2(.x[[2]], .x[[1]], ~ seq(from = ..1, to = ..2) )))))  )
    
    
    #### We need to reconstruct the same hierarchy as RegionMatrix Outpu (i.e : IdealMCC_Regionmat)
    compiled_RegionMat <<-   map(names(IdealMCC_RegionMat),
                                       ~ list( regions = makeGRangesFromDataFrame( bind_rows(Raw_Gap_filled_IRanges[[which(names(IdealMCC_RegionMat) == as.character(.x))]],
                                                                                             IdealMCC_RegionMat_tibble[[which(names(IdealMCC_RegionMat) == as.character(.x))]][-c(gap_filled_regions[[which(names(IdealMCC_RegionMat) == as.character(.x))]]),
                                                                                                                                                                               c("seqnames", "start","strand","end","width")]))) )
    names(compiled_RegionMat) <<- names(IdealMCC_RegionMat)
    compiled_RegionMat <<- map(names(compiled_RegionMat), ~ list(regions = compiled_RegionMat[[.x]][["regions"]],
                                                                       coverageMatrix = getRegionCoverage(fullCov = fullCov[.x],regions = compiled_RegionMat[[.x]][["regions"]]),
                                                                       bpCoverage = map(seq_len(NROW(compiled_RegionMat[[.x]][["regions"]])),
                                                                                               .f =  function(r) {as.matrix(fullCov[[.x]][as.integer(as_tibble(compiled_RegionMat[[.x]][["regions"]])[r, "start"]): as.integer(as_tibble(compiled_RegionMat[[.x]][["regions"]])[r, "end"]),])} )), .progress = TRUE)
    
    names(compiled_RegionMat) <<- names(IdealMCC_RegionMat)
    
    for (b in seq_len(length(compiled_RegionMat))){
      seqinfo(compiled_RegionMat[[b]][["regions"]]) <<- keepSeqlevels(OarSeqinfo, names(compiled_RegionMat)[b])
    }
    rm(b)
    
    ################ Calculating Delta for the created RegionMat (compiled_RegionMat) ################
    
    Deltas[[paste("RegionMats", as.character(MCC_MRG_Comb[c,1]), as.character(MCC_MRG_Comb[c,2]), sep = "_")]] <<-  bind_rows(map(seq_len(length(NonOverlappedPairs)),
                                                                                                                                        ~as_tibble(findOverlapPairs(ranges(NonOverlappedPairs[[..1]]),
                                                                                                                                                                    ranges(compiled_RegionMat[[..1]][['regions']])))), .id = "Chr") %>% transmute( DeltaVal = abs(first.start - second.start) + abs(first.end - second.end)) %>% summarise(IntactDeltaPercent = mean(DeltaVal == 0))
    
    
    
    
    ################ Calculating Delta for the created RegionMat (compiled_RegionMat) ################
    
    rm(list = objects()[which(objects() %in% c("IdealMCC_RegionMat","Gaps_Grange","Gaps_width","Gaps_filtered_length","IdealMCC_RegionMat_tibble","IdealMCC_RegionMat_Unified",
                "nearest_regions_gap_filtered", "CandidateRegionsIndices","S_and_E_Regions","Raw_Gap_filled_IRanges","gap_filled_regions","compiled_RegionMat"))], pos = 1)
    
    #  }
    gc()
    print(paste0("RegionMats data of cutoff value ", MCC_MRG_Comb[c,1], " and MRG value of ", MCC_MRG_Comb[c,2] , " has been added to Deltas."))
  })
  
  gc()
  #Deltas <- map(Deltas, ~ bind_rows(.x, .id = "Chr"))
  #Deltas <- map(Deltas, ~ transmute(.x, DeltaVal = abs(first.start - second.start) + abs(first.end - second.end))) %>% summarise(sum(DeltaVal)) 
  Deltas <<- bind_rows(Deltas, .id = "MCC_MRG")
  Deltas[["MCC_MRG"]] <<- as.factor(Deltas[["MCC_MRG"]])
  #Deltas[["Chr"]] <- as.factor(Deltas[["Chr"]])
  Deltas <<- arrange(Deltas, c(MCC_MRG))


### Saving Stage 
file.create(file.path(paste0(SavePath,"Deltas",sep = "/"), ".RData"))
save(list = c("Deltas") , file = file.path(paste0(paste(SavePath,sep = "/"), ".RData")))

}


debug(CalculateDeltas)
CalculateDeltas(MCC_MRG_Comb = MCC_MRG_Comb[2,],fullCov = fullCov, NonOverlappedPairs = NonOverlappedPairs,SavePath = "/home/animalscience/MPourjam/RegionMats" )








################################################################################################################################
RegionMatFilesPath <- list.files(path = "/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_V/SavedData/RegionMats", full.names = TRUE)

# Because there isn't enough space on the external  HDD I temporarily use the internal one
RegionMats_Dir <- file.path(path = "/home/animalscience/MPourjam/RegionMats/")
RegionMatFilesPath <- list.files(path = "/home/animalscience/MPourjam/RegionMats", full.names = TRUE)


if (!file.exists(file.path("/media/animalscience/nikookalam/MPourjam/RProject/RNA-Seq_Oar3.1_V/SavedData/Deltas.RData"))) {
  Deltas <<-  map(str_extract_all(RegionMatFilesPath,"RegionMats_[[:alpha:]]*_.*"), ~ .x)
  for (g in seq_len(length(RegionMatFilesPath))) {
    load(file = RegionMatFilesPath[g])
    regionMatVar <- get(ls()[which(str_detect(ls(), "RegionMats_[[:alpha:]]*_.*" ))])
    Deltas[[g]] <- map(seq_len(length(NonOverlappedPairs)), ~as_tibble(findOverlapPairs(ranges(NonOverlappedPairs[[..1]]), ranges(regionMatVar[[..1]][['regions']]))))
    rm(list = ls()[which(str_detect(ls(), "RegionMats_[[:alpha:]]*_.*" ))])
  }
  rm(g)
  names(Deltas) <- str_replace(string =  list.files(path = "/media/animalscience/nikookalam/MPourjam/RProject/RNA-Seq_Oar3.1_V/SavedData/RegionMats", full.names = FALSE),
                               pattern = "[[:alpha:]]*[[:punct:]]+[[:alpha:]]*[[:punct:]]+(.*)[[:punct:]]+[[:alpha:]]*", replacement = "\\1")
  Deltas <- modify(Deltas, ~ `names<-`(.x, names(NonOverlappedPairs)))
  Deltas <- map(Deltas, ~ bind_rows(.x, .id = "Chr"))
  Deltas <- map(Deltas, ~ mutate(.x, DeltaVal = abs(first.start - second.start) + abs(first.end - second.end)))
  Deltas <- bind_rows(Deltas, .id = "MCC")
  Deltas[["MCC"]] <- as.factor(Deltas[["MCC"]])
  Deltas[["Chr"]] <- as.factor(Deltas[["Chr"]])
  Deltas <- arrange(Deltas, c(MCC))
  save(Deltas,file = file.path("/media/animalscience/nikookalam/MPourjam/RProject/RNA-Seq_Oar3.1_V/SavedData/Deltas.RData"))
} else{
  load(file = file.path("/media/animalscience/nikookalam/MPourjam/RProject/RNA-Seq_Oar3.1_V/SavedData/Deltas.RData"))
}



