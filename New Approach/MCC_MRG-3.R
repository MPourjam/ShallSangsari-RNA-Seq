##### Prerequisites

### Libraries

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


### Path
SavePath <- "/HDD/DB/" ### This the root of our database ### The terminal '/' is necessary

### Functions
is.sequential <- function(vector, diff = 1){
  i <- seq_len(length(vector))
  gap <- numeric()
  for (iter in i[-length(i)]){
    if (vector[iter+ 1] == vector[iter] + diff){
      # it will go to upper iteration
    } else{
    c(gap, iter)
    }
  }
 c(gap, i[length(i)])
}


### Variables

MCC <- c(seq(1,20,0.5))
MRG <- c(seq(10,100, 10)


#### We need to show that regions don't usually lay on the edges of chromosomes ####


#### We need to show that regions don't usually lay on the edges of chromosomes ####

#### What if gaps lay on either one end or both of ends????? #####

gaps_Grange <- function(regionMat_region_Grange){
  gaps <- future_map(regionMat_region_Grange, ~ as.data.frame(gaps(.x))[-c(1,2)])
  gaps <- future_map(gaps,~makeGRangesFromDataFrame(.x[.x$start > 1,], keep.extra.columns = TRUE ))
  gaps <- future_map(gaps, ~ sort.GenomicRanges(.x))
  gaps
}

exclude_compiled_regions <- function(PrecedingRegion_each_Chr){  
rows <- nrow(PrecedingRegion_each_Chr)
exluding_regions <- unique(unlist(map(seq_len(rows), ~ seq(as.numeric(PrecedingRegion_each_Chr[.x,"FollowerRegion"]) ,
 as.numeric(PrecedingRegion_each_Chr[.x,"LeaderRegion"]) ))))
exluding_regions
}

NewRegions <- function(RegionMat_Path, MRG){  
  stopifnot(!file.exit(RegionMat_Path))         
  MRG_RegionSet <- vector(mode = 'list',length = 2)
  SaveDir <- paste0(dirname(RegionaMat_Path),"/")
  RegionMat <- load(file = RegionMat_Path)
  .subset2(MRG_RegionSet,1) <- future_map(RegionMat, ~ .subset2(.x, 1)) #### The Structure of RegionMat is RegionMat$Chr$regions
  .subset2(MRG_RegionSet,2) <- future_map(.subset2(MRG_RegionSet, 1), ~ gaps_Grange(.x))
  stopifnot(file.exit(paste0(SaveDir,"MRG.rds")))
  load(file.path(paste0(SaveDir,"MRG.rds")))
  stopifnot(MRG <= 100)

    gc()

    ## Filtering gaps
    gap_filtered <- future_map(.subset2(MRG_RegionSet,2), ~ .x[which(width(.x) <= as.integer(MRG)),])

    ## Finding receding regions
    preceding_region <- future_map2(gap_filtered, .subset2(MRG_RegionSet,1), ~ tibble(Leader_Region = precede(..1, ..2 )) %>% drop_na() )
    # preceding_region is a list of 27(number of chromosomes) elements

    ### 
    Leader_Follower_Index <- future_map(preceding_region, ~ tibble(LeaderIndex = is.sequential(.subset2(.x,1))) %>%
     mutate(FollowerIndex = dplyr::lag(LeaderIndex + 1, default = 1)) )  #??? What about the last one in preceding_region

    ### Retrieving corresponding regoins using indices from preceding regions vector (list of 27)
    Leader_Follower_PrecedingRegion <- future_map2(preceding_region, Leader_Follower_Index,
     ~tibble(LeaderRegion = ..1[.subset2(..2, 1)]+1 ), FollowerRegion = ..1[.subset2(..2,2)])
    names(Leader_Follower_PrecedingRegion) <- paste0("Chr", as.character(c(1:26,"X")))

    ### Constituting new regioins
    New_Regions_Seqnames <- future_map2(Leader_Follower_PrecedingRegion, .subset2(MRG_RegionSet, 1),
     ~ seqnames(..2[.subset2(..1, 2),]))
    New_Regions_Strand <- future_map2(Leader_Follower_PrecedingRegion, .subset2(MRG_RegionSet, 1),
      ~strand(..2[.subset2(..1, 2),])
    New_Regions_Ranges <- future_map2(Leader_Follower_PrecedingRegion, .subset2(MRG_RegionSet, 1),
      ~IRanges(start = start(..2[[.subset2(..1, 2)]]),
     end = end(..2[.subset2(..1, 1),])))
    New_Regions_GRange <- pmap(New_Regions_Seqnames, New_Regions_Strand, New_Regions_Ranges,
      ~ GRanges(seqnames = ..1, ranges = ..3, strand = ..2)

    ### Removing compiled regions from the list of regions
  exluding_regions_allChrs <- future_map(Leader_Follower_PrecedingRegion, ~ exclude_compiled_regions)

  compiled_regions_excluded_tibble <- future_map2(.subset2(MRG_RegionSet,1) , exluding_regions_allChr,
   ~ as_tibble(..1[-c(..2),c("seqnames","ranges","strand")])) ### Might incur problems due to nature of regions as a Grange (not subsettable)
  new_region_set <- future_map2(compiled_regions_excluded_tibble, New_Regions_GRange,
   ~ makeGRangesFromDataFrame(bind_rows(..1,..2)))
saveRDS(object= new_region_set,
  file = file.path(paste0(SaveDir, paste(basename(SaveDir), as.character(MRG), sep = "_" ), ".rds" )))

### Updating MRG
MRG <- MRG + 10
saveRDS( object = MRG, file = file.path(paste0(SaveDir,"MRG.rds")))

### Updating RegionMat_Path.rds
whichRegionMat <- which(str_detect(list.files(RegionMatsPath, recursive=TRUE, full.names=TRUE),
 pattern = RegionMat_Path))
  RegionMat_Path <- list.files(RegionMatsPath, recursive=TRUE, full.names=TRUE)[whichRegionMat+1]
  saveRDS(object = as.character(RegionMat_Path), file = paste0(SaveDir, "RegionMat_Path.rds"))
  }
########## Before we run we need to create initial RegionMat_Path.rds and MRG.rds
}





