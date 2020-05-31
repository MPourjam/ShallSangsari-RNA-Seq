### Firstly
# finding gaps
###  Secondly
# finding nearest regions to gaps with median size on both sides

############# Finding Ideal MRG #####################

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

### Creating RegionMats Direct
if (!"RegionMats" %in% list.dirs(SavedData, full.names = FALSE)) {
  dir.create(path = paste0(SavedData, "RegionMats"))
  RegionMatsPath <<- file.path(paste0(SavedData, "RegionMats/"))
} else {
  RegionMatsPath <<- file.path(paste0(SavedData, "RegionMats/"))
}
###




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

for (c in seq_len(nrow(MCC_MRG_Comb))){
  if (!any(str_detect(pattern = paste0(paste("RegionMats", MCC_MRG_Comb[c,1], MCC_MRG_Comb[c,2], sep = "_"),".RData"),string = list.files(RegionMatsPath, full.names = FALSE)))) {
    ### Creating RegjionMat using MCC 
    assign(paste("RegionMats", as.character(MCC_MRG_Comb[c,1]), as.character(MCC_MRG_Comb[c,2]), sep = "_"), 
           regionMatrix(fullCov, cutoff = as.numeric(MCC_MRG_Comb[c,1]), runfilter = FALSE ,L = c(rep(150,3), rep(100,2),rep(150, 4)),
                        chrsStyle = "Ensembl", species = "ovis_aries", currentStyle = "Ensembl"))
    
    ### Modifying RegionMat using MRG 
    IdealMCC_RegionMat <- get(ls()[which(str_detect(ls(), paste("RegionMats", as.character(MCC_MRG_Comb[c,1]), as.character(MCC_MRG_Comb[c,2]), sep = "_") ))])
    rm(list = c(paste("RegionMats", as.character(MCC_MRG_Comb[c,1]), as.character(MCC_MRG_Comb[c,2]), sep = "_") ))
    Gaps_Grange <- makeGRangesFromDataFrame(map_dfr(IdealMCC_RegionMat,~ as.data.frame(gaps(.x$regions)[-c(1,2)])[as.data.frame(gaps(.x$regions)[-c(1,2)])$start > 1,]), # To exclude initial gaps
                                            keep.extra.columns = TRUE) ### The first two records for gaps seems to conatain length of chromosome
    sort.GenomicRanges(Gaps_Grange)
    Gaps_width <- tibble(width= width(Gaps_Grange))
    
    Gaps_filtered_length <-  Gaps_Grange[which(width(Gaps_Grange) <= as.integer(MCC_MRG_Comb[1,2])),]
    
    # reformating the Regionamt regions would be much more useful  ### "getRegionCoverage" would be a good option 
    
    IdealMCC_RegionMat_tibble <- map_dfr(seq_along(IdealMCC_RegionMat), ~ as.tibble(IdealMCC_RegionMat[[.x]]$regions ))
    
    IdealMCC_RegionMat_Unified <- makeGRangesFromDataFrame(IdealMCC_RegionMat_tibble, keep.extra.columns = TRUE)
    
    nearest_regions_gap_filtered <- tibble(Leader_Region = precede(Gaps_filtered_length, IdealMCC_RegionMat_Unified))
    
    #map(seq_along(nearest_regions_gap_filtered), ~ tibble(LeaderIndex = is.sequential(as.vector(nearest_regions_gap_filtered[[.x]][[1]]))))
    CandidateRegionsIndices <- tibble(LeaderIndex = is.sequential(as.vector(nearest_regions_gap_filtered[[1]]))) %>% mutate( FollowerIndex = lag(LeaderIndex + 1, default = 1))
    # Retreiving the regions from the indices 
    S_and_E_Regions <- tibble(Leader = nearest_regions_gap_filtered[[1]][CandidateRegionsIndices[[1]]]
                              ,Follower = as.integer(nearest_regions_gap_filtered[[1]][CandidateRegionsIndices[[2]]]-1) )
    
    
    ### Creating candidate unifeid regions IRanges and Granges
    
    Raw_Gap_filled_IRanges <- bind_cols(IdealMCC_RegionMat_tibble[S_and_E_Regions[[2]],c("seqnames","start","strand")],
                                        IdealMCC_RegionMat_tibble[S_and_E_Regions[[1]], c("end")]) %>% mutate(width = (end - start)+1)
    
    gap_filled_regions <- as.integer(levels(as.factor(unlist(map2(S_and_E_Regions[[2]], S_and_E_Regions[[1]], ~ seq(from = ..1, to = ..2) )))))  
    
    
    
    compiled_RegionMat <- makeGRangesFromDataFrame(bind_rows(Raw_Gap_filled_IRanges,
                                                             IdealMCC_RegionMat_tibble[-c(gap_filled_regions), c("seqnames", "start","strand","end","width")])  )
    
    
    seqinfo(compiled_RegionMat) <- OarSeqinfo
    
    
    ############### Calculating Delata Value can be done here instead of being done seperately #######################
    
    
    ##################################################################################################################
    
    ### Saving Stage 
    file.create(file.path(paste0(paste(RegionMatsPath, paste("RegionMats", MCC_MRG_Comb[c,1], MCC_MRG_Comb[c,2], sep = "_"),sep = "/"), ".RData")))
    save(list = ls()[which(str_detect(ls(), paste("RegionMats", CF_Value_Type[c,1], CF_Value_Type[c,2], sep = "_")))],
         file = file.path(paste0(paste(RegionMatsPath, paste("RegionMats", CF_Value_Type[c,1], CF_Value_Type[c,2], sep = "_"),sep = "/"), ".RData")))
    rm(list = ls()[which(str_detect(ls(), paste("RegionMats", CF_Value_Type[c,1], CF_Value_Type[c,2], sep = "_")))]) 
    print(paste0("RegionMats data of cutoff value ", CF_Value_Type[c,2], " and (", CF_Value_Type[c,1] , ") filter policy has been saved."))
  }
}
rm(c)

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



