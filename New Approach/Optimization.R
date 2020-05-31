##### Optimizing the CalculateDeltas Function #####

library("microbenchmark")
library("ggplot2movies")
library("profvis")
library("Rcpp")

## Function 1
IdealMCC_RegionMat_tibble <- function(IdealMCC_RegionMat) future_map(seq_along(IdealMCC_RegionMat), ~ as_tibble(IdealMCC_RegionMat[[.x]]$regions ))

# Function 2
nearest_regions_gap_filtered <- function(fullCov = fullCov, MCC = MCC_MRG_Comb[c,1], MRG = MCC_MRG_Comb[c,2]){
  IdealMCC_RegionMat   <-   regionMatrix(fullCov, cutoff = as.numeric(MCC), runfilter = FALSE ,L = c(rep(150,3), rep(100,2),rep(150, 4)), ### L argument can be given in mother function
                                         chrsStyle = "Ensembl", species = "ovis_aries", currentStyle = "Ensembl")
  Gaps_filtered_length <- future_map(IdealMCC_RegionMat,~ makeGRangesFromDataFrame(as.data.frame(gaps(.x$regions)[-c(1,2)])[as.data.frame(gaps(.x$regions)[-c(1,2)])$start > 1,],
                                                                                   keep.extra.columns = TRUE ) ) %>% future_map( ~sort.GenomicRanges(.x)) %>% 
    future_map(~ .x[which(width(.x) <= as.integer(MRG)),])# To exclude initial gaps
  
  
  
 list(nearest_regions_gap_filtered = future_map(seq_along(Gaps_filtered_length), ~ tibble(Leader_Region = precede(Gaps_filtered_length[[.x]], IdealMCC_RegionMat[[.x]][["regions"]])) %>% drop_na()),
       IdealMCC_RegionMat_tibble = IdealMCC_RegionMat_tibble(IdealMCC_RegionMat))
}

# Function 3
S_and_E_Regions <- function(nearest_regions_gap_filtered) {
  list(S_and_E_Regions = future_map2(nearest_regions_gap_filtered[["nearest_regions_gap_filtered"]], future_map(nearest_regions_gap_filtered[["nearest_regions_gap_filtered"]], ~ tibble(LeaderIndex = is.sequential(as.vector(.x[[1]]))) %>%
                                                                                                                  mutate( FollowerIndex = lag(LeaderIndex + 1, default = 1))),
                                     ~  tibble(Leader = ..1[[1]][..2[[1]]], Follower = as.integer(..1[[1]][..2[[2]]]-1) )),
       IdealMCC_RegionMat_tibble = nearest_regions_gap_filtered[["IdealMCC_RegionMat_tibble"]])
}
# Function 4
Raw_Gap_filled_IRanges <- function(S_and_E_Regions) future_map2(S_and_E_Regions[["IdealMCC_RegionMat_tibble"]], S_and_E_Regions[["S_and_E_Regions"]],
                                                                ~ bind_cols(..1[..2[[2]],c("seqnames","start","strand")], ..1[..2[[1]], c("end")]) %>% mutate(width = (end - start)+1) )

# Function 5
gap_filled_regions <- function(S_and_E_Regions) future_map(S_and_E_Regions[["S_and_E_Regions"]], ~as.integer(levels(as.factor(unlist(map2(.x[[2]], .x[[1]], ~ seq(from = ..1, to = ..2) ))))))

# Function 6
compiled_RegionMat <- function(Chr_names, Raw_Gap_filled_IRanges, S_and_E_Regions, gap_filled_regions){

  future_map(Chr_names,
             ~ list( regions = makeGRangesFromDataFrame( bind_rows(Raw_Gap_filled_IRanges[[which(Chr_names == as.character(.x))]],
                                                                   S_and_E_Regions[["IdealMCC_RegionMat_tibble"]][[which(Chr_names == as.character(.x))]][-c(gap_filled_regions[[which(Chr_names == as.character(.x))]]),
                                                                                                                                                                          c("seqnames", "start","strand","end","width")]))) )
}

# Function 7 # This Function takes too much memory
Final_RegionMat <- function(compiled_RegionMat, Chr_names) {
  names(compiled_RegionMat) <- Chr_names
  future_map(Chr_names, ~ list(regions = compiled_RegionMat[[.x]][["regions"]],
                                               coverageMatrix = getRegionCoverage(fullCov = fullCov[.x],regions = compiled_RegionMat[[.x]][["regions"]]),
                                               bpCoverage = map(seq_len(NROW(compiled_RegionMat[[.x]][["regions"]])),
                                                                .f =  function(r) {as.matrix(fullCov[[.x]][as.integer(as_tibble(compiled_RegionMat[[.x]][["regions"]])[r, "start"]): as.integer(as_tibble(compiled_RegionMat[[.x]][["regions"]])[r, "end"]),])} )), .progress = TRUE)
  
}

# Function 8
Add_Delta <- function(Final_RegionMat, Chr_names, OarSeqinfo){
names(Final_RegionMat) <- Chr_names
  for (b in seq_len(length(compiled_RegionMat))){
    seqinfo(compiled_RegionMat[[b]][["regions"]]) <- keepSeqlevels(OarSeqinfo, Chr_names[b])
  }
  rm(b)
  bind_rows(map(seq_len(length(NonOverlappedPairs)), ~as_tibble(findOverlapPairs(ranges(NonOverlappedPairs[[..1]]), ranges(Final_RegionMat[[..1]][['regions']])))), .id = "Chr") %>% 
    transmute( DeltaVal = abs(first.start - second.start) + abs(first.end - second.end)) %>% summarise(IntactDeltaPercent = mean(DeltaVal == 0))
   
}

# Function 9
CalculateDeltas <- function(fullCov, MCC_MRG_Comb, Chr_names, NonOverlappedPairs, OarSeqinfo ) {
  Deltas <- as.list(as.data.frame(matrix(1:nrow(MCC_MRG_Comb), ncol = nrow(MCC_MRG_Comb) , dimnames = list(  NULL  , paste0("RegionMats_", as.character(unite(MCC_MRG_Comb,
                                                                                                                                                              col = "names", sep = "_")[[1]]))   )  )))   
  for (c in seq_len(nrow(MCC_MRG_Comb))) {
    S_and_E_Regions <- S_and_E_Regions(nearest_regions_gap_filtered = nearest_regions_gap_filtered(fullCov = fullCov, MCC = MCC_MRG_Comb[c,1], MRG = MCC_MRG_Comb[c,2]))
    Raw_Gap_filled_IRanges <- Raw_Gap_filled_IRanges(S_and_E_Regions = S_and_E_Regions)
    Deltas[[c]] <<- Add_Delta(Final_RegionMat = Final_RegionMat(compiled_RegionMat =  compiled_RegionMat(Chr_names = Chr_names, Raw_Gap_filled_IRanges = Raw_Gap_filled_IRanges(S_and_E_Regions = S_and_E_Regions), S_and_E_Regions = S_and_E_Regions, 
                                             gap_filled_regions = gap_filled_regions(S_and_E_Regions = S_and_E_Regions))), Chr_names , OarSeqinfo)
    rm(list = objects()[which(objects() %in% c("IdealMCC_RegionMat","Gaps_Grange","Gaps_width","Gaps_filtered_length","IdealMCC_RegionMat_tibble","IdealMCC_RegionMat_Unified",
                                               "nearest_regions_gap_filtered", "CandidateRegionsIndices","S_and_E_Regions","Raw_Gap_filled_IRanges","gap_filled_regions","compiled_RegionMat"))])
    gc()
  }
  Deltas
}

CalculateDeltas(fullCov, MCC_MRG_Comb[1:3,],Chr_names = as.character(c(seq(1:26), "X")),NonOverlappedPairs = NonOverlappedPairs, OarSeqinfo = OarSeqinfo )


