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


CalculateDeltas <- function(MCC_MRG_Comb = MCC_MRG_Comb, fullCov , NonOverlappedPairs, SavePath) {
  Deltas <<- as.list(as.data.frame(matrix(1:nrow(MCC_MRG_Comb), ncol = nrow(MCC_MRG_Comb) , dimnames = list(  NULL  , paste0("RegionMats_", as.character(unite(MCC_MRG_Comb,
                                                                                                                                                               col = "names", sep = "_")[[1]]))   )  )))
  # The "future" package does not work properly when a global variables are refered so that I am using map instead of future_map
  

    map(seq_len(nrow(MCC_MRG_Comb)), .f = function(c) {
    
     Final_RegionMat <- function() {
       compiled_RegionMat <- function(){
        
        S_and_E_Regions <- function(){
          
                                                                ######### Maybe I need to define specific data objects by various slots ##################################
          nearest_regions_gap_filtered <- function(){
            IdealMCC_RegionMat   <-   regionMatrix(fullCov, cutoff = as.numeric(MCC_MRG_Comb[c,1]), runfilter = FALSE ,L = c(rep(150,3), rep(100,2),rep(150, 4)), ### L argument can be given in mother function
                                                   chrsStyle = "Ensembl", species = "ovis_aries", currentStyle = "Ensembl")
            Gaps_filtered_length <- future_map(IdealMCC_RegionMat,~ makeGRangesFromDataFrame(as.data.frame(gaps(.x$regions)[-c(1,2)])[as.data.frame(gaps(.x$regions)[-c(1,2)])$start > 1,],
                                                                                      keep.extra.columns = TRUE ) ) %>% # To exclude initial gaps
              future_map( ~sort.GenomicRanges(.x)) %>%
              future_map(~ .x[which(width(.x) <= as.integer(MCC_MRG_Comb[c,2])),])
              future_map(seq_along(Gaps_filtered_length), ~tibble(Leader_Region = precede(Gaps_filtered_length[[.x]], IdealMCC_RegionMat[[.x]][["regions"]])) %>% drop_na())
          }
          
          
          nearest_regions_gap_filtered <- nearest_regions_gap_filtered()
          
          
          future_map2(nearest_regions_gap_filtered, future_map(nearest_regions_gap_filtered, ~ tibble(LeaderIndex = is.sequential(as.vector(.x[[1]]))) %>% mutate( FollowerIndex = lag(LeaderIndex + 1, default = 1))), ~  tibble(Leader = ..1[[1]][..2[[1]]]
                                                                                ,Follower = as.integer(..1[[1]][..2[[2]]]-1) ))
        }
        
        
        S_and_E_Regions <- S_and_E_Regions()
        IdealMCC_RegionMat_tibble <- function() future_map(seq_along(IdealMCC_RegionMat), ~ as.tibble(IdealMCC_RegionMat[[.x]]$regions ))
        
        Raw_Gap_filled_IRanges <- function() future_map2(IdealMCC_RegionMat_tibble, S_and_E_Regions, ~ bind_cols(..1[..2[[2]],c("seqnames","start","strand")], ..1[..2[[1]], c("end")]) %>% mutate(width = (end - start)+1) )
        
        gap_filled_regions <- function() future_map(S_and_E_Regions, ~as.integer(levels(as.factor(unlist(map2(.x[[2]], .x[[1]], ~ seq(from = ..1, to = ..2) ))))))
        IdealMCC_RegionMat <- seq(1:27)
        names(IdealMCC_RegionMat) <- as.character(c(seq(1:26),"X"))
        
        future_map(names(IdealMCC_RegionMat),
            ~ list( regions = makeGRangesFromDataFrame( bind_rows(Raw_Gap_filled_IRanges[[which(names(IdealMCC_RegionMat) == as.character(.x))]],
                                                                  IdealMCC_RegionMat_tibble[[which(names(IdealMCC_RegionMat) == as.character(.x))]][-c(gap_filled_regions[[which(names(IdealMCC_RegionMat) == as.character(.x))]]),
                                                                                                                                                    c("seqnames", "start","strand","end","width")]))) )
          
      }
      
      compiled_RegionMat <- compiled_RegionMat()
      names(compiled_RegionMat) <- as.character(c(seq(1:26),"X"))
      future_map(names(compiled_RegionMat), ~ list(regions = compiled_RegionMat[[.x]][["regions"]],
                                                                coverageMatrix = getRegionCoverage(fullCov = fullCov[.x],regions = compiled_RegionMat[[.x]][["regions"]]),
                                                                  bpCoverage = map(seq_len(NROW(compiled_RegionMat[[.x]][["regions"]])),
                                                                                   .f =  function(r) {as.matrix(fullCov[[.x]][as.integer(as_tibble(compiled_RegionMat[[.x]][["regions"]])[r, "start"]): as.integer(as_tibble(compiled_RegionMat[[.x]][["regions"]])[r, "end"]),])} )), .progress = TRUE)
      
     }
     
     Final_RegionMat <- Final_RegionMat()
    names(Final_RegionMat) <-  as.character(c(seq(1:26),"X"))
    for (b in seq_len(length(compiled_RegionMat))){
      seqinfo(compiled_RegionMat[[b]][["regions"]]) <- keepSeqlevels(OarSeqinfo, names(compiled_RegionMat)[b])
    }
    rm(b)
    Deltas[[paste("RegionMats", as.character(MCC_MRG_Comb[c,1]), as.character(MCC_MRG_Comb[c,2]), sep = "_")]] <<-  bind_rows(map(seq_len(length(NonOverlappedPairs)),
                                                                                                                                 ~as_tibble(findOverlapPairs(ranges(NonOverlappedPairs[[..1]]),
                                                                                                                                                             ranges(Final_RegionMat[[..1]][['regions']])))), .id = "Chr") %>%
      transmute( DeltaVal = abs(first.start - second.start) + abs(first.end - second.end)) %>% summarise(IntactDeltaPercent = mean(DeltaVal == 0))
    
    })
    
    Deltas <<- bind_rows(Deltas, .id = "MCC_MRG")
    Deltas[["MCC_MRG"]] <<- as.factor(Deltas[["MCC_MRG"]])
    #Deltas[["Chr"]] <- as.factor(Deltas[["Chr"]])
    Deltas <<- arrange(Deltas, c(MCC_MRG))
    
}

debug(CalculateDeltas)
CalculateDeltas(MCC_MRG_Comb = MCC_MRG_Comb[2:3,],fullCov = fullCov, NonOverlappedPairs = NonOverlappedPairs,SavePath = "/home/animalscience/MPourjam/RegionMats" )





