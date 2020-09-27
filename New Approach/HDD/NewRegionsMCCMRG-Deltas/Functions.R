### Functions
is.sequential <- function(vector, diff = 1){
  i <- seq_len(length(vector))
  gap <- numeric()
  for (iter in i[-length(i)]){
    if (vector[iter+ 1] == vector[iter] + diff){
      # it will go to upper iteration
    } else{
    gap <- c(gap, iter)
    }
  }
 gap <- c(gap, i[length(i)])
}


gaps_Grange <- function(regionMat_region_Grange) {
  gaps <- as.data.frame(gaps(regionMat_region_Grange))
  gaps <- makeGRangesFromDataFrame(gaps[gaps$start > 1,], keep.extra.columns = TRUE )
  gaps <-  sort.GenomicRanges(gaps)
  gaps
}

exclude_compiled_regions <- function(PrecedingRegion_each_Chr){  
rows <- nrow(PrecedingRegion_each_Chr)
exluding_regions <- unique(unlist(map(seq_len(rows), ~ seq(as.numeric(PrecedingRegion_each_Chr[.x,"FollowerRegion"]) ,
 as.numeric(PrecedingRegion_each_Chr[.x,"LeaderRegion"]) ))))
exluding_regions
}

# This is the core function
NewRegions <- function(MCC_MRG_Grid, fullCov_Path, DBDir_Path, currentSample, bamfileslist, CalclulateDelta = FALSE ) {
  
  SavePath_Sample <- paste0(SavePath, currentSample,"/")
  stopifnot(file.exists(file.path(paste0(SavePath_Sample,"RegionMat_Path.RData")))) 
  load( file = file.path(paste0(SavePath_Sample,"RegionMat_Path.RData"))) ### This will load the variable "RegionMat_Path" containing the path to RegionMat_MCC file
  stopifnot(file.exists(file.path(paste0(SavePath_Sample,"MRG.RData"))))
  load( file = file.path(paste0(SavePath_Sample,"MRG.RData")))### This will load the variable "MRG"        
  MRG_RegionSet <- vector(mode = 'list',length = 2)
  if (!"fullCov" %in% ls(envir = globalenv()) ){
    stopifnot(file.exists(fullCov_Path))
    load(file = file.path(fullCov_Path))
  }
  
  bamfileslist_1Samp_Path <- file.path(paste0(DBDir_Path,"bamfilelist.RData"))
  stopifnot(file.exists(file.path(paste0(DBDir_Path,"bamfilelist.RData"))))
  load(bamfileslist_1Samp_Path)
  Deltas_file_path <- file.path(paste0(DBDir_Path,"Deltas.RData"))
  SaveDir <- paste0(dirname(RegionMat_Path),"/")
  load(file = file.path(RegionMat_Path)) ### This loads RegionMat (13.72 Gb) to the memory
  MRG_RegionSet[[1]] <- map(RegionMat, ~ base::.subset2(.x, 1)) #### The Structure of RegionMat is RegionMat$Chr$regions
  MRG_RegionSet[[2]] <- map(base::.subset2(MRG_RegionSet, 1), ~ gaps_Grange(.x))
  rm("RegionMat")
  #gc(verbose=T)
    gc()

if ( !file.exists(file.path(paste0(SaveDir, paste(basename(SaveDir), as.character(MRG), sep = "_" ), ".RData" ))) ) {
    ## Filtering gaps
    gap_filtered <- map(.subset2(MRG_RegionSet,2), ~ .x[which(width(.x) <= as.integer(MRG)),])

    ## Finding receding regions
    preceding_region <- map2(gap_filtered, .subset2(MRG_RegionSet,1), ~ tibble(Leader_Region = precede(..1, ..2 )) %>% drop_na() )
    # preceding_region is a list of 27(number of chromosomes) elements

    ### 
    Leader_Follower_Index <- map(preceding_region, ~ tibble(LeaderIndex = is.sequential(.subset2(.x,1))) %>%
     mutate(FollowerIndex = dplyr::lag(LeaderIndex + 1, default = 1)) )  #??? What about the last one in preceding_region

    ### Retrieving corresponding regoins using indices from preceding regions vector (list of 27)
   Leader_Follower_PrecedingRegion <- map(names(preceding_region), ~tibble(LeaderRegion = as_tibble(unlist(ifelse(any(preceding_region[[.x]][.subset2(Leader_Follower_Index[[.x]], 1),]+1 > IRanges::NROW(.subset2(MRG_RegionSet, 1)[[.x]])),
                                                                                                               rbind.data.frame(as_tibble(preceding_region[[.x]][.subset2(Leader_Follower_Index[[.x]], 1)[-c(length(.subset2(Leader_Follower_Index[[.x]], 1)))],]+1),
                                                                                                                                preceding_region[[.x]][.subset2(Leader_Follower_Index[[.x]], 1)[c(length(.subset2(Leader_Follower_Index[[.x]], 1)))],]),
                                                                                                     preceding_region[[.x]][.subset2(Leader_Follower_Index[[.x]], 1),]+1)))[[1]] , FollowerRegion = preceding_region[[.x]][.subset2(Leader_Follower_Index[[.x]],2),][[1]]))
    names(Leader_Follower_PrecedingRegion) <- seqnames(OarSeqinfo)

    ### Constituting new regioins
    New_Regions_Seqnames <- map2(Leader_Follower_PrecedingRegion, .subset2(MRG_RegionSet, 1),
     ~ seqnames(..2[.subset2(..1, 2),]))
    New_Regions_Strand <- map2(Leader_Follower_PrecedingRegion, .subset2(MRG_RegionSet, 1),
      ~strand(..2[.subset2(..1, 2),]))
    New_Regions_Ranges <- map2(Leader_Follower_PrecedingRegion, .subset2(MRG_RegionSet, 1), 
     ~IRanges(start = start(..2[.subset2(..1, 2),]), end = end(..2[.subset2(..1, 1),])))
    New_Regions_GRange <- map(seqnames(OarSeqinfo), ~GRanges(seqnames = New_Regions_Seqnames[[.x]],
                                                             ranges = New_Regions_Ranges[[.x]],
                                                             strand = New_Regions_Strand[[.x]]))
    names(New_Regions_GRange) <- seqnames(OarSeqinfo)

    ### Removing compiled regions from the list of regions
  excluding_regions_allChrs <- map(Leader_Follower_PrecedingRegion, ~ exclude_compiled_regions(.x))

  compiled_regions_excluded_tibble <- map2(.subset2(MRG_RegionSet,1) , excluding_regions_allChrs,
     ~ as_tibble(..1[-c(..2),])[,c("seqnames","start","end","width","strand")]) ### Might incur problems due to nature of regions as a Grange (not subsettable)
  
  #RegionMat_MCC_MRG <- map2(compiled_regions_excluded_tibble, New_Regions_GRange, ### fitering for regions longer than 3 to exclude microexons
  # ~ list(regions = sort.GenomicRanges(makeGRangesFromDataFrame(dplyr::filter(bind_rows(..1,..2), width(ranges) > 3) ) %>% 'Seqinfo<-'(OarSeqinfo)), #### Retrieving regionCoverage for 
  #  bpCoverage = getRegionCoverage(fullCov, sort.GenomicRanges(makeGRangesFromDataFrame(dply::filter(bind_rows(..1,..2), width(ranges) > 3) ))) ) )
  
 RegionMat_MCC_MRG <- makeGRangesFromDataFrame( map2_dfr(compiled_regions_excluded_tibble, New_Regions_GRange, ### It Seems 'makeGRangesListFromDataFrame' needs to be changed to 'makeGRangesFromDataFrame'
                                                              ~ as_tibble(sort.GenomicRanges(makeGRangesFromDataFrame(dplyr::filter(bind_rows(..1,as_tibble(..2)), width > 3) )))))
  
  seqinfo(RegionMat_MCC_MRG) <- OarSeqinfo
  #### Doing all at once may put such a heavy burden on the RAM (pay attention to above)
  
  # Re-creating RegionMat with new set of regions
    # Regions
      ReformedRegionMat_regions <- vector("list", 27)
      names(ReformedRegionMat_regions) <- c(rep(1:26), "X")
      for (i in seq(1,27)){
        ReformedRegionMat_regions[[i]]$regions <- split(RegionMat_MCC_MRG, factor(RegionMat_MCC_MRG@seqnames))[[i]]
      } 
      rm(i)
      RegionMat_MCC_MRG <- ReformedRegionMat_regions
      gc()
    # CoverageMatrix
      # They are not needed now.
    
    # bpCoverage
     # They are not needed now.
      
      
  # Saving NewRegions
save(RegionMat_MCC_MRG,
  file = file.path(paste0(SaveDir, paste(basename(SaveDir), as.character(MRG), sep = "_" ), ".RData" )))
}

### Updating MRG.RData and RegionMat_Path.RData
MCC <- as.numeric(stringr::str_remove(stringr::str_remove(stringr::str_extract(string = basename(RegionMat_Path),
 pattern = "[^[[:alpha:]]]+"), ".$"), "^_"))

 ### Calcualting Deltas
if (CalclulateDelta){
  if (!c("NonOverlappedExons") %in% ls() ){
  warning(paste0("NonOverlappedExons was not loaded trying to load it from the", as.character(DBDir_Path),"!"),
   call. = FALSE, immediate. = TRUE)
  stopifnot(file.exists(paste0(DBDir_Path,"NonOverlappedExons.RData")) )
  # if (!file.exists(paste0(DBDir_Path,"NonOverlappedExons.RData") )) {
  #   warning(paste0("NonOverlappedExons.RData does not exist in the following directory:",
  #    as.character(DBDir_Path)) ,call. = FALSE, immediate. = TRUE)
  # } else {
    load(paste0(DBDir_Path,"NonOverlappedExons.RData"))
  # }
}
# } else {
    if (!file.exists(Deltas_file_path)) {
      warning(paste0("Deltas dataframe file was not found in the following directory: ",
       as.character(DBDir_Path), ". Attemping to create one."), immediate. = TRUE)
      Deltas_names <- tidyr::unite(MCC_MRG_Grid, "MCC_MRG_Sample", sep= "_", remove = TRUE )
      Deltas <- vector(mode = "list", length = nrow(MCC_MRG_Grid))
      names(Deltas) <- as.character(Deltas_names[["MCC_MRG_Sample"]])

      ## Filtering for regions loonger than 3bp-long has been done on line 80 (RegionMat_MCC_MRG)
      LaidPairs <- map(seq_len(length(NonOverlappedExons)),
        ~findOverlapPairs(granges(RegionMat_MCC_MRG[[.x]][['regions']]),
          granges(NonOverlappedExons[[.x]])), .id = "Chr")
      ### Filtering Er's laid on multiple Exons
      Ers_Laid_One_Exon <- map(LaidPairs,~ unite(data.frame(.x@first), col = "Region", sep="-", remove = TRUE) %>% table() %>% as_tibble() %>%
          dplyr::filter(n == 1) %>% tidyr::separate(col = ".", into = c("seqnames", "start", "end", "width", "strand"),sep = "-", remove = TRUE) %>% dplyr::select(-c(n)) )
     Ers_Laid_One_Exon <- map(Ers_Laid_One_Exon, ~ base::`$<-`(.x, "seqnames", as.factor(.x[[1]]) ))
     Ers_Laid_One_Exon <- map(Ers_Laid_One_Exon, ~ base::`$<-`(.x, "start", as.integer(.x[["start"]]) ))
     Ers_Laid_One_Exon <- map(Ers_Laid_One_Exon, ~ base::`$<-`(.x, "end", as.integer(.x[["end"]]) ))
     Ers_Laid_One_Exon <- map(Ers_Laid_One_Exon, ~ base::`$<-`(.x, "width", as.integer(.x[["width"]]) ))
     Ers_Laid_One_Exon <- map(Ers_Laid_One_Exon, ~ base::`$<-`(.x, "strand", as.factor(.x[["strand"]]) ))
     names(Ers_Laid_One_Exon) <- OarSeqinfo@seqnames 
     gc()
     LaidPairs <- map(seq_len(length(NonOverlappedExons)),
        ~findOverlapPairs(makeGRangesFromDataFrame(Ers_Laid_One_Exon[[.x]]), granges(NonOverlappedExons[[.x]])), .id = "Chr")
      names(LaidPairs) <- OarSeqinfo@seqnames 
     gc()
     
      Deltas[[paste0(as.character(MCC),"_", as.character(MRG), "_", currentSample)]] <- map_dfr(LaidPairs,
       ~ tibble(Delta = abs(as_tibble(.x@first)$start - as_tibble(.x@second)$start) + abs(as_tibble(.x@first)$end - as_tibble(.x@second)$end)),.id = "Chr" )
     # %>% transmute( DeltaVal = abs(first.start - second.start) + abs(first.end - second.end))
      
      save(Deltas, file = Deltas_file_path)
    } else {
      load(Deltas_file_path)
      load(file.path(paste0(SaveDir, paste(basename(SaveDir), as.character(MRG), sep = "_" ), ".RData" ))) ## This will load "RegionMat_MCC_MRG"
      LaidPairs <- map(seq_len(length(NonOverlappedExons)),
        ~findOverlapPairs(granges(RegionMat_MCC_MRG[[.x]][['regions']]),
          granges(NonOverlappedExons[[.x]])), .id = "Chr")
      Deltas[[paste0(as.character(MCC),"_", as.character(MRG))]] <- map_dfr(LaidPairs,
       ~ tibble(Delta = abs(as_tibble(.x@first)$start - as_tibble(.x@second)$start) + abs(as_tibble(.x@first)$end - as_tibble(.x@second)$end)),.id = "Chr" )
      save(Deltas, file = Deltas_file_path)
    }
  # }
}

NextMCC_MRG <- MCC_MRG_Grid[which(MCC_MRG_Grid[[2]] == MRG & MCC_MRG_Grid[[1]] == MCC & MCC_MRG_Grid[[3]] == currentSample) + 1 , ]

if (is.na(NextMCC_MRG[1,2])) {
  stop(call. = FALSE, "Out of MCC_MRG_Grid range!")
} else {
  MRG <- as.numeric(NextMCC_MRG[1,2])
 save( MRG , file = file.path(paste0(DBDir_Path,"MRG.RData")))
 
 RegionMat_Path <- stringr::str_subset(list.files(RegionMatsPath, recursive=TRUE, full.names=TRUE),
    pattern = paste0("RegionMats_",as.numeric(NextMCC_MRG[1,1]),".RData"))
  
bamfileslist_1Samp <- bamfileslist[[which(names(bamfileslist) == NextMCC_MRG[[3]])]]
  save(bamfileslist_1Samp,bamfileslist_1Samp_Path) 
 if (length(RegionMat_Path) != 1) { 
 #### Because MCC iteration is slower than MRG
    stop(paste0(" The initial ", paste0("RegionMats_",as.numeric(NextMCC_MRG[1,1]),".RData"),
      " was not found or appears more than once!!!"), call. = FALSE)
  } else{
  save(RegionMat_Path, file = file.path(paste0(DBDir_Path,"RegionMat_Path.RData")))
    }
  }

}



