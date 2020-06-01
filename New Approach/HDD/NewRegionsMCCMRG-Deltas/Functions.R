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

# This is the core function
NewRegions <- function(MCC_MRG_Grid, fullCov_Path, DBDir_Path, CalclulateDelta = FALSE ){  
  stopifnot(!file.exist(file.path(paste0(DBDir_Path,"RegionMat_Path.RData"))))
  load( file = file.path(paste0(DBDir_Path,"RegionMat_Path.RData"))) ### This will load the variable "RegionMat_Path" containing the path to RegionMat_MCC file
  stopifnot(file.exist(file.path(paste0(DBDir_Path,"MRG.RData"))))
  load( file = file.path(paste0(DBDir_Path,"MRG.RData")))### This will load the variable "MRG"        
  MRG_RegionSet <- vector(mode = 'list',length = 2)
  stopifnot(file.exist(fullCov_Path))
  fullCov <- load(file = file.path(fullCov_Path))
  Deltas_file_path <- file.path(paste0(DBDir_Path,"Deltas.RData"))
  SaveDir <- paste0(dirname(RegionMat_Path),"/")
  RegionMat <- load(file = file.path(RegionMat_Path))
  .subset2(MRG_RegionSet,1) <- future_map(RegionMat, ~ .subset2(.x, 1)) #### The Structure of RegionMat is RegionMat$Chr$regions
  .subset2(MRG_RegionSet,2) <- future_map(.subset2(MRG_RegionSet, 1), ~ gaps_Grange(.x))
  
  
    gc()

if (!file.exist(file.path(paste0(SaveDir, paste(basename(SaveDir), as.character(MRG), sep = "_" ), ".RData" )) )) {
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
      ~strand(..2[.subset2(..1, 2),]))
    New_Regions_Ranges <- future_map2(Leader_Follower_PrecedingRegion, .subset2(MRG_RegionSet, 1), 
      ~IRanges(start = start(..2[[.subset2(..1, 2)]]), end = end(..2[.subset2(..1, 1),])))
    New_Regions_GRange <- pmap(New_Regions_Seqnames, New_Regions_Strand, New_Regions_Ranges,
      ~ GRanges(seqnames = ..1, ranges = ..3, strand = ..2))

    ### Removing compiled regions from the list of regions
  excluding_regions_allChrs <- future_map(Leader_Follower_PrecedingRegion, ~ exclude_compiled_regions)

  compiled_regions_excluded_tibble <- future_map2(.subset2(MRG_RegionSet,1) , excluding_regions_allChrs,
   ~ as_tibble(..1[-c(..2),c("seqnames","ranges","strand")])) ### Might incur problems due to nature of regions as a Grange (not subsettable)
  RegionMat_MCC_MRG <- future_map2(compiled_regions_excluded_tibble, New_Regions_GRange, ### fitering for regions longer than 3 to exclude microexons
   ~ list(regions = sort.GenomicRanges(makeGRangesFromDataFrame(dplyr::filter(bind_rows(..1,..2), width(ranges) > 3) ) %>% 'Seqinfo<-'(OarSeqinfo)), #### Retrieving regionCoverage for 
    bpCoverage = getRegionCoverage(fullCov, sort.GenomicRanges(makeGRangesFromDataFrame(dplyr::filter(bind_rows(..1,..2), width(ranges) > 3) ))) ) )
  #### Doing all at once may put such a heavy burden on the RAM (pay attention to above)
  names(RegionMat_MCC_MRG) <- seqnames(OarSeqinfo)

  # Saving NewRegions
save(RegionMat_MCC_MRG,
  file = file.path(paste0(SaveDir, paste(basename(SaveDir), as.character(MRG), sep = "_" ), ".RData" )))
}

### Updating MRG.RData and RegionMat_Path.RData
MCC <- as.numeric(stringr::str_remove(stringr::str_extract(string = basename(RegionMat_Path),
 pattern = "[^[[:alpha:]]]+"), ".$"))

 ### Calcualting Deltas
if (CalclulateDelta){
  if (!c("NonOverlappedExons") %in% ls() ){
  warning(paste0("NonOverlappedExons was not loaded trying to load it from the", as.character(DBDir_Path),"!"),
   call. = FALSE, immediate. = TRUE)
  if (!file.exist(paste0(DBDir_Path,"NonOverlappedExons.RData") )) {
    warning(paste0("NonOverlappedExons.RData does not exist in the following directory:",
     as.character(DBDir_Path)) ,call. = FALSE, immediate. = TRUE)
  } else {
    load(paste0(DBDir_Path,"NonOverlappedExons.RData"))
  }
} else {
    if (!file.exist(Deltas_file_path)) {
      warning(paste0("Deltas dataframe file was not found in the following directory: ",
       as.character(DBDir_Path), ". Attemping to create one."), immediate. = TRUE)
      Deltas_names <- tidyr::unite(MCC_MRG_Grid, "MCC_MRG", sep= "_", remove = TRUE )
      Deltas <- vector(mode = "list", length = nrow(MCC_MRG_Grid))
      names(Deltas) <- as.character(Deltas_names[["MCC_MRG"]])

      ## Filtering for regions loonger than 3bp-long has been done on line 80 (RegionMat_MCC_MRG)
      LaidPairs <- future_map(seq_len(length(NonOverlappedExons)),
        ~findOverlapPairs(granges(RegionMat_MCC_MRG[[.x]][['regions']]),
          granges(NonOverlappedExons[[.x]])), .id = "Chr")
      ### Filtering Er's laid on multiple Exons
      Ers_Laid_One_Exon <- future_map(LaidPairs,~ unite(data.frame(.x@first), col = "Region", sep="-", remove = TRUE) %>% table() %>% as_tibble() %>%
          dplyr::filter(n == 1) %>% tidyr::separate(col = ".", into = c("seqnames", "start", "end", "width", "strand"),sep = "-", remove = TRUE) %>% dplyr::select(-c(n)) )
     Ers_Laid_One_Exon <- map(Ers_Laid_One_Exon, ~ base::`$<-`(.x, "seqnames", as.factor(.x[[1]]) ))
     Ers_Laid_One_Exon <- map(Ers_Laid_One_Exon, ~ base::`$<-`(.x, "start", as.integer(.x[["start"]]) ))
     Ers_Laid_One_Exon <- map(Ers_Laid_One_Exon, ~ base::`$<-`(.x, "end", as.integer(.x[["end"]]) ))
     Ers_Laid_One_Exon <- map(Ers_Laid_One_Exon, ~ base::`$<-`(.x, "width", as.integer(.x[["width"]]) ))
     Ers_Laid_One_Exon <- map(Ers_Laid_One_Exon, ~ base::`$<-`(.x, "strand", as.factor(.x[["strand"]]) ))
      
      LaidPairs <- future_map2(LaidPairs, Ers_Laid_One_Exon, ~ dplyr::semi_join(x = as_tibble(..1@first), y = ..2))
      ############################################### Yet to be done ####################################################
      ### We need to prune the ERs laying on multiple exons.Temporarily, we omit these steps due to uncertainty regarding
      ###  applying width() function on a output of findOverlapPairs().
      ###################################################################################################################
      Deltas[[paste0(as.character(MCC),"_",as.character(MRG))]] <-
     # %>% transmute( DeltaVal = abs(first.start - second.start) + abs(first.end - second.end))
      
      save(Deltas, file = Deltas_file_path)
    } else {
      load(Deltas_file_path)
      Deltas[[paste0(as.character(MCC),"_",as.character(MRG))]] <- future_map_dfr(seq_len(length(NonOverlappedExons)),
       ~as_tibble(findOverlapPairs(ranges(NonOverlappedExons[[..1]]), # We may need to change the query and subject 
        ranges(RegionMat_MCC_MRG[[..1]][['regions']]))), .id = "Chr") %>% 
      transmute( DeltaVal = abs(first.start - second.start) + abs(first.end - second.end))
      save(Deltas, file = Deltas_file_path)
    }
  }
}

NextMCC_MRG <- MCC_MRG_Grid[which(MCC_MRG_Grid[,2] == MRG & MCC_MRG_Grid[1,] == MCC) + 1 , ]

if (is.na(NextMCC_MRG[1,2])) {
  stop(call. = FALSE, "Out of MCC_MRG_Grid range!")
} else {
  MRG <- as.numeric(NextMCC_MRG[1,2])
 save( MRG , file = file.path(paste0(DBDir_Path,"MRG.RData")))
 
 RegionMat_Path <- stringr::str_subset(list.files(RegionMatsPath, recursive=TRUE, full.names=TRUE),
    pattern = paste0("RegionMats_",as.numeric(NextMCC_MRG[1,1]),".RData"))

 if (length(RegionMat_Path) != 1) { 
 #### Because MCC iteration is slower than MRG
    stop(paste0(" The initial ", paste0("RegionMats_",as.numeric(NextMCC_MRG[1,1]),".RData"),
      " was not found or appears more than once!!!"), call. = FALSE)
  }else{
  save(RegionMat_Path, file = file.path(paste0(DBDir_Path,"RegionMat_Path.RData")))
  }
  
  }
}



