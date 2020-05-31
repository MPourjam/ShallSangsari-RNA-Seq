##### Finding Nonoverlapping coding exons

# Down are the numbers of overlaps between different exon pairs on each chromosome
codingexons <- genomicState$codingGenome[mcols(genomicState$codingGenome)[,1] == "exon", ]
#Next function used to utilize 'ranges' function which incurrd loss of 'strand' so that the better "granges" was substituted.
codingexonsByChr <- map(levels(seqnames(codingexons)), ~ granges(codingexons[seqnames(codingexons) == .x,])) 
names(codingexonsByChr) <-  as.character(seqnames(OarSeqinfo))
OverlappedPairs <- map(codingexonsByChr, ~ findOverlapPairs(.x,.x))
OverlappedExonsNames <- map(seq_len(27),~ as.data.frame(table(as.data.frame(OverlappedPairs[[..1]])[,4]))[as.data.frame(table(as.data.frame(OverlappedPairs[[..1]])[,4]))[,2] > 1,1])
NonOverlappedPairs <- map(seq_len(length(codingexonsByChr)), ~ codingexonsByChr[[..1]][! as.data.frame(codingexonsByChr[[..1]])[,4] %in% OverlappedExonsNames[[..1]]] )
names(NonOverlappedPairs) <-  seqnames(OarSeqinfo)


# RegionMat exploration

# Setting Delta for ER's alignment on each RegionMat(MCC)
regionMat <- load(file = file.path(paste(SavedData,"SavedData/RegionMats/RegionMats_mean_10.RData", sep = "")))


##  Step5-1 Save path of big dataframes
for (i in seq_len(length(dirNamesPer))){
  
  if (! dirNamesPer[i] %in% list.dirs(path = "/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_V/Results/", full.names = FALSE)) {
    dir.create(path = file.path(paste0("/home/MPourjam/RNA-Seq/RProject", "/","RNA-Seq_Oar3.1_V", "/", "Results","/", dirNamesPer[i])), recursive = TRUE)
  }
  
} 
rm(i)
SavedData <- file.path("/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_V/SavedData/")
if (!"SavedData" %in% list.dirs("/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_V/", full.names = FALSE)) {
  dir.create(SavedData, recursive = TRUE)
}




##  Step2-3 Filter coverage

# Filter coverage
#if (any(grepl("^filteredCov\\.RData$", dir(SavedData)))) {
#  load(paste0(SavedData, "filteredCov.RData", sep = ""))
#} else {
#  filteredCov <- lapply(fullCov, filterData, cutoff = 2, filter =  , TotalMapped = rep(1, length(bamfileslist)), targetSize = 1)
#  save(filteredCov, file = paste0(SavedData, ls()[which(str_detect(ls(),"^filteredCov$"))], ".RData"))
#}

CF_Value_Type <- unique(str_split_fixed(dirNamesPer, "_",4)[,3:4])
CF_Value_Type <- CF_Value_Type[CF_Value_Type[,1] == "mean",] # Make it in case you only need "mean" filter policy # Critical as deciding wether to use 'RegionMat' or 'filteredCov'
if (!"filteredCovs" %in% list.dirs(SavedData, full.names = FALSE)) {
  dir.create(path = paste0(SavedData, "filteredCovs"))
  filteredCovsPath <<- file.path(paste0(SavedData, "filteredCovs"))
} else {
  filteredCovsPath <<- file.path(paste0(SavedData, "filteredCovs"))
}

for (x in seq_len(nrow(CF_Value_Type))){
  if (!any(str_detect(pattern = paste0(paste("filteredCov", CF_Value_Type[x,1], CF_Value_Type[x,2], sep = "_"),".RData"),string = list.files(filteredCovsPath, full.names = FALSE)))) {
    assign(paste("filteredCov", CF_Value_Type[x,1], CF_Value_Type[x,2], sep = "_"), 
           lapply(fullCov, filterData, cutoff = as.numeric(CF_Value_Type[x,2]), filter = CF_Value_Type[x,1] , TotalMapped = rep(1, length(bamfileslist)), targetSize = 1))
    file.create(file.path(paste0(paste(filteredCovsPath, paste("filteredCov", CF_Value_Type[x,1], CF_Value_Type[x,2], sep = "_"),sep = "/"), ".RData")))
    save(list = ls()[which(str_detect(ls(), paste("filteredCov", CF_Value_Type[x,1], CF_Value_Type[x,2], sep = "_")))],
         file = file.path(paste0(paste(filteredCovsPath, paste("filteredCov", CF_Value_Type[x,1], CF_Value_Type[x,2], sep = "_"),sep = "/"), ".RData")))
    rm(list = ls()[which(str_detect(ls(), paste("filteredCov", CF_Value_Type[x,1], CF_Value_Type[x,2], sep = "_")))]) 
    print(paste0("Filtered Coverage data of cutoff value ", CF_Value_Type[x,2], " and (", CF_Value_Type[x,1] , ") filter policy has been saved."))
  }
}
rm(x)

str_glue_data(.x = map(seq_len(nrow(xc)), ~xc[.x,]))

map(seq_len(nrow(xc)), ~ as.character(xc[.x,]))


map(ncol(xc), ~list(as.character(xc[,.x])))

apply(xc, 1, str_flatten, collapse = "_")





for (b in names(models)) {print(b)}
