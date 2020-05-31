					##### Installing and/or loading Packages
source("Packages.R")
					##### Loading Preliminaries
source("Preliminaries.R")
					##### Definfing Functions
source("Functions.R")

					##### Step1- Complining some non-overlapping exons as a milestone

# It is possible and beneficial to define this step as a functioin. However,
# it is time consuming rigth now (12th April 2020).
if (!file.exist(paste0(SavePath, "NonOverlappedPairs.RData"))){
codingexons <- genomicState$codingGenome[mcols(genomicState$codingGenome)[,1] == "exon", ]
codingexonsByChr <- future_map(levels(seqnames(codingexons)), ~ ranges(codingexons[seqnames(codingexons) == .x,]))
names(codingexonsByChr) <-  seqnames(OarSeqinfo)
OverlappedPairs <- future_map(codingexonsByChr, ~ findOverlapPairs(.x,.x))
OverlappedExonsNames <- future_map(seq_len(27), ### Unless a reason could be found, it seems futile finding overlappedExonNames first!!!!
	 ~ as.data.frame(table(as.data.frame(OverlappedPairs[[..1]])[,4]))[as.data.frame(table(as.data.frame(OverlappedPairs[[..1]])[,4]))[,2] > 1,1])
NonOverlappedPairs <- future_map(seq_len(length(codingexonsByChr)),
 ~ codingexonsByChr[[..1]][! as.data.frame(codingexonsByChr[[..1]])[,4] %in% OverlappedExonsNames[[..1]]] )
names(NonOverlappedPairs) <-  seqnames(OarSeqinfo)
save(NonOverlappedPairs, file = file.path(paste0(SavePath, "NonOverlappedPairs.RData")))
}else{
load(file.path(paste0(SavePath, "NonOverlappedPairs.RData")))
}


					##### Step2- Finding overlaps between newregions and NonOverlappedPairs

if(length(list.files(SavePath)))

bind_rows(future_map(seq_len(length(NonOverlappedPairs)), 
	~as_tibble(findOverlapPairs(ranges(NonOverlappedPairs[[..1]]),
	 ranges(Final_RegionMat[[..1]][['regions']])))), .id = "Chr")

	  %>% 
transmute( DeltaVal = abs(first.start - second.start) + abs(first.end - second.end)) %>% summarise(IntactDeltaPercent = mean(DeltaVal == 0))