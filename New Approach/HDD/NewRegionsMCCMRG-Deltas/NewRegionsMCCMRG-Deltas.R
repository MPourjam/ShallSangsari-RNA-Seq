          ##### Installing and/or loading Packages #####
source("Packages.R")
					#####	Setting preliminary variables	#####
source("Preliminaries.R")
          ### Functions
source("Functions.R")

		### Step1- Creating the RegionsMat's folder

if (!"RegionMats" %in% list.dirs(SavePath, full.names = FALSE)) {
  dir.create(path = paste0(SavePath, "RegionMats"))
  RegionMatsPath <- file.path(paste0(SavePath, "RegionMats/"))
} else {
  RegionMatsPath <- file.path(paste0(SavePath, "RegionMats/"))
}


		### Step2- Creating the RegionMats


for (c in seq_along(MCC)){   ##### Create seperate dirs for each RegionMat_MCC to store RegionMat_MCC.RData and MRG derived regions.RData 
  if (!any(str_detect(pattern = paste0(paste("RegionMats", MCC[c], sep = "_"),".RData"),
    string = list.files(RegionMatsPath, full.names = FALSE, recursive=TRUE)))) {

    file.create(file.path(paste0(paste0(paste0(RegionMatsPath,
     paste("RegionMats", MCC[c], sep = "_"),"/"), paste("RegionMats", MCC[c], sep = "_")), ".RData")))

    New_RegionMat <- regionMatrix(fullCov, cutoff = MCC[c], runfilter = FALSE ,L = c(rep(150,3), rep(100,2),rep(150, 4)),
                        chrsStyle = "Ensembl", species = "ovis_aries", currentStyle = "Ensembl", returnBP = FALSE)
    
    save(New_RegionMat,
         file = file.path(paste0(paste0(paste0(RegionMatsPath,
     paste("RegionMats", MCC[c], sep = "_"),"/"), paste("RegionMats", MCC[c], sep = "_")), ".RData")))
    rm(New_RegionMat)
    gc()
    print(paste0("RegionMats data of cutoff value: ", MCC[c], " and mean filter policy has been saved to ",
      file.path(paste0(paste0(RegionMatsPath, paste("RegionMats", MCC[c], sep = "_")), ".RData"))))
  }
}
rm(c)


		### Step3- Initializing  MRG.RData and RegionMat_Path.RData

MRGfile_Path <- file.path(paste0(SavePath,"MRG.RData"))
RegionMat_Path_file <- file.path(paste0(SavePath,"RegionMat_Path.RData"))

if (!file.exist(MRGfile_Path)){
  file.create(MRGfile_Path)
  MRG <- as.numeric(MCC_MRG_Grid[1,2])
  save(MRG, file = MRGfile_Path)
} 
if (!file.exist(RegionMat_Path_file)){
  file.create(RegionMat_Path_file)
  ########### NEED THE PATH TO RegionMat Matrice
  RegionMat_Path <- stringr::str_subset(list.files(RegionMatsPath, recursive=TRUE, full.names=TRUE),
    pattern = paste0("RegionMats_",as.numeric(MCC_MRG_Grid[1,1]),".RData"))
  if (length(RegionMat_Path) == 0){
    stop(paste0(" The initial ", paste0("RegionMats_",as.numeric(MCC_MRG_Grid[1,1]),".RData"),
      " was not found or appears more than once!!!"), call. = FALSE)
  }else{
  save(RegionMat_Path, file = RegionMat_Path_file)  
   }
}
    ### Step4- Creating the NonOverlappedExons

if (!file.exist(paste0(SavePath, "NonOverlappedExons.RData"))){
codingexons <- genomicState$codingGenome[mcols(genomicState$codingGenome)[,1] == "exon", ]
#Next function used to utilize 'ranges' function which incurrd loss of 'strand' so that the better "granges" was substituted.
codingexonsByChr <- future_map(levels(seqnames(codingexons)), ~ granges(codingexons[seqnames(codingexons) == .x,]))
names(codingexonsByChr) <-  seqnames(OarSeqinfo)
OverlappedPairs <- future_map(codingexonsByChr, ~ findOverlapPairs(.x,.x))
OverlappedExonsNames <- future_map(seq_len(27), ### Unless a reason could be found, it seems futile finding overlappedExonNames first!!!!
   ~ as.data.frame(table(as.data.frame(OverlappedPairs[[..1]])[,4]))[as.data.frame(table(as.data.frame(OverlappedPairs[[..1]])[,4]))[,2] > 1,1])
NonOverlappedExons <- future_map(seq_len(length(codingexonsByChr)),
 ~ codingexonsByChr[[..1]][! as.data.frame(codingexonsByChr[[..1]])[,4] %in% OverlappedExonsNames[[..1]]] )
names(NonOverlappedExons) <-  as.character(seqnames(OarSeqinfo))
save(NonOverlappedExons, file = file.path(paste0(SavePath, "NonOverlappedExons.RData")))
}else{
load(file.path(paste0(SavePath, "NonOverlappedExons.RData")))
}


    ### Step5- Performing the core function, creating the new regions and Deltas

# Examine the load and save behaviour
NewRegions( MCC_MRG_Grid = MCC_MRG_Grid, DBDir_Path = SavePath, CalclulateDelta = TRUE)
