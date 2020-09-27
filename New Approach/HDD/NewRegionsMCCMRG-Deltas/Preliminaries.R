          ##### Setting preliminary variables #####


### Adjuctments
plan("multisession")
options(future.globals.maxSize= 8912896000)
options(future.globals.onMissing= "ignore")


### BAM files
bamfilespath <- list.files(path = paste0(cwd,"/DATA/Analyze/Out/Sheep") #### Why do we need bamfilespath...
                            ## in cwd might be anywhere pay attention to addess of bamfilespath
                            , pattern = ".Aligned.out.bam", all.files = TRUE, full.names = TRUE
                            , recursive = TRUE )

readsnames <- gsub(".*/(.*)_Aligned.out.bam", "\\1", bamfilespath)
names(bamfilespath) <- readsnames

Coordbamfilespath <- list.files(path = paste0(cwd,"/DATA/Analyze/Out/Sheep")
                                 , pattern = "sortedByCoord.out.bam$", full.names = TRUE
                                , recursive = TRUE)

Coordbamfilespathindex <- paste0(Coordbamfilespath, ".bai")

Bamfiles <- map(seq_along(Coordbamfilespath), function(i) {BamFile(file = Coordbamfilespath[i],
index = Coordbamfilespathindex[i])})

bamfileslist <-  BamFileList(Bamfiles)
names(bamfileslist) <- readsnames


### Path
cwd <- getwd()
SavePath <- paste0(cwd,"/DB/") ### This the root of our database ### The terminal '/' is necessary ### R session must be initiated in the folder containing /DB folder
SavePath_Sample_FilePath <- file.path(paste0(SavePath,"SavePath_Sample.RData")) ### Path to the RData file which saves tha root path for each samples sets of ERs

if (!file.exists(SavePath_Sample_FilePath)){
SavePath_Sample <- paste0(SavePath, readsnames[1],"/")
save(SavePath_Sample, file = SavePath_Sample_FilePath)
}

bamfileslist_1Samp_Path <- file.path(SavePath,"bamfilelist.RData")
if (!file.exists(bamfileslist_1Samp_Path)){
file.create(bamfileslist_1Samp_Path)
bamfileslist_1Samp <- bamfileslist[[1]]
save(bamfileslist_1Samp, file = bamfileslist_1Samp_Path)
} else {
load(bamfileslist_1Samp_Path)
}

fullCov_Path <- file.path(paste0(SavePath_Sample,"fullCov.RData"))
if (!file.exists(fullCov_Path)){
  
  ### Totalmapped
  
  TotalMapped <- getTotalMapped(bamfileslist_1Samp$path)
  
  ### fullCov 

  fullCov <- fullCoverage(bamfileslist, chrs = c(rep(1:26), "X"), totalMapped = TotalMapped, returnMean = TRUE, runFilter = FALSE) ### leaving targetsize to default
  save(fullCov, file=fullCov_Path)
}

### Variables
MCC_Set <- c(seq(0.5,20,0.5))
MRG_Set <- c(seq(10,100, 10))
MCC_MRG_Grid <- expand.grid(MCC_Set,MRG_Set, readsnames)
if (length(fullCov_Path) == 1){
if (str_detect(fullCov_Path,".*.rds$")  & !"fullCov" %in% ls(envir = globalenv()) ){  ###### What if there are more than one file containing "fullcov" and ".rds
  load(fullCov_Path) #### How is it possible to have fullCov data in a directory which is newly built (SavedData)
} else if(stringr::str_detect(fullCov_Path, ".*.RData")  & !"fullCov" %in% ls(envir = globalenv()) ) {
 load(fullCov_Path)
}
} else {
  stop("There are more than one fulCov files in our working directory(cwd) or niehter rds nor RData files were found." )
}

OarSeqinfo <- Seqinfo(seqnames = c(1:26, "X"),
 seqlengths = c( 275612895, 248993846, 224283230, 119255633, 107901688, 117031472, 100079507,
  90695168, 94726778, 86447213, 62248096, 79100223, 83079144, 62722625, 80923592, 71719816, 
  72286588, 68604602, 60464314, 51176841, 50073674, 50832532, 62330649, 42034648, 45367442,
  44077779, 135437088),isCircular = rep(FALSE, 27), genome = "Oar3.1")

ReadLength <- c(rep(150,3), rep(100,2),rep(150, 4))

Fecuntrt <- c(c(rep("LF",3), "HF",rep("LF",2), rep("HF",2),"LF"))
Breedtrt <- c(c(rep("San",3), rep("Shall",6)))
ThreeGroup <- c(rep("San", 3), rep("Shall_OUT",2), rep("Shall",3), "Shall_OUT")
pheno <- as.data.frame( cbind(Breed = Breedtrt, Fecun = Fecuntrt,
 ThreeGroup = ThreeGroup, ReadLength = ReadLength))
pheno <-  map_dfc(colnames(pheno), ~ as.factor(pheno[[as.character(.x)]]))
colnames(pheno) <- c("Breedtrt","Fecuntrt","ThreeGroup","ReadLength")
rownames(pheno) <- NULL

Breeds <- c("Shall","San")

# pheno$Samples <- names(bamfileslist)

# # Creating Metadata columns and TxDb Object
# met <- data.frame( IllSeqSheep = c(paste("LF_San_", rep(1:3,1), sep = "" )
#                                    , c("HF_Shall_1", "LF_Shall_1","LF_Shall_2",
#                                     "HF_Shall_2","HF_Shall_3", "LF_Shall_3")))
# metadata(bamfileslist) <- met

if (!file.exists(paste0(SavePath, "txdb.RData"))){
  txdb <- TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1 <- makeTxDbPackageFromBiomart(taxonomyId = 9940, port = 80, biomart = "ENSEMBL_MART_ENSEMBL"
                                                                                       , host = "www.ensembl.org", dataset = "oaries_gene_ensembl"
                                                                                       , filter = NULL, circ_seqs = "MT", id_prefix="ensembl_", miRBaseBuild = NA
                                                                                       , version = "1.0"
                                                                                       , maintainer = "Mohsen Pourjam <Pourjam.cs@hotmail.com>"
                                                                                       , author = "Mohsen Pourjam")
  install.packages("TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1",
                   repos = NULL, type = "source") #### Possible Errors due to the path
  library( "TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1")
  keepSeqlevels(TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1::TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1,c(1:26, "X"))
  save(txdb, file = paste0(SavePath, "txdb.RData"))
} else {
  keepSeqlevels(TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1::TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1,c(1:26, "X"))
  load(file = paste0(SavePath, "txdb.RData"))
}

if (file.exists(file.path(paste0(SavePath,"genomicState.RData")))){
  load(file = paste0(SavePath,"genomicState.RData"))
} else {
  genomicState <- makeGenomicState(txdb =  TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1::TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1
                                   , chrs = OarSeqinfo@seqnames
                                   , style = getOption("chrsStyle", "Ensembl"), currentStyle = "Ensembl")
  seqinfo(genomicState) <- OarSeqinfo
  save(object = genomicState,  file = paste0(SavePath,"genomicState.RData"))
}
