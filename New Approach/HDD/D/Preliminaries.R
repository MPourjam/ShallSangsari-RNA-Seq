          ##### Setting preliminary variables #####


### Adjuctments
plan("multisession")
options(future.globals.maxSize= 8912896000)
options(future.globals.onMissing= "ignore")


### Path
cwd <- getwd()
SavePath <- paste0(cwd,"/DB/") ### This the root of our database ### The terminal '/' is necessary
fullCov_Path <- list.files(path = SavePath, pattern = "full(c|C)ov.(RData|rds)$", recursive = TRUE)


### Variables
MCC <- c(seq(1,20,0.5))
MRG <- c(seq(10,100, 10))
MCC_MRG_Grid <- expand.grid(MCC,MRG)
if (length(fullCov_Path) == 1){
if (str_detect(fullCov_Path,".*.rds$")){  ##########?????? Do we really need fullCov for this step
  load(fullCov_Path)
} else if(stringr::str_detect(fullCov_Path, ".*.RData")) {
  fullCov <<- load(fullCov_Path)
}
} else {
  stop("There are more than one fulCov files in our working directory(cwd) or /n niehter rds nor RData files was found." )
}

OarSeqinfo <- Seqinfo(seqnames = c(1:26, "X"),
 seqlengths = c( 275612895, 248993846, 224283230, 119255633, 107901688, 117031472, 100079507,
  90695168, 94726778, 86447213, 62248096, 79100223, 83079144, 62722625, 80923592, 71719816, 
  72286588, 68604602, 60464314, 51176841, 50073674, 50832532, 62330649, 42034648, 45367442,
  44077779, 135437088),isCircular = rep(FALSE, 27), genome = "Oar3.1")

ReadLength <- c(rep(150,3), rep(100,2),rep(150, 4))

pheno <- as.data.frame( cbind(Breed = Breedtrt, Fecun = Fecuntrt,
 ThreeGroup = ThreeGroup, ReadLength = ReadLength))
pheno <- map_dfc(colnames(pheno), ~ as.factor(pheno[,.]))
colnames(pheno) <- c("Breedtrt","Fecuntrt","ThreeGroup","ReadLength")
rownames(pheno) <- NULL

Breeds <- c("Shall","San")

bamfilespath <- list.files(path = paste0(SavePath,"DATA/Analyze/Out/Sheep/")
                           , pattern = ".Aligned.out.bam", all.files = TRUE, full.names = TRUE
                           , recursive = TRUE )

readsnames <- gsub(".*/(.*)_Aligned.out.bam", "\\1", bamfilespath)
names(bamfilespath) <- readsnames

Coordbamfilespath <- list.files(path = paste0(SavePath,"DATA/Analyze/Out/Sheep/")
                                , pattern = "sortedByCoord.out.bam$", full.names = TRUE
                                , recursive = TRUE)

Coordbamfilespathindex <- paste0(Coordbamfilespath, ".bai")

Bamfiles <- map(seq_along(Coordbamfilespath), function(i) {BamFile(file = Coordbamfilespath[i],
 index = Coordbamfilespathindex[i])})

bamfileslist <-  BamFileList(Bamfiles)
names(bamfileslist) <- readsnames

pheno$Samples <- names(bamfileslist)

# Creating Metadata columns and TxDb Object
met <- data.frame( IllSeqSheep = c(paste("LF_San_", rep(1:3,1), sep = "" )
                                   , c("HF_Shall_1", "LF_Shall_1","LF_Shall_2",
                                    "HF_Shall_2","HF_Shall_3", "LF_Shall_3")))
metadata(bamfileslist) <- met

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

genomicState <- makeGenomicState(txdb =  TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1::TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1
                                 , chrs = names(seqlengths(bamfileslist)[1:27])
                                 , style = getOption("chrsStyle", "Ensembl"), currentStyle = "Ensembl")
seqinfo(genomicState) <- OarSeqinfo
save(object = genomicState,  file = paste0(SavePath,"genomicState.RData"))
