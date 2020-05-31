library(purrr)
library(ggplot2)

Y <- rnorm(n = 1000000, mean = 98,sd = 25)
Samples <- map(1:100,~ sample(Y,size = 100*(.^2),replace = FALSE))
Means <- data.frame(Mean = unlist(map(Samples, ~mean(.))))
SDs <- data.frame(SD = unlist(map(Samples, ~sd(.))))
SEs <- data.frame(SE = unlist(map(Samples, ~ sd(.)/length(Samples))))

ggplot(Means, aes(x = 1:nrow(Means), y = Means[,1])) + geom_smooth(se = FALSE) + ylab("Means") + xlab("Sample Size") +
geom_hline(aes(yintercept = mean(Y)))

ggplot(SEs, aes(x = 1:nrow(SEs), y = SEs[,1])) + geom_smooth(se = FALSE) + ylab("SE") + xlab("Sample Size")

ggplot(SDs, aes(x = 1:nrow(SDs), y = SDs[,1])) + geom_smooth(se = FALSE) + ylab("SD") + xlab("Sample Size") +
geom_hline(aes(yintercept = sd(Y)))

sd(Y)/length(Y)

########################################################################################### Code repository ########################################################################################################
Testpath <- file.path(path = "/home/MPourjam/RNA-Seq/RProject/RNA-Seq_Oar3.1_V/Results/TestFilteredCov2/")

if (any(grepl("^fullRegions.Rdata$", dir(path = Testpath )))){
  load(paste0(SavedData, "results.RData"))
} else {
  system.time(results <- map(seq_along(filteredCov), ~ analyzeChr(chr = names(filteredCov[.]), coverageInfo = filteredCov[[.]], models[[3]],
                                                                  groupInfo = as.factor(ThreeGroup), writeOutput = TRUE, cutoffFstat = 5e-02, cutoffType = "theoretical",
                                                                  nPermute = 50, seeds = 20190312 + seq_len(50), returnOutput = TRUE, txdb = TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1,
                                                                  runAnnotation = TRUE, genomicState = genomicState$fullGenome, annotationPackage = NULL, mc.cores = 7)))
  names(results) <- names(regionMat)
  analyzeChrWarnings <- warnings()
  save(results, file = paste0(SavedData, ls()[which(str_detect(ls(), "^results$"))], ".RData", sep = ""))
}

###################################################

# Working with FilteredCov = 2 
setwd("./TestFilteredCov2/")
mergeResults(chrs=c(1:26,"X"), prefix=".",
             genomicState = genomicState$fullGenome, 
             optionsStats = optionsStats)

load("./fullRegions.Rdata")
load("./fullNullSummary.Rdata")
load("./fullAnnotatedRegions.Rdata")
load("./Results_Fecuntrt_NULL_mean_11.RData")
setwd("~/RNA-Seq/RProject/RNA-Seq_Oar3.1_V/Results")

report <- derfinderReport(prefix='TestFilteredCov2', browse=FALSE,
                          nBestRegions=15, makeBestClusters=TRUE, outdir='html',
                          fullCov= fullCov, optionsStats = optionsStats, hg19 = FALSE,txdb = TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1, p.ideos = FALSE)
?plotOverview()

for (w in list.dirs("~/RNA-Seq/RProject/RNA-Seq_Oar3.1_V/Results", recursive = FALSE)) {
  setwd(w)
  load("./fullRegions.Rdata")
  load("./fullNullSummary.Rdata")
  load("./fullAnnotatedRegions.Rdata")
  load("./optionsMerge.Rdata")
  derfinderReport(prefix='.', browse=FALSE,
                  nBestRegions=25, makeBestClusters=TRUE, outdir='analysisResults',
                  fullCov= fullCov, optionsStats = optionsMerge$optionsStats, hg19 = FALSE,txdb = TxDb.Oaries.BioMart.ENSEMBLMARTENSEMBL.Oarv3.1, p.ideos = FALSE)
  
}
rm(w)



