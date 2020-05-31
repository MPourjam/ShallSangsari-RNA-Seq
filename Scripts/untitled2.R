#install.packages("ggfortify")
library(ggfortify)
#install.packages(c("FactoMineR", "factoextra"))
library("FactoMineR")
library("factoextra")
#install.packages("pca3d")
source("http://bioconductor.org/biocLite.R")
#BiocManager::install("rgl")
#install.packages("rgl")
library("rgl")


regionMat$`6`$coverageMatrix

RegionMatrice <- map(names(regionMat), ~ as_tibble(regionMat[[.]][["coverageMatrix"]]))
names(RegionMatrice) <- names(regionMat)
RegionMatrice <- bind_rows(RegionMatrice, .id = "Chr")[-1]
as.matrix(RegionMatrice)
RegionMatriceT<- t(RegionMatrice)



regionPCA<- prcomp(RegionMatriceT, scale. = TRUE)
regionPCAT <- princomp(RegionMatrice, cor = TRUE)

#SamplePCA <- princomp(RegionMatrice, scores = FALSE, cor = FALSE)

autoplot(regionPCA, label = TRUE)
biplot(regionPCAT)
plot3d(regionPCAT$scores)
################################################################################################################
names(counts) <- names(regionMat)
COU <- map_dfr(counts, ~as.data.frame(.))
COUTmeanNOT0 <- t(as.matrix(COU))[,colMeans(t(as.matrix(COU))) != 0 & apply(t(as.matrix(COU)), 2, sd) != 0]
COUTmeanNOT0<- scale(COUTmeanNOT0)
 countPCA <- prcomp(COUTmeanNOT0, scale. = TRUE)
 countPCAT <- princomp(COUTmeanNOT0, scores = TRUE)
 
 autoplot(countPCA, label = TRUE)
 biplot(countPCAT)
 
################################################################################################################

str(counts)
names(counts) <- names(regionMat)
countsCombination <- map_dfr(counts, ~as_tibble(.))
countsCombinationT <- t(as.matrix(countsCombination))


countsPCA <- prcomp(countsCombination, scale. = TRUE)

countsPCAT <- princomp(countsCombinationT)

