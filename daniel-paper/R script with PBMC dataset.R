source("https://bioconductor.org/biocLite.R")
biocLite("FlowSOM")


library(FlowSOM)
setwd("/Users/danielliu/Documents/Irish Lab/Automated Population Identification/FlowSOM/Test")

fileName <-  dir(pattern="*.fcs")

fSOM <- FlowSOM(fileName,
                # Input options:
                pattern = ".fcs", compensate = FALSE, transform = TRUE, toTransform=c(13:53),
                scale = TRUE,
                # SOM options:
                colsToUse = c(13:17,19:22,24:26,28,30:32,34:40,42,43,46:48,52,53), xdim = 7, ydim = 7,
                # Metaclustering options:
                nClus = 30)
                #maxMeta = 30)


#Create a FlowSOM plot
PlotStars(fSOM[[1]])

#Get a table of clustering IDs
flowSOM.clustering <- as.matrix(fSOM[[2]][fSOM[[1]]$map$mapping[,1]])

#Stuff from the example shell (couldn't get their example version to work)
#PlotStars(fSOM$FlowSOM,
#backgroundValues = as.factor(fSOM$metaclustering))
