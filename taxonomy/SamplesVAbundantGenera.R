# Samples vs abundant Genera for plot to send to Viral team
# Read in data
setwd("~/Documents/Postgrad/CF Hackathon/Day3minIONReports/PavianProcessedMiniIONReportsDay3") #change to your directory
GenusTable <- na.omit(as.matrix(read.csv("genus_bacteria_corrected.csv", check.names = FALSE, row.names = 1)))

library(ape)
#RelAbund
RelabundGenera <- sweep(GenusTable,2,colSums(GenusTable),"/") #Sweeps cols, divides by col sum
#Filter
#Selects only the rows where the maximum value is greater than 10%
FilteredGeneraRel <- RelabundGenera[which(apply(RelabundGenera,1,max)>0.1),]
FGR <- as.data.frame(t(FilteredGeneraRel))
write.table(FGR, file = "FilteredGenera.csv", col.names=NA)

#Upload full non-filt
GR <- as.data.frame(RelabundGenera)
write.table(GR, file = "RelGenera_NoFilter.csv", col.names=NA)
