library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Mmusculus.v79)
library(BSgenome.Mmusculus.UCSC.mm10)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
library(tidyr)
set.seed(1234)



df1 <- read.csv("merged3samplesDEGs.csv", head = TRUE, sep=",")
df2 <- read.csv("merged3samplesannotatedDEGs.csv", head = TRUE, sep=",") 

r <- merge(df1, df2,by=c("pos","pos"),all.x=F)

write.csv(r,"DEGs.csv" , all(T) )


