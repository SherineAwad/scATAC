library(valr)
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


args <- commandArgs(trailingOnly = TRUE)

mysample = args[1]
myRDS <- paste(mysample, "_diffPeaks.rds", sep="")

mysample
myRDS

myObject <- readRDS(myRDS)

DefaultAssay(myObject) <- "ATAC"

df <- readfile <- read.csv("merged3samplesDEGs.csv", head = TRUE, sep=",")
pos <- separate(df, col="gene", into=c('chr', 'start', 'end'), sep='-')
head(pos)
pos = subset(pos, select = -c(1:7) )
write.csv(pos, file="pos.csv")
mybed <- makeGRangesFromDataFrame(pos,ignore.strand=FALSE)
write.csv(mybed, file="mybed.csv")
rtracklayer::export.bed(mybed, "fdiffPeaks.bed", "bed")
#bed <- gr_to_bed(mybed)
#write.table(bed, file="diffPeaks2.bed", sep="\t", col.names = F, row.names = F)
