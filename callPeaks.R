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
set.seed(1234)

args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]

myRDS <- paste(mysample, "_FinalClusters.rds", sep="") 
myRDS

myObject <- readRDS(myRDS)

library(chromVAR)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)   ## For plotting graph, similar  ggplot2
library(qlcMatrix)   ### For linking gene
library('ggforce')


#####
DefaultAssay(myObject) <- "ATAC"
peaks <- CallPeaks(object = myObject,  outdir ="/nfs/turbo/umms-thahoang/sherine/scATAC", cleanup = FALSE)
# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(myObject),
  features = peaks,
  cells = colnames(myObject)
)


file_name <- paste(mysample, "peaks.csv", sep="")
write.csv(macs2_counts, file=file_name)
write.csv(peaks,"peaks.txt", sep="") 

print("Macs2 counts done")

rtracklayer::export.bed(peaks, "my_peaks.bed")

myRDS <- paste(mysample, "_testpeaks.rds", sep="")
saveRDS(myObject, file = myRDS)
