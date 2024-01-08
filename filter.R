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

mysample = args[1]
myRDS <- paste(mysample, ".rds", sep="") 

myObject <- readRDS(myRDS)


high_ns <- quantile(myObject[["nucleosome_signal"]]$nucleosome_signal, probs = 0.95) 
 
low_ts <- quantile(myObject[["TSS.enrichment"]]$TSS.enrishment, probs =0.09) 

high_ns 
low_ts

myObject <- subset(x = myObject,subset = nCount_ATAC < 100000 & nCount_RNA < 30000 & nCount_ATAC > 1000 & nCount_RNA > 1000 & nFeature_RNA > 500 & nucleosome_signal < high_ns & TSS.enrichment > low_ts  & percent.mt < 15)

myRDS <- paste(mysample, "_filtered.rds", sep="") 
saveRDS(myObject, file = myRDS)




