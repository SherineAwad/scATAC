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
nATAC1 = as.double(args [2]) 
nRNA1 = as.double(args [3]) 
nATAC2 = as.double(args [4])
nRNA2 = as.double(args[5])
features = as.double(args[6])
nucl = as.double(args [7]) 
tss = as.double(args [8]) 
mt = as.double(args [9])

nATAC1
nRNA1
nATAC2
nRNA2
features 
nucl 
tss 
mt

myRDS = paste(mysample, "_preprocessed.rds", sep="") 
myRDS
myObject <- readRDS(myRDS)


#myObject <- subset(x = myObject,subset = nCount_ATAC < 100000 & nCount_RNA < 30000 & nCount_ATAC > 1000 & nCount_RNA > 1000 & nFeature_RNA > 500 & nucleosome_signal < 1.2 & TSS.enrichment > 2 & percent.mt < 15)
myObject <- subset(x = myObject,subset = nCount_ATAC < nATAC1 & nCount_RNA < nRNA1 & nCount_ATAC > nATAC2 & nCount_RNA > nRNA2 & nFeature_RNA > features & nucleosome_signal < nucl & TSS.enrichment > tss & percent.mt < mt)


#myObject <- subset(x = myObject,subset = nCount_ATAC < 100000 & nCount_RNA < 30000 & nCount_ATAC > 1000 & nCount_RNA > 1000 & nFeature_RNA > 500 & nucleosome_signal < 1.0  & TSS.enrichment > 5 & percent.mt < 20) 
myRDS <- paste(mysample, "_filtered.rds", sep="") 
saveRDS(myObject, file = myRDS)



