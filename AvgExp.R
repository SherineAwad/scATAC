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
library(ChIPpeakAnno)
set.seed(1234)

args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]
myRDS <- paste(mysample, "_diffPeaks.rds", sep="")
myRDS

myObject <- readRDS(myRDS)

DefaultAssay(myObject) <- "RNA"

myObject <- NormalizeData(object = myObject)
Avg <- AverageExpression(object = myObject, assays = "RNA")
write.csv(Avg,file ="AvgExp.csv")

