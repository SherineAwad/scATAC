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
myRDS <- paste(mysample, "_diffPeaks.rds", sep="")

mysample
myRDS

myObject <- readRDS(myRDS)

DefaultAssay(myObject) <- "ATAC"

head(myObject)

#head(colnames(myObject))
#head(Cells(myObject))


Idents(myObject)
table(Idents(myObject))
prop.table(table(Idents(myObject)))

