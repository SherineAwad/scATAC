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

mysample1 <- args[1] 
mysample2 <- args[2] 
mysample3 <- args[3] 
mysample4 <- args[4] 
myRDS <- args[5] 

myObject1 <- paste(mysample1,".rds", sep="")
myObject2 <- paste(mysample2,".rds", sep="") 
myObject3 <- paste(mysample3,".rds", sep="")
myObject4 <- paste(mysample4,".rds", sep="")
myRDS <- paste(myRDS, ".rds", sep="")

myObject1
myObject2
myObject3
myObject4 
myRDS 

myObject1 <- readRDS(myObject1) 
myObject2 <- readRDS(myObject2) 
myObject3 <- readRDS(myObject3) 
myObject4 <- readRDS(myObject4) 

myObject <- merge(myObject1, y = c(myObject2, myObject3, myObject4), add.cell.ids = c("mysample1", "mysample2", "mysample3", "mysample4"), project = "mergedObject")

saveRDS(myObject, file = myRDS)






