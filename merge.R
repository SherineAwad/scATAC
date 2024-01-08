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


mysample1 <- paste(args[1],".rds", sep="")
mysample2 <- paste(args[2],".rds", sep="") 
mysample3 <- paste(args[3],".rds", sep="")
mysample4 <- paste(args[4],".rds", sep="")

myRDS <- paste(args[5], ".rds", sep="")

myObject1 <- readRDS(mysample1) 
myObject2 <- readRDS(mysample2) 
myObject3 <- readRDS(mysample3) 
myObject4 <- readRDS(mysample4) 


myObject <- merge(myObject1, y = c(myObject2, myObject3, myObject4), add.cell.ids = c("myObject1", "myObject2", "myObject3", "myObject4"), project = "myObject")

saveRDS(myObject, file = myRDS)






