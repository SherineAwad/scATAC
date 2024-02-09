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

myRDS <- paste(mysample, "_renamedClusters.rds", sep="")
myRDS

myObject <- readRDS(myRDS)


cell_values <- c("Muller/AC doublets", "Muller/BC doublets", "Muller/Cone doublets","Muller/rod doublets")
mySubset <- subset(myObject, idents = cell_values, invert = TRUE)


figure_name <- ""
figure_name <- paste(mysample, "FinalClustersUMAP.pdf", sep="")
pdf(file =figure_name, width =12)
DimPlot(mySubset, reduction = "umap.wnn", group.by = "orig.ident",  repel = TRUE) + ggtitle("WNN")
DimPlot(mySubset, reduction = "umap.wnn", label=TRUE, repel = TRUE) + ggtitle("WNN")
dev.off()

head(Idents(mySubset))
myRDS <- paste(mysample, "_FinalClusters.rds", sep="")
saveRDS(mySubset, file = myRDS)

