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

mysample
myRDS

myObject <- readRDS(myRDS)

myObject$Sample <- myObject$orig.ident

myObject <- FindClusters(myObject, graph.name = "wsnn", resolution = 1.2, algorithm = 3, verbose = FALSE)

figure_name <- "" 
figure_name <- paste(mysample, "Clusters.pdf", sep="") 
pdf(file =figure_name, width =12) 
DimPlot(myObject, reduction = "umap.wnn", group.by = "orig.ident",  repel = TRUE) + ggtitle("WNN")
DimPlot(myObject, reduction = "umap.wnn", label=TRUE, repel = TRUE) + ggtitle("WNN")
dev.off()

figure_name <- "" 
figure_name <- paste(mysample, "CLustersGenes.pdf", sep="") 
pdf(file =figure_name, width =12)
FeaturePlot(myObject , features =  c( "sct_Prdm1", "sct_Meis2", "sct_Lhx4", "sct_Ankrd33b"), reduction = "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject , features =  c( "sct_Dscam", "sct_Otor", "sct_Gnat2", "sct_Gngt2"), reduction = "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject , features =  c( "sct_Alk", "sct_Rgr", "sct_Slit2", "sct_Lrat"), reduction = "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject , features = c("sct_Rbfox3", "sct_Sebox", "sct_Gad1", "sct_Elavl3"), reduction = "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject , features = c("sct_Sox9", "sct_Glul",  "sct_Rlbp1", "sct_Ascl1"), reduction =  "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject , features = c( "sct_Otx2", "sct_Olig2", "sct_Crx", "sct_Neurog2"), reduction = "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject , features = c("sct_Rho", "sct_Arr3", "sct_Tfap2b", "sct_Vsx1"), reduction =  "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject , features = c("sct_Insm1", "sct_Prdm1", "sct_GFP", "sct_Elavl4"), reduction =  "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
dev.off()


myRDS <- paste(mysample, "_cleaned.rds", sep="")
saveRDS(myObject, file = myRDS)





