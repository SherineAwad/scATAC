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

myRDS <- paste(mysample, ".rds", sep="")
myRDS

myObject <- readRDS(myRDS)

#### Analyze RNA part

DefaultAssay(myObject) <- "RNA"

myObject <- SCTransform(myObject, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'UMAP_') #, useNames=FALSE)
#### Analyze ATAC part

DefaultAssay(myObject) <- "ATAC"
myObject <- RunTFIDF(myObject)
myObject <- FindTopFeatures(myObject, min.cutoff = 20)
myObject <- RunSVD(myObject)
myObject <- RunUMAP(myObject, dims = 2:30, reduction = 'lsi',reduction.name = "umap.atac", reduction.key = "atacUMAP_")

figure_name <- "" 
figure_name <- paste(mysample, "depth_cor.pdf", sep="")
pdf(file=figure_name, width=12)
DepthCor(myObject)
dev.off()

### Combine RNA and ATAC with WNN analysis

myObject <- FindMultiModalNeighbors(myObject, reduction.list = list("pca", "lsi"),dims.list = list(1:30, 2:30))

myObject <- RunUMAP(myObject, nn.name = "weighted.nn", reduction.name = "umap.wnn", reduction.key ="wnnUMAP_")

myObject <- FindClusters(myObject, graph.name = "wsnn", resolution = 1.2, algorithm = 3, verbose = FALSE)


myObject$Sample <- myObject$orig.ident


figure_name <- ""
figure_name <- paste(mysample, "Clusters.pdf", sep="")
pdf(file =figure_name, width =12)
DimPlot(myObject, reduction = "umap.wnn", group.by = "orig.ident",  repel = TRUE) + ggtitle("WNN")
DimPlot(myObject, reduction = "umap.wnn", label=TRUE, repel = TRUE) + ggtitle("WNN")
dev.off()


head(myObject)
figure_name <- ""
figure_name <- paste(mysample, "dim_plot.pdf", sep="")
pdf(file=figure_name, width=12)
DimPlot(object = myObject, label = TRUE) + NoLegend() 
dev.off()


p1 <- DimPlot(myObject, reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("RNA")

p2 <- DimPlot(myObject, reduction = "umap.atac",group.by = "seurat_clusters", label = TRUE,  repel = TRUE) + ggtitle("ATAC")

p3 <- DimPlot(myObject, reduction = "umap.wnn", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + ggtitle("WNN")

figure_name <- ""
figure_name <- paste(mysample, "Cluster_VlnPlot1.pdf", sep="")
pdf(file=figure_name, width=12)
p1+p3 & NoAxes()
dev.off()

figure_name <- ""
figure_name <- paste(mysample, "Cluster_VlnPlot2.pdf", sep="")
pdf(file=figure_name, width=12)
p2+p3 & NoAxes()
dev.off()


p1 <- DimPlot(myObject, reduction = "umap.rna", group.by = "orig.ident", repel = TRUE) + ggtitle("RNA")

p3 <- DimPlot(myObject, reduction = "umap.wnn", group.by = "orig.ident",  repel = TRUE) + ggtitle("WNN")

figure_name <- ""
figure_name <- paste(mysample, "ident_VlnPlot.pdf", sep="")

pdf(file=figure_name, width=12) 
p1+p3 & NoAxes()
dev.off()

Idents(myObject)
head(Idents(myObject) ) 

DefaultAssay(myObject) <- "ATAC"
myRDS <- paste(mysample, "_analysed.rds", sep="")
saveRDS(myObject, file = myRDS)

DefaultAssay(myObject) <- "RNA"
figure_name <-"" 
figure_name <- paste(mysample, "_RNA_WNN.pdf", sep="") 
pdf(file =figure_name, width=12) 
FeaturePlot(myObject, features = c("sct_Rbfox3", "sct_Sebox", "sct_Gad1", "sct_Elavl3"), reduction = "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject, features = c("sct_Sox9", "sct_Glul",  "sct_Rlbp1", "sct_Ascl1"), reduction =  "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject, features = c( "sct_Otx2", "sct_Olig2", "sct_Crx", "sct_Neurog2"), reduction = "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject, features = c("sct_Rho", "sct_Arr3", "sct_Tfap2b", "sct_Vsx1"), reduction =  "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject, features = c("sct_Insm1", "sct_Prdm1", "sct_GFP", "sct_Elavl4"), reduction =  "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject, features = c("sct_Gnat1", "sct_Pcp2", "sct_Prkca","sct_Sebox"), reduction = "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject, features = c("sct_Cabp5", "sct_Isl1", "sct_Slc6a9","sct_Gad2"), reduction =  "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject, features = c("sct_Chat", "sct_Pou4f2", "sct_Rbpms", "sct_Lhx1"), reduction = "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject, features = c("sct_Csf1r", "sct_Ccr2", "sct_Pax2", "sct_Kcnj8"), reduction =  "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject, features = c("sct_Rpe65", "sct_Acta2", "sct_Tie1", "sct_Klf4"), reduction =  "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject, features = c("sct_Hbb-bt", "sct_Grm6", "sct_Grik1"), reduction =  "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
dev.off()

myRDS <- paste(mysample, "_analysed.rds", sep="")
saveRDS(myObject, file = myRDS)

