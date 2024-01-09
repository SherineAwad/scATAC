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
mytsv <- paste(mysample, "_atac_fragments.tsv.gz", sep="")

myRDS <- paste(mysample, ".rds", sep="")
myRDS

myObject <- readRDS(myRDS)

head(myObject)

#### Analyze RNA part

DefaultAssay(myObject) <- "RNA"

myObject <- SCTransform(myObject, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'UMAP_')

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

myObject <- FindClusters(myObject, graph.name = "wsnn", resolution = 0.8, algorithm = 3, verbose = FALSE)

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

DefaultAssay(myObject) <- "RNA"
figure_name <-"" 
figure_name <- paste(mysample, "Genes_merge_RNA_WNN.pdf", sep="") 
pdf(file =figure_name, width=12) 
FeaturePlot(myObject, features = c("sct_Rbfox3", "sct_Sebox", "sct_Gad1", "sct_Elavl3"), reduction = "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject, features = c("sct_Sox9", "sct_Glul",  "sct_Rlbp1", "sct_Ascl1"), reduction =  "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject, features = c( "sct_Otx2", "sct_Olig2", "sct_Crx", "sct_Neurog2"), reduction = "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject, features = c("sct_Rho", "sct_Arr3", "sct_Tfap2b", "sct_Vsx1"), reduction =  "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
FeaturePlot(myObject, features = c("sct_Insm1", "sct_Prdm1", "sct_GFP", "sct_Elavl4"), reduction =  "umap.wnn", cols = c("lightgrey", "red"), pt.size = 0.01)
dev.off()

myObject$sample <- myObject$orig.ident 

# Finder all markers
DefaultAssay(myObject) <- "ATAC"
myObject.atac.markers <- FindAllMarkers(myObject, assay = "ATAC", test.use = "LR", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
myObject.atac.markers %>%
  group_by(cluster) %>%
  top_n(n = 2, wt = avg_log2FC)

file_name <- paste(mysample, "MG_DEGs_CtrlvsKO.csv", sep="")
write.csv(myObject.atac.markers, file=file_name)

DefaultAssay(myObject) <- "ATAC"

saveRDS(myObject, file = myRDS)


