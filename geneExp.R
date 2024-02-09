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
myRDS <- paste(mysample, "_analysed_renamedClusters.rds", sep="")
myRDS

myObject <- readRDS(myRDS)


DefaultAssay(myObject) <- "ATAC"
gene.activities <- GeneActivity(myObject)

# add the gene activity matrix to the Seurat object as a new assay and normalize it
myObject[['RNA']] <- CreateAssayObject(counts = gene.activities)
myObject <- NormalizeData(
  object = myObject,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(myObject$nCount_RNA)
)

DefaultAssay(myObject) <- 'RNA'
figure_name <- ""
figure_name <- paste(mysample, "gene_activities.pdf", sep="")
pdf(file=figure_name, width=12) 
FeaturePlot(
  object = myObject,
  features = c('sct_Rbfox3', 'sct_Sebox', 'sct_Gad1', 'sct_Elavl3'), 
  pt.size = 0.1,
  max.cutoff = 'q95',
  ncol = 3
)
dev.off() 

myRDS <- paste(mysample, "_geneExp.rds", sep="")
saveRDS(myObject, file = myRDS)



