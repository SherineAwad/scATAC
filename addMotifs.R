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

library(chromVAR)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)   ## For plotting graph, similar  ggplot2
library(qlcMatrix)   ### For linking gene
library('ggforce')


# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
myObject <- AddMotifs(
  object = myObject,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

da_peaks <- FindMarkers(
  object = myObject,
  ident.1 ="7",
  only.pos = TRUE,
  test.use = 'LR',
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)
head(da_peaks) 
# get top differentially accessible peaks
top.da.peak <- rownames(da_peaks[da_peaks$p_val < 0.005, ])

summary(top.da.peak)


# test enrichment
enriched.motifs <- FindMotifs(
  object = myObject,
  features = top.da.peak
)

enriched_motifs 

figure_name <- "" 
figure_name <- paste(mysample, "_enriched_motifs.pdf", sep="")
pdf(file=figure_name, width =12) 
MotifPlot(
  object = myObject,
  motifs = head(rownames(enriched.motifs))
)
dev.off() 

