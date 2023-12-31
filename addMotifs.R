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
library(patchwork) 
library(qlcMatrix) 
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
print(da_peaks)
head(da_peaks) 


atac_small = RegionStats(myObject, genome = BSgenome.Mmusculus.UCSC.mm10)
print(atac_small) 
head(atac_small@assays$peaks@meta.features)
class(atac_small@assays$peaks@meta.features)
de.motif <- head(rownames(atac_small))
bg.peaks <- tail(rownames(atac_small))
enriched.motifs <- FindMotifs(
  object = atac_small,
  features = de.motif,
  background = bg.peaks
)

figure_name <- "" 
figure_name <- paste(mysample, "_enriched_motifs.pdf", sep="")
pdf(file=figure_name, width =12) 
MotifPlot(
  object = myObject,
  motifs = head(rownames(enriched.motifs))
)
dev.off() 

