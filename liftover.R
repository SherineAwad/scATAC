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
library(ChIPpeakAnno)
library(rtracklayer)

grObject <- GRanges(seqnames=c("chr1"), ranges=IRanges(start=10037618, end=10038530))

# import the chain file
chainObject <- import.chain("mm10ToHg38.over.chain")

# run liftOver
results <- as.data.frame(liftOver(grObject, chainObject))
head(results)

