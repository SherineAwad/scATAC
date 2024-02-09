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
myRDS <- paste(mysample, "_diffPeaks.rds", sep="")
myRDS

myObject <- readRDS(myRDS)


if (FALSE) { 
granges <- keepStandardChromosomes(granges(myObject), pruning.mode='coarse')
print("granges")
head (granges)
main.chroms <- standardChromosomes(BSgenome.Mmusculus.UCSC.mm10)
print("main.chroms")
head(main.chroms)
keep.peaks <- which(as.character(seqnames(granges(myObject))) %in% main.chroms)
print("keep.peaks")
head(keep.peaks)
myObject[["ATAC"]] <- subset(myObject[["ATAC"]], features = rownames(myObject[["ATAC"]])[keep.peaks])

## Linking peaks to genes ###
myObject <- RegionStats(myObject, genome = BSgenome.Mmusculus.UCSC.mm10)
print("done region stats")
#myObject <- LinkPeaks(object = myObject, peak.assay = "ATAC", expression.assay = "SCT", genes.use = c("Ascl2", "Otx2"))
myObject <- LinkPeaks(object = myObject, peak.assay = "ATAC", expression.assay = "SCT")
print("done linkpeaks")

myRDS <- paste(mysample, "_linkPeaks.rds", sep="")
saveRDS(myObject, file = myRDS)

head(myObject) 

} 



