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

#head(Idents(myObject), 5)

aPeaks <- AccessiblePeaks(
  myObject,
  min.cells = 10
)
head(aPeaks) 


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
myObject <- LinkPeaks(object = myObject, peak.assay = "ATAC", expression.assay = "SCT", genes.use = c("Ascl1", "Otx2"))
print("done linkpeaks") 

figure_name <- "" 
figure_name <- paste(mysample, "linkpeaks_coverage.pdf", sep="") 
pdf(file =figure_name, width =12)
p1 <- CoveragePlot(myObject,
                   region = "Ascl1",
                   features = "Ascl1",
                   group.by = "seurat_clusters",
                   extend.upstream = 1000,
                   extend.downstream = 1000)
p2 <- CoveragePlot(myObject,
                   region = "Otx2",
                   features = "Otx2",
                   group.by = "seurat_clusters",
                   extend.upstream = 1000,
                   extend.downstream = 1000)
patchwork::wrap_plots(p1, p2, ncol = 1)
dev.off()


