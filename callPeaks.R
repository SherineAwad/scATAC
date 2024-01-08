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
myObject <- readRDS(myRDS)

library(chromVAR)
library(motifmatchr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(patchwork)   ## For plotting graph, similar  ggplot2
library(qlcMatrix)   ### For linking gene
library('ggforce')



#####
DefaultAssay(myObject) <- "ATAC"

if (FALSE){
peaks <- CallPeaks(object = myObject)
print("Peaks")
head(peaks)

# quantify counts in each peak
macs2_counts <- FeatureMatrix(
  fragments = Fragments(myObject),
  features = peaks,
  cells = colnames(myObject)
)
print("Macs2 counts done")
# create a new assay using the MACS2 peak set and add it to the Seurat object
myObject[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  fragments = mytsv,
  annotation = annotation
)
print("peaks chromassay")
}


figure_name <- ""
figure_name <- paste(mysample, "Link_peak.genes.pdf", sep="")
pdf(file=figure_name, width=12)
p1 <- CoveragePlot(myObject,
                   region = "Ascl1",
                   features = "Ascl1",
                   #group.by = "seurat_clusters",
                   extend.upstream = 1000,
                   extend.downstream = 1000)
p2 <- CoveragePlot(myObject,
                   region = "Otx2",
                   features = "Otx2",
                   #group.by = "seurat_clusters",
                   extend.upstream = 1000,
                   extend.downstream = 1000)
patchwork::wrap_plots(p1, p2, ncol = 1)
dev.off()





