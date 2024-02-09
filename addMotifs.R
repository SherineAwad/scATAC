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
head(myObject)

dapeaks <- FindAllMarkers(myObject,only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

top.da.peak <- rownames(dapeaks[dapeaks$p_val < 0.005, ])
write.csv(top.da.peak, file ="topMotifs.csv")



enriched.motifs <- FindMotifs(
  object = myObject,
  features = top.da.peak[1:200]
)

write.csv(enriched.motifs, file ="enrichedMotifs.csv")


DefaultAssay(myObject) <- "ATAC"

#Motifs =c("MA1513.1", "MA1653.1","MA0467.1","MA1564.1","MA0685.1","MA0685.1","MA0599.1")  
#figure_name <- ""
#figure_name <- paste(mysample, "enrichedMotifs.pdf",sep="")
#pdf(file = figure_name, width=20)
#FeaturePlot(myObject, features = Motifs, cols = c("lightgrey", "red"),order=TRUE, min.cutoff = 'q10',max.cutoff = 'q90',pt.size = 0.01, reduction = 'umap.wnn') & NoAxes() 
#dev.off()

figure_name <- ""
figure_name <- paste(mysample, "enrichedMotifs.pdf",sep="")
pdf(file = figure_name, width=20)
MotifPlot(
  object = myObject,
  motifs = head(rownames(enriched.motifs))
)
dev.off() 


DefaultAssay(myObject) <- "ATAC"
myRDS <- paste(mysample, "_Motifs.rds", sep="")
saveRDS(myObject, file = myRDS)
