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
set.seed(1234)

args <- commandArgs(trailingOnly = TRUE)
mysample <- args[1]
myRDS <- paste(mysample, "_diffPeaks.rds", sep="")
myRDS

myObject <- readRDS(myRDS)
head(myObject)

DefaultAssay(myObject) <- "ATAC"

figure_name <- ""
figure_name <- paste(mysample, "_linkCoverage1.pdf", sep="")
pdf(file =figure_name, width=20, height= 60)
p1 <- CoveragePlot(
  object = myObject,
  region = "chr18-74025498-74026454",
  annotation = TRUE,
  peaks = TRUE,
)

p2 <- CoveragePlot(
  object = myObject,
  region = "chr12-112743042-112743995",
  annotation = TRUE,
  peaks = TRUE,
)

p3 <- CoveragePlot(
  object = myObject,
  region = "chr18-80986389-80987280",
  annotation = TRUE,
  peaks = TRUE,
)

p4 <- CoveragePlot(
  object = myObject,
  region = "chr3-108092272-108093221",
  annotation = TRUE,
  peaks = TRUE,
)
p5 <- CoveragePlot(
  object = myObject,
  region = "chr19-56721560-56722478",
  annotation = TRUE,
  peaks = TRUE,
)

patchwork::wrap_plots(p1, p2, p3, p4, p5, ncol=1) 
dev.off() 

figure_name <- ""
figure_name <- paste(mysample, "_linkCoverage2.pdf", sep="")
pdf(file =figure_name, width=20, height= 60)
p6 <- CoveragePlot(object = myObject,
  region = "Rbfox3",
  annotation = TRUE,
  peaks = TRUE,
)

p7 <- CoveragePlot(object = myObject,
  region = "Sox9",
  annotation = TRUE,
  peaks = TRUE,
)

p8 <- CoveragePlot(object = myObject,
  region = "Otx2",
  annotation = TRUE,
  peaks = TRUE,
)

p9 <- CoveragePlot(object = myObject,
  region = "Rho",
  annotation = TRUE,
  peaks = TRUE,
)

p10 <- CoveragePlot(object = myObject,
  region = "Insm1",
  annotation = TRUE,
  peaks = TRUE,
)
patchwork::wrap_plots(p6, p7, p8, p9, p10, ncol=1)
dev.off() 

figure_name <- ""
figure_name <- paste(mysample, "_linkCoverage3.pdf", sep="")
pdf(file =figure_name, width=20, height= 60)
p11 <- CoveragePlot(object = myObject,
  region = "Ascl1",
  annotation = TRUE,
  peaks = TRUE,
)

p12 <- CoveragePlot(object = myObject,
  region = "Sebox",
  annotation = TRUE,
  peaks = TRUE,
)

p13 <- CoveragePlot(object = myObject,
  region = "Glul",
  annotation = TRUE,
  peaks = TRUE,
)

p14 <- CoveragePlot(object = myObject,
  region = "Olig2",
  annotation = TRUE,
  peaks = TRUE,
)

p15 <- CoveragePlot(object = myObject,
  region = "Arr3",
  annotation = TRUE,
  peaks = TRUE,
)


patchwork::wrap_plots(p11, p12, p13, p14, p15, ncol =1) 
dev.off()

figure_name <- ""
figure_name <- paste(mysample, "_linkCoverage4.pdf", sep="")
pdf(file =figure_name, width=20, height= 60)
p16 <- CoveragePlot(object = myObject,
  region = "Gad1",
  annotation = TRUE,
  peaks = TRUE,
)

p17 <- CoveragePlot(object = myObject,
  region = "Rlbp1",
  annotation = TRUE,
  peaks = TRUE,
)

p18 <- CoveragePlot(object = myObject,
  region = "Crx",
  annotation = TRUE,
  peaks = TRUE,
)

p19 <- CoveragePlot(object = myObject,
  region = "Elavl3",
  annotation = TRUE,
  peaks = TRUE,
)

p20 <- CoveragePlot(object = myObject,
  region = "Neurog2",
  annotation = TRUE,
  peaks = TRUE,
)

patchwork::wrap_plots(p16, p17, p18, p19, p20, ncol=1)
dev.off()
