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
myCell <- args[2] 
myRDS <- paste(mysample, "_diffPeaks.rds", sep="")
myRDS

myObject <- readRDS(myRDS)

DefaultAssay(myObject) <- "ATAC"

myfile = paste(myCell, "_sorted.csv", sep="")
myfile 
df <- readfile <- read.csv(myfile, head = TRUE, sep=",")

regions = df[1] 
top_regions = regions[1:15,]


figure_name <- ""
figure_name <- paste(myCell, "_coverage.pdf", sep="")
pdf(file =figure_name, width=20, height= 60)
p1 <- CoveragePlot(
  object = myObject,
  region = top_regions[1],
  annotation = TRUE,
  peaks = TRUE,
  extend.upstream = 2000,
  extend.downstream = 2000,
)

p2 <- CoveragePlot(
  object = myObject,
  region = top_regions[2],
  annotation = TRUE,
  peaks = TRUE,
  extend.upstream = 2000,
  extend.downstream = 2000,
)

p3 <- CoveragePlot(
  object = myObject,
  region = top_regions[3],
  annotation = TRUE,
  peaks = TRUE,
  extend.upstream = 2000,
  extend.downstream = 2000,
)

p4 <- CoveragePlot(
  object = myObject,
  region = top_regions[4],
  annotation = TRUE,
  peaks = TRUE,
  extend.upstream = 2000,
  extend.downstream = 2000,
)

p5 <- CoveragePlot(
  object = myObject,
  region = top_regions[5],
  annotation = TRUE,
  peaks = TRUE,
  extend.upstream = 2000,
  extend.downstream = 2000,
)

p6 <- CoveragePlot(
  object = myObject,
  region = top_regions[6],
  annotation = TRUE,
  peaks = TRUE,
  extend.upstream = 2000,
  extend.downstream = 2000,
)

p7 <- CoveragePlot(
  object = myObject,
  region = top_regions[7],
  annotation = TRUE,
  peaks = TRUE,
  extend.upstream = 2000,
  extend.downstream = 2000,
)

p8 <- CoveragePlot(
  object = myObject,
  region = top_regions[8],
  annotation = TRUE,
  peaks = TRUE,
  extend.upstream = 2000,
  extend.downstream = 2000,
)
p9 <- CoveragePlot(
  object = myObject,
  region = top_regions[9],
  annotation = TRUE,
  peaks = TRUE,
  extend.upstream = 2000,
  extend.downstream = 2000,
)
p10 <- CoveragePlot(
  object = myObject,
  region = top_regions[10],
  annotation = TRUE,
  peaks = TRUE,
  extend.upstream = 2000,
  extend.downstream = 2000,
)

p11 <- CoveragePlot(
  object = myObject,
  region = top_regions[11],
  annotation = TRUE,
  peaks = TRUE,
  extend.upstream = 2000,
  extend.downstream = 2000,
)

p12 <- CoveragePlot(
  object = myObject,
  region = top_regions[12],
  annotation = TRUE,
  peaks = TRUE,
  extend.upstream = 2000,
  extend.downstream = 2000,
)

p13 <- CoveragePlot(
  object = myObject,
  region = top_regions[13],
  annotation = TRUE,
  peaks = TRUE,
  extend.upstream = 2000,
  extend.downstream = 2000,
)

p14 <- CoveragePlot(
  object = myObject,
  region = top_regions[14],
  annotation = TRUE,
  peaks = TRUE,
  extend.upstream = 2000,
  extend.downstream = 2000,
)

p15 <- CoveragePlot(
  object = myObject,
  region = top_regions[15],
  annotation = TRUE,
  peaks = TRUE,
  extend.upstream = 2000,
  extend.downstream = 2000,
)

patchwork::wrap_plots(p1, p2, p3,p4,p5,p6,p7,p8,p9,p10,p11, p12,p13,p14,p15, ncol=1)
dev.off()

