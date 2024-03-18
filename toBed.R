library(valr)
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
library(tidyr)
set.seed(1234)

df <- readfile <- read.csv("merged3samplesDEGs.csv", head = TRUE, sep=",")
pos <- separate(df, col="gene", into=c('chr', 'start', 'end'), sep='-')
cone <- subset(pos, cluster == 'Cone')
Rod <-  subset(pos, cluster == 'Rod') 
Mullerglia <-  subset(pos, cluster == 'Muller glia')
RodBC <-  subset(pos, cluster == 'Rod BC')
ONconeBC <-  subset(pos, cluster == 'ON cone BC')
OFFconeBC <-  subset(pos, cluster == 'OFF cone BC')
GABAergic <-  subset(pos, cluster == 'GABAergic AC')
Glycinergic <-  subset(pos, cluster == 'Glycinergic AC')
Starbust <-  subset(pos, cluster == 'Starbust AC')
RGC <-  subset(pos, cluster == 'RGC')
HC <-  subset(pos, cluster == 'HC')
Microglia <-  subset(pos, cluster == 'Microglia')
Endo <-  subset(pos, cluster == 'Endothelial cells')


pos_cone = subset(cone, select = -c(1:7) )
pos_Rod <-  subset(Rod, select = -c(1:7) )  
pos_Mullerglia <-  subset(Mullerglia, select = -c(1:7) )
pos_RodBC <-  subset(RodBC, select = -c(1:7) )
pos_ONconeBC <-  subset(ONconeBC, select = -c(1:7) )
pos_OFFconeBC <-  subset(OFFconeBC, select = -c(1:7) )
pos_GABAergic <-  subset(GABAergic, select = -c(1:7) )
pos_Glycinergic <-  subset(Glycinergic, select = -c(1:7) )
pos_Starbust <-  subset(Starbust, select = -c(1:7) )
pos_RGC <-  subset(RGC, select = -c(1:7) )
pos_HC <-  subset(HC, select = -c(1:7) )
pos_Microglia <- subset(Microglia, select = -c(1:7) )
pos_Endo <- subset(Endo, select = -c(1:7) )


head(pos_cone)
head(pos_Rod)
head(pos_Mullerglia)
head(pos_RodBC)
head(pos_ONconeBC)
head(pos_OFFconeBC)
head(pos_GABAergic)
head(pos_Glycinergic)
head(pos_Starbust)
head(pos_RGC)
head(pos_HC)
head(pos_Microglia)
head(pos_Endo)

mybed <- makeGRangesFromDataFrame(pos,ignore.strand=FALSE, starts.in.df.are.0based=TRUE)
rtracklayer::export.bed(mybed, "diffPeaks.bed", "bed")

mybed <- makeGRangesFromDataFrame(pos_cone,ignore.strand=FALSE, starts.in.df.are.0based=TRUE)
rtracklayer::export.bed(mybed, "cone_diffPeaks.bed", "bed")


mybed <- makeGRangesFromDataFrame(pos_Rod,ignore.strand=FALSE, starts.in.df.are.0based=TRUE)
rtracklayer::export.bed(mybed, "Rod_diffPeaks.bed", "bed")

mybed <- makeGRangesFromDataFrame(pos_Mullerglia,ignore.strand=FALSE, starts.in.df.are.0based=TRUE)
rtracklayer::export.bed(mybed, "Mullerglia_diffPeaks.bed", "bed")

mybed <- makeGRangesFromDataFrame(pos_RodBC,ignore.strand=FALSE, starts.in.df.are.0based=TRUE)
rtracklayer::export.bed(mybed, "RodBC_diffPeaks.bed", "bed")

mybed <- makeGRangesFromDataFrame(pos_ONconeBC,ignore.strand=FALSE, starts.in.df.are.0based=TRUE)
rtracklayer::export.bed(mybed, "ONconeBC_diffPeaks.bed", "bed")

mybed <- makeGRangesFromDataFrame(pos_OFFconeBC,ignore.strand=FALSE, starts.in.df.are.0based=TRUE)
rtracklayer::export.bed(mybed, "OFFconeBC_diffPeaks.bed", "bed")

mybed <- makeGRangesFromDataFrame(pos_GABAergic,ignore.strand=FALSE, starts.in.df.are.0based=TRUE)
rtracklayer::export.bed(mybed, "GABAergic_diffPeaks.bed", "bed")

mybed <- makeGRangesFromDataFrame(pos_Glycinergic,ignore.strand=FALSE, starts.in.df.are.0based=TRUE)
rtracklayer::export.bed(mybed, "Glycinergic_diffPeaks.bed", "bed")


mybed <- makeGRangesFromDataFrame(pos_Starbust,ignore.strand=FALSE, starts.in.df.are.0based=TRUE)
rtracklayer::export.bed(mybed, "Starbust_diffPeaks.bed", "bed")


mybed <- makeGRangesFromDataFrame(pos_RGC,ignore.strand=FALSE, starts.in.df.are.0based=TRUE)
rtracklayer::export.bed(mybed, "RGC_diffPeaks.bed", "bed")


mybed <- makeGRangesFromDataFrame(pos_HC,ignore.strand=FALSE, starts.in.df.are.0based=TRUE)
rtracklayer::export.bed(mybed, "HC_diffPeaks.bed", "bed")


mybed <- makeGRangesFromDataFrame(pos_Microglia,ignore.strand=FALSE, starts.in.df.are.0based=TRUE)
rtracklayer::export.bed(mybed, "Microglia_diffPeaks.bed", "bed")

mybed <- makeGRangesFromDataFrame(pos_Endo,ignore.strand=FALSE, starts.in.df.are.0based=TRUE)
rtracklayer::export.bed(mybed, "Endo_diffPeaks.bed", "bed")
