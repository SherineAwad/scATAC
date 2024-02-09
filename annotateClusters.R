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

myRDS <- paste(mysample, "_analysed.rds", sep="")
myRDS

myObject <- readRDS(myRDS)

#head(myObject) 
#levels(myObject)
#head(Idents(myObject) 


myObject <- RenameIdents(
  object = myObject,
  "0" = 'Cone', 
  "1" = 'Rod', 
  "2" = 'Rod',
  "3" = 'Muller glia', 
  "4" = 'Muller glia', 
  "5" = 'OFF cone BC',
  "6" = 'Rod', 
  "7" = 'Glycinergic AC', 
  "8" = 'Rod', 
  "9" = 'RGC', 
  "10" = 'Rod BC', 
  "11" = 'Glycinergic AC', 
  "12" = 'OFF cone BC', 
  "13" = 'ON cone BC', 
  "14" = 'GABAergic AC', 
  "15" = 'ON cone BC', 
  "16" = 'GABAergic AC',
  "17" = 'Rod BC', 
  "18" = 'ON cone BC', 
  "19" = 'Cone', 
  "20" = 'Rod BC', 
  "21" = 'Starbust AC', 
  "22" = 'Glycinergic AC', 
  "23" = 'Muller glia',
  "24" = 'OFF cone BC', 
  "25" = 'Rod BC',
  "26" = 'Glycinergic AC',
  "27" = 'ON cone BC',
  "28" = 'OFF cone BC',
  "29" = 'OFF cone BC',
  "30" = 'Glycinergic AC',
  "31" = 'HC',
  "32" = 'Glycinergic AC',
  "33" = 'Glycinergic AC',
  "34" = 'GABAergic AC',
  "35" = 'Rod BC',
  "36" = 'GABAergic AC',
  "37" = 'GABAergic AC',
  "38" = 'GABAergic AC',
  "39" = 'ON cone BC',
  "40" = 'Glycinergic AC',
  "41" = 'GABAergic AC',
  "42" = 'Glycinergic AC',
  "43" = 'GABAergic AC',
  "44" = 'GABAergic AC',
  "45" = 'Muller/rod doublets',
  "46" =  'GABAergic AC',
  "47" = 'Glycinergic AC',
  "48" = 'ON cone BC',
  "49" = 'OFF cone BC',
  "50" = 'GABAergic AC',
  "51" = 'Cone',
  "52" = 'GABAergic AC',
  "53" = 'Endothelial cells',
  "54" = 'Glycinergic AC',
  "55" = 'GABAergic AC',
  "56" = 'Glycinergic AC',
  "57" = 'GABAergic AC',
  "58" = 'GABAergic AC',
  "59" = 'GABAergic AC',
  "60" = 'Microglia',
  "61" = 'Glycinergic AC',
  "62" = 'GABAergic AC', 
  "63" = 'GABAergic AC',
  "64" = 'Glycinergic AC',
  "65" = 'ON cone BC',
  "66" = 'Muller/BC doublets',
  "67" = 'Muller/Cone doublets',
  "68" = 'HC',
  "69" = 'Starbust AC',
  "70" = 'GABAergic AC',
  "71" = 'Muller/AC doublets',
  "72" = 'Muller/BC doublets',
  "73" = 'GABAergic AC') 
  
head(Idents(myObject))
myRDS <- paste(mysample, "_renamedClusters.rds", sep="")
saveRDS(myObject, file = myRDS)

