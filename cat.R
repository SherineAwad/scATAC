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


#rename column query region in merged3samplesannotatedDEGs.csv to gene to be the same as the gene column in merged3samplesDEGs.csv

df1 <- read.csv("merged3samplesDEGs.csv", head = TRUE, sep=",")
df2 <- read.csv("merged3samplesannotatedDEGs.csv", head = TRUE, sep=",") 

dim(df1)
dim(df2)

unique_df1 <- unique(df1[c("gene", "cluster","p_val","avg_log2FC","p_val_adj")]) 
head(unique_df1) 

unique_df2 <- unique(df2[c("query_region","tx_id","gene_name","gene_id","gene_biotype","type","closest_region")]) 
head(unique_df2) 

r <- merge(unique_df1, unique_df2,  by.x = "gene", by.y ="query_region",all=F)
dim(r)
write.csv(r,"allDEGs.csv", all(F), row.names= FALSE)


