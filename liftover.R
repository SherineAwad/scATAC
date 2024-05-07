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


args <- commandArgs(trailingOnly = TRUE)
myCell <- args[1]

myfile = paste(myCell, "_sorted.csv", sep="")
myfile
outfile = paste(myCell, "mtohliftover.csv", sep="_") 
outfile 

df <- readfile <- read.csv(myfile, head = TRUE, sep=",")

regions = df[1]
top_regions = regions[1:20,]
top_regions 

results  <- data.frame(Date=as.Date(character()),
                 File=character(), 
                 User=character(), 
                 stringsAsFactors=FALSE) 
for (i in 1:20)
{
s <- unlist(strsplit(top_regions[i], split = "-") )
chr <- s[1]
st <- strtoi(s[2]) 
en <- strtoi(s[3])

grObject <- GRanges(seqnames=c(chr), ranges=IRanges(start=st, end=en))

# import the chain file
chainObject <- import.chain("mm10ToHg38.over.chain")

# run liftOver

newR<- as.data.frame(liftOver(grObject, chainObject))
results <- rbind(results, newR)
}
#write.table(results, outfile, sep = ",", append = T,col.names=F)write.csv(results, outfil=, row.names = FALSE) 
write.csv(results, file=outfile, row.names = FALSE) 

