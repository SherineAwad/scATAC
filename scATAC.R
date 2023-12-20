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

mysample = args[1]
myh5 <- paste(mysample,"_filtered_feature_bc_matrix.h5", sep="")
mytsv <- paste(mysample, "_atac_fragments.tsv.gz", sep="")
myRDS <- paste(mysample, ".rds", sep="") 

mysample 
myh5
mytsv 
myRDS 
## Prepare genome annotation 
annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotation) <- 'UCSC'
genome(annotation) <- "mm10"
# Check annotation seqnames
unique(seqnames(annotation))


###load sample, the RNA and ATAC data
counts <- Read10X_h5(myh5)
head(counts)
# Check counts
lapply(counts, dim)
lapply(counts, class)

fragpath <-(mytsv)

# create a Seurat object containing the RNA adata
myObject <- CreateSeuratObject(counts = counts$`Gene Expression`,assay = "RNA", project=mysample)

# create ATAC assay and add it to the object
myObject[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,sep = c(":", "-"),fragments = fragpath,annotation = annotation)

## Quality control
DefaultAssay(myObject) <- "ATAC"
myObject <- NucleosomeSignal(myObject)

head(Annotation(myObject))
head(Fragments(myObject)[[1]])

myObject <- TSSEnrichment(myObject)

DefaultAssay(myObject) <- "RNA"
myObject[["percent.mt"]] <- PercentageFeatureSet(myObject, pattern = "^mt-")

### Check cell quality
figure_name <- ""
figure_name <- paste(mysample, "_QC_vlnplot.pdf", sep="") 
pdf(file=figure_name)
VlnPlot(object = myObject, features = c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "percent.mt"), size.title.use =0.1,  pt.size=0.1, ncol = 6)
dev.off()

myObject <- subset(x = myObject,subset = nCount_ATAC < 100000 & nCount_RNA < 30000 & nCount_ATAC > 1000 & nCount_RNA > 1000 & nFeature_RNA > 500 & nucleosome_signal < 1.2 & TSS.enrichment > 2 & percent.mt < 15)


#myObject$pct_reads_in_peaks <- myObject$peak_region_fragments / myObject$passed_filters * 100
#myObject$blacklist_ratio <- myObject$blacklist_region_fragments / myObject$peak_region_fragments

figure_name <- ""
figure_name <- paste(mysample, "_ATAC_DensityScatter.pdf", sep="") 
pdf(file=figure_name)
DensityScatter(myObject, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off()

figure_name <- ""
figure_name <- paste(mysample, "_RNA_DensityScatter.pdf", sep="")
pdf(file=figure_name) 
DensityScatter(myObject, x ='nCount_RNA', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
dev.off() 

figure_name <- "" 
figure_name <- paste(mysample, "nucleosome_signalDensityScatter.pdf", sep="") 
pdf(file=figure_name) 
DensityScatter(myObject, x = 'nucleosome_signal' , y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE) 
dev.off() 

figure_name <- ""
figure_name <- paste(mysample, "_vlpot.pdf", sep="")
pdf(file="vlplot.pdf")
VlnPlot(
  object = myObject,
  features = c('nCount_RNA', 'nFeature_RNA', 'nCount_ATAC', 'TSS.enrichment', 'nucleosome_signal'),
  pt.size = 0.1,
  ncol = 5, 
  size.title.use = 0.1
)
dev.off()

figure_name <- ""
figure_name <- paste(mysample, "_FragmentHist.pdf", sep="")
pdf(file=figure_name)
myObject$nucleosome_group <- ifelse(myObject$nucleosome_signal > 4, 'NS > 4', 'NS < 4')
FragmentHistogram(object = myObject, group.by = 'nucleosome_group')
dev.off()

figure_name <- ""
figure_name <- paste(mysample, "_TSS.pdf", sep="")
pdf(file=figure_name)
myObject$high.tss <- ifelse(myObject$TSS.enrichment > 3, 'High', 'Low')
table(myObject$high.tss)
TSSPlot(myObject, group.by = 'High') + NoLegend()
dev.off()

saveRDS(myObject, file = myRDS)



