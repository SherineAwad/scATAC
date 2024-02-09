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

annotation <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79, standard.chromosomes = TRUE)
seqlevelsStyle(annotation) <- "UCSC"
unique(seqnames(annotation)) 

###load sample, the RNA and ATAC data
counts <- Read10X_h5(myh5)
# Check counts
lapply(counts, dim)
lapply(counts, class)

# create a Seurat object containing the RNA adata
myObject <- CreateSeuratObject(counts = counts$`Gene Expression`,assay = "RNA", project=mysample)

# create ATAC assay and add it to the object
myObject[["ATAC"]] <- CreateChromatinAssay(counts = counts$Peaks,sep = c(":", "-"),fragments = mytsv,annotation = annotation)



## Quality control

DefaultAssay(myObject) <- "ATAC"
myObject <- NucleosomeSignal(myObject)

head(Annotation(myObject))
head(Fragments(myObject)[[1]])

myObject <- TSSEnrichment(myObject)

myObject[["percent.mt"]] <- PercentageFeatureSet(myObject, pattern = "^mt-", assay="RNA") 
### Check cell quality
figure_name <- ""
figure_name <- paste(mysample, "_QC_vlnplot.pdf", sep="")
pdf(file=figure_name, width=12)
VlnPlot(object = myObject, features = c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "percent.mt"),pt.size=0.1, ncol = 6)
dev.off()

myRDS <- paste(mysample, "_preprocessed.rds", sep="")
saveRDS(myObject, file = myRDS)





