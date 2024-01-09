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

#install.packages("Matrix", type = "source",repos = "http://cran.us.r-project.org")
#install.packages("irlba", type = "source",repos = "http://cran.us.r-project.org")

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
genome(annotation) <- "mm10"
# Check annotation seqnames
unique(seqnames(annotation))

#edb <- EnsDb.Mmusculus.v79
#seqlevelsStyle(edb) <- "UCSC"
#annotation <- GetGRangesFromEnsDb(edb)

###load sample, the RNA and ATAC data
counts <- Read10X_h5(myh5)
#No need to print a long list now 
#head(counts)
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

DefaultAssay(myObject) <- "RNA"
myObject[["percent.mt"]] <- PercentageFeatureSet(myObject, pattern = "^mt-")

### Check cell quality
figure_name <- ""
figure_name <- paste(mysample, "_QC_vlnplot.pdf", sep="") 
pdf(file=figure_name, width=12)
VlnPlot(object = myObject, features = c("nCount_RNA", "nFeature_RNA", "nCount_ATAC", "TSS.enrichment", "nucleosome_signal", "percent.mt"),pt.size=0.1, ncol = 6)
dev.off()

myObject <- subset(x = myObject,subset = nCount_ATAC < 100000 & nCount_RNA < 30000 & nCount_ATAC > 1000 & nCount_RNA > 1000 & nFeature_RNA > 500 & nucleosome_signal < 1.2 & TSS.enrichment > 2 & percent.mt < 15)


#### Analyze RNA part

DefaultAssay(myObject) <- "RNA"

myObject <- SCTransform(myObject, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:30, reduction.name = 'umap.rna', reduction.key = 'UMAP_')

saveRDS(myObject, file = myRDS)


