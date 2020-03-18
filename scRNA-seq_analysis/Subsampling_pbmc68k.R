######################
# Subsampling cells
######################

# pbmc68k dataset
library(Seurat)
mydata <- Read10X(data.dir="/data/hg19/")
fname <- "pbmc68k"

# load cell labels
label <- read.table("zheng17-PBMC68K-cell-labels.txt", header=T, sep="\t")

##############################
# set seed for reproducibility
myseed <- 311
set.seed(myseed)

size <- 10000
size <- 50000

subcol <- sample(c(1: ncol(mydata)), size, replace=TRUE)
head(subcol)

subdata <- mydata[, subcol]
dim(subdata)

#########################################################
library(Seurat)
seu <- CreateSeuratObject(counts=subdata, project="pbmc68k", 
                          min.cells = 1, min.features = 1)


# Visualize QC metrics as a violin plot
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# Filter cells based on QC
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)


# Normalization
seu <- NormalizeData(seu, verbose=FALSE, normalization.method = "LogNormalize", scale.factor = 10000)

# Scaling
all.genes <- rownames(seu)
seu <- ScaleData(seu, features=all.genes)

# Find HVG for dim reduction
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)


# PCA
seu <- RunPCA(seu, features = VariableFeatures(object = seu))

# Run Clustering
ndim <- 30
resol <- 0.8

# create an empty vector to save the run time
seu_speed <- c()
for (i in 1: 10) {
  print(i)
  seu_speed[i] <- system.time( {
    seu <- FindNeighbors(seu,
                         reduction = "pca", # default is pca
                         k.param = 20, # default is 20
                         # default is 10 dimensions
                         dims = 1: ndim) ;
    seu <- FindClusters(seu, resolution = resol)
  })[3]
}
seu_speed

# get PC matrix
ep <- seu@reductions$pca@cell.embeddings[, 1: ndim]
saveRDS(ep, file=paste("pbmc68k.", nrow(seu@meta.data), "cells.30PC.rds", sep=""))


# save Seurat label
saveRDS(seu@meta.data, file=paste("pbmc68k.", "SeuLabel.", nrow(seu@meta.data), "cells.rds", sep=""))

# save Seurat runtime for the 10 runs
saveRDS(seu_speed, file=paste("pbmc68k.", nrow(seu@meta.data), "cells.30PC.Seu.runtime.10runs.rds", sep=""))




