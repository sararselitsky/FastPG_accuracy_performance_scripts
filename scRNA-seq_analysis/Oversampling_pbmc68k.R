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

size <- 100000
size <- 200000
# save the indices for the samples
overcol <- sample(c(1: ncol(mydata)), size, replace=TRUE)
# create data matrix for the oversampled samples
overdata <- mydata[, overcol]
# make oversampled matrix colnames unique
colnames(overdata) <- paste(paste0("o", c(1:size)), colnames(overdata), sep="-")


################################################################
library(Seurat)
seu <- CreateSeuratObject(counts=overdata, project="pbmc68k", 
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

seu_speed <- system.time( {
  seu <- FindNeighbors(seu,
                       reduction = "pca", # default is pca
                       k.param = 20, # default is 20
                       # default is 10 dimensions
                       dims = 1: ndim) 
  seu <- FindClusters(seu, resolution = resol)
}) [1:3]
seu_speed


# get PC matrix
ndim <- 30
ep <- seu@reductions$pca@cell.embeddings[, 1: ndim]
saveRDS(ep, file=paste("pbmc68k.", nrow(seu@meta.data), "cells.30PC.rds", sep=""))

# save Seurat labels
saveRDS(seu@meta.data, file=paste("pbmc68k.", "SeuLabel.", nrow(seu@meta.data), "cells.rds", sep=""))




