# scRNA-seq datasets preproessing

# download Duo et al. datasets
#BiocManager::install("DuoClustering2018")
library(DuoClustering2018)
library(SingleCellExperiment)

# retrieve data
sce <- sce_full_Zhengmix4eq()
fname <- "Zhengmix4eq"

sce <- sce_full_Zhengmix8eq()
fname <- "Zhengmix8eq"

# extract count matrix
count <- counts(sce)

# extract labels
label <- as.data.frame(colData(sce)[c("barcode","phenoid")])

####################################################
# Run Seurat to do all the preprocessing steps

library(Seurat)

# Create Seurat object
seu <- CreateSeuratObject(count,
                          meta.data = label, 
                          min.cells=1,
                          min.features=1,
                          project=fname)

# Visualize QC metrics as a violin plot
VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)

# Filter cells based on QC
seu <- subset(seu, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)
seu

# Normalization
seu <- NormalizeData(seu, verbose=FALSE, normalization.method = "LogNormalize", scale.factor = 10000)

# Scaling
all.genes <- rownames(seu)
seu <- ScaleData(seu, features=all.genes)

# Find HVG for dim reduction
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

# PCA
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
DimPlot(seu, reduction = "pca")
ElbowPlot(seu, ndims = 50)

# Run Clustering
ndim <- 30
resol <- 0.8

system.time( {
  seu <- FindNeighbors(seu,
                       reduction = "pca", # default is pca
                       k.param = 20, # default is 20
                       # default is 10 dimensions
                       dims = 1: ndim) 
  seu <- FindClusters(seu, resolution = resol)
})

# save Seurat cell labels
saveRDS(seu@meta.data, file=paste(fname, ".Seurat.res0.8.label", ".rds", sep=""))

# get PC matrix for input data matrix for FastPG
ep <- seu@reductions$pca@cell.embeddings[, 1: ndim]
saveRDS(ep, file=paste(fname, ".30PC.rds", sep=""))




