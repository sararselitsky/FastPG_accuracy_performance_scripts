######################
# Subsampling cells
######################

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

############################################################
# set seeds for reproducibility
myseed <- 1234
set.seed(myseed)

# set # of cells to subsample
size <- 2000

# subsample 10 times
subcol.list <- list()
for (i in 1: 10) {
  print(i)
  subcol.list[[i]] <- sample(c(1:nrow(label)), 2000, replace=FALSE)
}
subcol.list

#############################################################
library(Seurat)

# set # of PCs
ndim <- 30
# set Seurat resolution parameter
resol <- 0.8

# # of iterations
rep <- 10

# create empty lists to store output
ep <- list() # for the PCs
meta.data.list <- list() # for Seurat labels

for (i in 1: rep) {
  print("rep")
  tmp <- count[, subcol.list[[i]] ]
  seu <- CreateSeuratObject(tmp,
                            meta.data = label, 
                            min.cells=1,
                            min.features=1,
                            project=fname)
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
  ep[[i]] <- seu@reductions$pca@cell.embeddings[, 1: ndim]
  
  
  # Run Clustering
  seu <- FindNeighbors(seu,
                       reduction = "pca", # default is pca
                       k.param = 20, # default is 20
                       # default is 10 dimensions
                       dims = 1: ndim) 
  seu <- FindClusters(seu, resolution = resol)
  
  meta.data.list[[i]] <- seu@meta.data
}

# save output
saveRDS(subcol.list, file=paste(fname, ".2k.30PC.list.rds"), sep="")
saveRDS(meta.data.list, file=paste(fname, ".2k.30PC.Seu.meta.data.list.rds", sep=""))


