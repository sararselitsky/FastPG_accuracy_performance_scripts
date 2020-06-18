
######## libraries ###########
# FastPG
library(reticulate)
use_python('/usr/bin/python3',required = TRUE)
library(FastPG)

# Seurat
library(Seurat)

# read in data
ldf_p68<-readRDS(file="/home/Sara/pbmc68k.30PC.rds")

# test times of FastPG and Seurat for 1000 to 68K cells
results<- NULL
size<-c(1000,5000,10000,50000,68000)

for(i in 1:length(size)){
  s<-sample(1:dim(ldf_p68)[1],size[i],replace=FALSE)
  dat<-as.matrix(ldf_p68[s,])
  
  # fastPG
  all_start = Sys.time()
  start = Sys.time()
  cluster<-fastCluster(dat,k = 30, num_threads = 4)
  nms<-cluster$communities
  end = Sys.time()
  fpgTime = end - start
  fpgTime<-as.numeric(fpgTime,units="secs")
  
  all_start = Sys.time()
  start = Sys.time()
  seu <- FindNeighbors(dat, k.param = 30,reduction = "pca") 
  seu <- FindClusters(seu$snn,resolution = 0.8)
  seu<-NULL
  end = Sys.time()
  sTime = end - start
  sTime<-as.numeric(sTime,units="secs")
  
  # record results
  results<-rbind(results,data.frame(size[i],fpgTime,sTime))
}

write.table(results,"/home/Sara/time_tests/results/time_results_scRNA_seq.txt",sep="\t",quote = F,row.names = F)
