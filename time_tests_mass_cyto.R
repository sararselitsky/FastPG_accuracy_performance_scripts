# Author: Sara Selitsky
# SaraRSelitsky@gmail.com

######## libraries ###########
# for FCS
library(flowCore)
library(data.table)

# metrics
library(mclust)
library(FlowSOM)

# for fastPG
library(FastPG)

# PhenoGraph
library(igraph)
library(RANN)
library(Rcpp)
find_neighbors <- function(data, k){
  nearest <- nn2(data, data, k, treetype = "kd", searchtype = "standard")
  return(nearest[[1]])
}
sourceCpp(file = "/home/Sara/jaccard_coeff.cpp")

# PARC
library(reticulate)
use_python('/usr/bin/python3',required = TRUE)
parc = import("parc")

# load data #

# fr_fcm_zyt6
fr<-fread("/home/Sara/time_tests/fr_fcm_zyt6_combined.txt",header = T)
fr<-as.data.frame(fr)

# test times for 10K to 5M cells
results<- NULL
size<-c(10000,50000,100000,1000000,5000000)

for(i in 1:length(size)){
  s<-sample(1:dim(fr)[1],size[i],replace=TRUE)
  dat<-as.matrix(fr[s,])
  # phenoGraph
  all_start = Sys.time()
  start = Sys.time()
  neighborMatrix <- find_neighbors(dat, k=30)[,-1]
  links <- jaccard_coeff(neighborMatrix)
  links <- links[links[,1]>0, ]
  relations <- as.data.frame(links)
  colnames(relations)<- c("from","to","weight")
  g <- graph.data.frame(relations, directed=FALSE)
  community <- cluster_louvain(g)
  bd<-community$membership
  end = Sys.time()
  pgTime = end - start
  pgTime<-as.numeric(pgTime,units="secs")
  
  # PARC
  all_start = Sys.time()
  start = Sys.time()
  parc1<-parc$PARC(dat)
  parc1$run_PARC()
  parc1_labels <- unlist(parc1$labels)
  end = Sys.time()
  pTime = end - start
  pTime<-as.numeric(pTime,units="secs")
  
  # flowSOM
  all_start = Sys.time()
  start = Sys.time()
  l<-flowFrame(exprs=as.matrix(dat))
  out <- FlowSOM::ReadInput(l, transform = FALSE, scale = FALSE)
  out <- FlowSOM::BuildSOM(out)
  out <- FlowSOM::BuildMST(out)
  labels_pre <- out$map$mapping[, 1]
  out <- FlowSOM::metaClustering_consensus(out$map$codes) 
  labels <- out[labels_pre]
  end = Sys.time()
  fsTime = end - start
  fsTime<-as.numeric(fsTime,units="secs")
  
  # fastPG
  all_start = Sys.time()
  start = Sys.time()
  cluster<-fastCluster(dat,k = 30, num_threads = 4)
  nms<-cluster$communities
  end = Sys.time()
  fpgTime = end - start
  fpgTime<-as.numeric(fpgTime,units="secs")
  
  # record results
  results<-rbind(results,data.frame(size[i],pgTime,fpgTime,pTime,fsTime))
}

write.table(results,"/home/Sara/time_tests/results/time_results.txt",sep="\t",quote = F,row.names = F)
