
### read in arguments
args <- commandArgs(TRUE)

# number of cells for testing
cells <- args[1]
# how many iterations
iteration <- args[2]
# output name
outName <- args[3]
# working directory
wd<-args[4]

setwd(wd)

print(cells)

######## load libraries ###########
# for FCS
library(flowCore)

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
sourceCpp(file = "jaccard_coeff.cpp")

# PARC
library(reticulate)
use_python('/usr/bin/python3',required = TRUE)
parc = import("parc")

############# load gold standard data #################
## Levine 13
ldf<-as.data.frame(exprs(read.FCS("gs/Levine_13dim.fcs")))
ldf<-ldf[which(!is.na(ldf$label)),]
lab_l13<-ldf$label
Levine13<-ldf[,1:13]
# Levine 32
ldf<-as.data.frame(exprs(read.FCS("gs/Levine_32dim.fcs")))
ldf<-ldf[which(!is.na(ldf$label)),]
lab_l32<-ldf$label
Levine32<-ldf[,5:36]
# Samusik
ldf<-as.data.frame(exprs(read.FCS("gs/Samusik_all.fcs")))
ldf<-ldf[which(!is.na(ldf$label)),]
lab_samAll<-ldf$label
SamusikAll<-ldf[,10:47]
# Samusik 01
ldf<-as.data.frame(exprs(read.FCS("gs/Samusik_01.fcs")))
ldf<-ldf[which(!is.na(ldf$label)),]
lab_sam01<-ldf$label
Samusik01<-ldf[,10:47]
gs<-c("Levine13","Levine32","SamusikAll","Samusik01")
list <- lapply(gs, get)
lab<-c("lab_l13","lab_l32","lab_samAll","lab_sam01")
list_lab <- lapply(lab, get)

########## run comparisons ###########
results<- NULL
count=0
for(model in list){
  count=count+1
  for(i in 1:iteration){
    model<-as.data.frame(model)
    s<-sample(1:dim(model)[1],cells)
    orig<-as.data.frame(list_lab[count])
    colnames(orig)<-"label"
    origs<-orig[s,]
    dat<-as.matrix(model[s,])
    
    # fastPG
    cluster<-fastCluster(dat,k = 30, num_threads = 4)
    nms<-cluster$communities
    fpgf<-FMeasure(origs,nms)
    
    # phenoGraph
    neighborMatrix <- find_neighbors(dat, k=30)[,-1]
    links <- jaccard_coeff(neighborMatrix)
    links <- links[links[,1]>0, ]
    relations <- as.data.frame(links)
    colnames(relations)<- c("from","to","weight")
    g <- graph.data.frame(relations, directed=FALSE)
    community <- cluster_louvain(g)
    bd<-community$membership
    pgf<-FMeasure(origs,bd)
    
    # PARC
    parc1<-parc$PARC(dat)
    parc1$run_PARC()
    parc1_labels <- unlist(parc1$labels)
    pf<-FMeasure(origs,parc1_labels)
    
    # flowSOM
    l<-flowFrame(exprs=as.matrix(dat))
    out <- FlowSOM::ReadInput(l, transform = FALSE, scale = FALSE)
    out <- FlowSOM::BuildSOM(out)
    out <- FlowSOM::BuildMST(out)
    labels_pre <- out$map$mapping[, 1]
    out <- FlowSOM::metaClustering_consensus(out$map$codes) 
    labels <- out[labels_pre]
    fsf<-FMeasure(origs,labels)
    
    
    # record results
    results<-rbind(results,data.frame(gs[count],fpgf,pgf,pf,fsf))
  }
}

write.table(results,outName,sep="\t",quote = F,row.names = F)
