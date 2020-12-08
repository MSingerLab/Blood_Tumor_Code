bootstrap_auc <- function(resol, markers, neg){
  auc = vector(mode="numeric", length=length(markers))
  md$X = factor(md$X, levels=as.character(colnames(marker_data)))
  md = md[order(md$X),]
  aucs = foreach (i=1:length(markers), .combine = c, .inorder=TRUE) %dopar% {
    select = sample(1:nrow(md), nrow(md), replace=TRUE)
    mds= md[select,]
    mds = add_pass(mds, markers[i])
    col = which(colnames(mds)==markers[i])
    mds$marker = mds[,col]
    md_match = mds[mds$TM==1,]
    md_nomatch = mds[mds$TM!=1,]
    dist = c(seq(min(mds[,col]), max(mds[,col]), length=(1/resol)), as.numeric(quantile(mds[,col], probs=seq(resol, (1-resol), by=resol))))
    ax = sapply(1:length(dist), function(x) 1- (nrow(md_nomatch[md_nomatch$marker>dist[x],])/nrow(md_nomatch)))
    if (neg[i]==FALSE){
      ax = 1 - ax
    }
    ay = sapply(1:length(dist), function(x) nrow(md_match[md_match$marker<=dist[x],])/nrow(md_match))
    if (neg[i]==FALSE){
      ay = 1 - ay
    }
    ax = c(ax, 0, 1)
    ay = c(ay, 0, 1)
    auch = AUC(ax, ay)
    md_null = mds
    md_null$TM = sample(md_null$TM, nrow(md_null), replace=FALSE)
    md_match = md_null[md_null$TM==1,]
    md_nomatch = md_null[md_null$TM!=1,]
    nx = sapply(1:length(dist), function(x) 1- (nrow(md_nomatch[md_nomatch$marker>dist[x],])/nrow(md_nomatch)))
    if (neg[i]==FALSE){
      nx = 1 - nx
    }
    ny = sapply(1:length(dist), function(x) nrow(md_match[md_match$marker<=dist[x],])/nrow(md_match))
    if (neg[i]==FALSE){
      ny = 1 - ny
    }
    nx = c(0, nx, 1)
    ny = c(0, ny, 1)
    aucn = AUC(nx, ny)
    c(auch, aucn)
  }
  return(aucs)
}

library(doParallel)
library(doSNOW)
library(DescTools)
library(dplyr)
library(Seurat)

source("bootstrap_sens_spec_1128.R")

pack = c("doParallel", "doSNOW", "DescTools", "dplyr", "Seurat")
samples = c("k409c1", "integrated", "k411Ac1", "k411Bc1", "k468Ac1", "k468Bc1", "k484c1")

set.seed(27)
args = commandArgs(trailingOnly=TRUE)
id = samples[as.numeric(as.character(args[1]))]
md <- readRDS(paste(id, "blood_metadata.rds", sep="_"))
marker_data <- readRDS(paste(id, "blood_markers.rds", sep="_"))
markers = c("CCR7", "GYPC", "FLT3LG", "LTB")
for (i in 1:length(markers)){
  ind = which(rownames(marker_data)==markers[i])
  data = as.numeric(as.character(marker_data[ind,]))
  md <- data.frame(cbind(md,data))
  colnames(md)[ncol(md)] = markers[i]
}
neg = rep(TRUE, 4)
#markers = c(markers, "combo", "ands", "ors")
cl <- parallel::makeCluster(parallel::detectCores())
doParallel::registerDoParallel(cl)
parallel::clusterExport(cl, c("marker_data", "md", "markers", "pack"))
df = foreach (i=1:10000, .combine=rbind, .packages=pack, .inorder=FALSE) %dopar% {
  aucs = bootstrap_auc(0.05, markers, neg)
}
cnm = vector(mode="character")
for (i in 1:length(markers)){
  cnm[(2*i)-1] = markers[i]
  cnm[(2*i)] = paste(markers[i], "null", sep="_")
}
colnames(df) = cnm
write.csv(df, paste(id, "aucs_bootstrapped.csv", sep="_"))