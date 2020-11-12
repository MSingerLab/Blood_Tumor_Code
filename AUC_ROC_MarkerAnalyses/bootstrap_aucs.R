
# Load libraries
library(Seurat)
library(DescTools)
library(doParallel)
library(doSNOW)
library(dplyr)

# Load the function 'Keep_AB_Cells' which subsets the rds object
# to keep only cells with at least one alpha and one beta chain
source("Keep_AB_Cells.R")
source("bootstrap_sens_spec.R")


bootstrap_auc <- function(resol, markers, neg){
  print('Started function')
  auc = vector(mode="numeric", length=length(markers))
  md$Barcode = factor(md$Barcode, levels=as.character(colnames(marker_data)))
  md = md[order(md$Barcode),]
  aucs = foreach (i=1:length(markers), .combine = c, .inorder=TRUE) %dopar% {
    marker = markers[i]
	  select = sample(1:nrow(md), nrow(md), replace=TRUE)
	  mds = md[select,]
    mds = add_pass(mds, marker)
	  md_match = mds[mds$TM==1,]
    md_nomatch = mds[mds$TM!=1,]
	  col = which(colnames(md)==markers[i])
	  dist = c(seq(min(mds[,col]), max(mds[,col]), length=(1/resol)), as.numeric(quantile(mds[,col], probs=seq(resol, (1-resol), by=resol))))
	x = foreach (j=1:length(dist), .combine=c, .inorder=TRUE) %dopar% {
      		xj = 1 - (nrow(md_nomatch[md_nomatch[, col]>dist[j],])/nrow(md_nomatch))
      		if (neg[i]==FALSE){
			xj = 1 - xj
      		}
		xj
    	}
    	y = foreach (j=1:length(dist), .combine=c, .inorder=TRUE) %dopar% {
      		yj = nrow(md_match[md_match[, col]<=dist[j],])/nrow(md_match)
      		if (neg[i]==FALSE){
			yj = 1 - yj
      		}
		yj
    	}
	x = c(x, 0, 1)
    	y = c(y, 0, 1)
    	auch = AUC(x, y)
	md_null = mds
	md_null$TM = sample(md_null$TM, nrow(md_null), replace=FALSE)
	md_match = md_null[md_null$TM==1,]
    	md_nomatch = md_null[md_null$TM!=1,]
	nx = foreach (j=1:length(dist), .combine=c, .inorder=TRUE) %dopar% {
      		xj = 1 - (nrow(md_nomatch[md_nomatch$marker>dist[j],])/nrow(md_nomatch))
      		if (neg[i]==FALSE){
			xj = 1 - xj
      		}
		xj
    	}
    	ny = foreach (j=1:length(dist), .combine=c, .inorder=TRUE) %dopar% {
      		yj = nrow(md_match[md_match$marker<=dist[j],])/nrow(md_match)
      		if (neg[i]==FALSE){
			yj = 1 - yj
      		}
		yj
    	}
        nx = c(0, nx, 1)
        ny = c(0, ny, 1)
	aucn = AUC(nx, ny)
	c(auch, aucn)               
  }
  return(aucs)
}



pack = c("doParallel", "doSNOW", "DescTools", "dplyr", "Seurat")
samples = c("K409_Blood",
            "K411_Blood", "K411_Blood_Longitudinal",
            "K468_Blood", "K468_Blood_Longitudinal",
            "K484_Blood", "HumanIntegratedBlood_wLongitudinal")

set.seed(27)
for (s in 1:length(samples)){
  id = samples[s]
  
  print(id)
  
  # Load data, keep only AB cells
  rds <- Keep_AB_Cells(paste0('data/',  id, ".rds"))
  
  # Extract metadata and count matrix
  md <- rds@meta.data
  marker_data <- rds@assays$RNA@counts
  

  markers = c("CCR7", "GYPC", "FLT3LG", "LTB")
  for (i in 1:length(markers)){
    ind = which(rownames(marker_data)==markers[i])
    data = as.numeric(as.character(marker_data[ind,]))
    md <- data.frame(cbind(md,data))
    colnames(md)[ncol(md)] = markers[i]
  }
 neg = rep(TRUE, 4)
  cl <- parallel::makeCluster(parallel::detectCores())
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl, c("marker_data", "md", "markers", "pack"))
  
  df = foreach (i=1:10, .combine=rbind, .packages=pack, .inorder=FALSE) %dopar% {
    aucs = bootstrap_auc(0.05, markers, neg)
  }
  cnm = vector(mode="character")
  for (i in 1:length(markers)){
    cnm[(2*i)-1] = markers[i]
    cnm[(2*i)] = paste(markers[i], "null", sep="_")
  } 
  colnames(df) = cnm
  write.csv(df, paste0('outputs/', id, "_aucs_bootstrapped.csv"))
}

