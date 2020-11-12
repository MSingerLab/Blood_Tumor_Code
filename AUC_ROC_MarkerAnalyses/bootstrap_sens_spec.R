
# Load libraries
library(Seurat)
library(DescTools)
library(doParallel)
library(doSNOW)
library(dplyr)

# Load the function 'Keep_AB_Cells' which subsets the rds object
# to keep only cells with at least one alpha and one beta chain
source("Keep_AB_Cells.R")

add_pass <- function(md, gate){
  # md = metadata
  md$pass = 0
  if (gate=="combo"){ #this is the combinatorial gate previously selected to have the best (balanced) performance
    md$pass[md$LTB <= 0.001] = 1
    md$subpass = 1
    md$subpass[md$CCR7 > 0.001] = 0
    md$subpass[md$FLT3LG > 0.001] = 0
    md$subpass[md$GYPC > 0.001] = 0
    md$pass[md$subpass == 1] = 1
    md$subpass = NULL
  } else if (gate=="ands"){
    md$pass = 1
    md$pass[md$CCR7 > 0.001] = 0
    md$pass[md$FLT3LG > 0.001] = 0
    md$pass[md$GYPC > 0.001] = 0
    md$pass[md$LTB > 0.001] = 0
  } else if (gate=="ors"){
    md$pass[md$CCR7 <= 0.001] = 1
    md$pass[md$FLT3LG <= 0.001] = 1
    md$pass[md$GYPC <= 0.001] = 1
    md$pass[md$LTB <= 0.001] = 1
  } else if (gate=="CCR7"){
    md$pass[md$CCR7 <= 0.001] = 1
  } else if (gate=="FLT3LG"){
    md$pass[md$FLT3LG <= 0.001] = 1
  } else if (gate=="GYPC"){
    md$pass[md$GYPC <= 0.001] = 1
  } else if (gate=="LTB"){
    md$pass[md$LTB <= 0.001] = 1
  }
  return(md)
}

bootstrap_ss <- function(markers){
  md$Barcode = factor(md$Barcode, levels=as.character(colnames(marker_data)))
  md = md[order(md$Barcode),]
  md$marker = 0
  spec_sens = foreach (i=1:length(markers), .combine = c, .inorder=TRUE) %dopar% {
	select = sample(1:nrow(md), nrow(md), replace=TRUE)
	mds= md[select,]
    	mds = add_pass(mds, markers[i])
	md_match = mds[mds$TM==1,]
    	md_nomatch = mds[mds$TM!=1,]
	sens = nrow(md_match[md_match$pass == 1,])/nrow(md_match)
	spec = nrow(md_nomatch[md_nomatch$pass == 0,])/nrow(md_nomatch)	
	md_null = mds
	md_null$TM = sample(md_null$TM, nrow(md_null), replace=FALSE)
	md_match = md_null[md_null$TM==1,]
    	md_nomatch = md_null[md_null$TM!=1,]
	nsens = nrow(md_match[md_match$pass == 1,])/nrow(md_match)
	nspec = nrow(md_nomatch[md_nomatch$pass == 0,])/nrow(md_nomatch)
	c(sens, spec, nsens, nspec)               
  }
  return(spec_sens)
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
 markers = c(markers, "combo", "ors", "ands")
  cl <- parallel::makeCluster(parallel::detectCores())
  doParallel::registerDoParallel(cl)
  parallel::clusterExport(cl, c("marker_data", "md", "markers", "pack"))
  df = foreach (i=1:10000, .combine=rbind, .packages=pack, .inorder=FALSE) %dopar% {
    values = bootstrap_ss(markers)
  }
  cnm = vector(mode="character")
  for (i in 1:length(markers)){
    cnm[(4*i)-3] = paste(markers[i], "sens", sep="_")
    cnm[(4*i)-2] = paste(markers[i], "spec", sep="_")
    cnm[(4*i)-1] = paste(markers[i], "sens_null", sep="_")
    cnm[(4*i)] = paste(markers[i], "spec_null", sep="_")
  } 
  colnames(df) = cnm
  write.csv(df, paste0('outputs/', id, "_specsens_bootstrapped.csv"))
}
