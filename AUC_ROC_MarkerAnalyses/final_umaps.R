
# Load libraries
library(Seurat)
library(ggplot2)

# Load the function 'Keep_AB_Cells' which subsets the rds object
# to keep only cells with at least one alpha and one beta chain
source("Keep_AB_Cells.R")

# Define functions
processCometFile <- function(id){
  d = read.csv(paste(id, "singletons_HS_full.csv", sep="_"))
  d$Sample = id
  return(d)
}

set_up_cutoffs <- function(markers, samples){
  cutoffs = data.frame()
  if (all(samples == c("M4_Blood", "M5_Blood"))){
    for (i in 1:length(markers)){
      for (j in 1:length(samples)){
        comet_data = read.csv(paste0('data/', samples[j], "_singletons_full.csv"))
        ind = which(comet_data$gene_1 == markers[i])
        cutoffs[j,i] = max(0.001, comet_data$cutoff_val[ind])
      }
    }
  } else {
    for (i in 1:length(markers)){
      for (j in 1:length(samples)){
        comet_data = read.csv(paste0('data/', samples[length(samples)], "_singletons_HS_full.csv"))
        ind = which(comet_data$gene_1 == markers[i])
        cutoffs[j,i] = max(0.001, comet_data$cutoff_val[ind])
      }
    }
  }
  return(cutoffs)
}

# Define samples, markers for CITE-seq mice
samples = c("M4_Blood", "M5_Blood")
markers = c("NKG2D_TotalSeqC", "CX3CR1_TotalSeqC", "CD39_TotalSeqC")
cutoffs = set_up_cutoffs(markers, samples)

# For each mouse,
for (s in 1:length(samples)){
  
  # Load Seurat object (so), metadata (md)
  so = Keep_AB_Cells(paste0('data/', samples[s], '.rds'))
  md = so@meta.data
  md = md[order(md$TM),]
  
  # Plot tumor-matching and non-tumor matching cells
  pdf(paste0('outputs/', samples[s], "_answerkey_UMAP.pdf"), height=12, width=12)
  g = ggplot(md, aes(UMAP_1, UMAP_2, color=TM))
  g = g + geom_point(size = 1) + scale_color_manual(values=c("darkgreen", "gray60")) + theme_classic()
  print(g)
  dev.off()
  
  # Plot UMAP for each marker after binarizing counts
  for (m in 1:length(markers)){
    col = which(colnames(md)==markers[m])
    md$pass = 0
    for (i in 1:nrow(md)){
      if (md[i,col]  >= cutoffs[s,m] ) { md$pass[i] = 1 }
    }
    md$pass = as.factor(md$pass)
    md = md[order(md$pass),]
    pdf(paste0('outputs/', samples[s], '_', markers[m], "_UMAP.pdf"), height=12, width=12)
    g = ggplot(md, aes(UMAP_1, UMAP_2, color=pass))
    g = g + geom_point(size = 1) + scale_color_manual(values=c("red", "gray60")) + theme_classic() 
    print(g)
    dev.off()
  }
}

# Now work on integrated object
so = Keep_AB_Cells('data/HumanIntegratedBlood_wLongitudinal.rds')
md = so@meta.data
counts = so@assays$RNA

samples = c('K409_Blood', 'K411_Blood', 'K411_Blood_Longitudinal', 'K468_Blood', 'K468_Blood_Longitudinal', 'K484_Blood')

setwd('data/')
all_comet_results <- processCometFile(samples[1])
for (i in 2:length(samples)){
  all_comet_results = data.frame(rbind(all_comet_results, processCometFile(samples[i])))
}
setwd('..')

markers = c("LTB_negation", "GYPC_negation", "FLT3LG_negation", "CCR7_negation")
markers = sort(markers)
cd = all_comet_results[all_comet_results$gene_1 %in% markers,]
cd$Sample = as.character(cd$Sample)
for (i in 1:nrow(cd)){
  cd$cutoff_val[i] = max(0.001, cd$cutoff_val[i])
}


l = which(rownames(counts)=="LTB")
f = which(rownames(counts)=="FLT3LG")
c = which(rownames(counts)=="CCR7")
g = which(rownames(counts)=="GYPC")
p = which(rownames(counts)=="PDCD1")

md$LTB = as.vector(as.numeric(as.character(counts[l,])))
md$FLT3LG = as.vector(as.numeric(as.character(counts[f,])))
md$CCR7 = as.vector(as.numeric(as.character(counts[c,])))
md$GYPC = as.vector(as.numeric(as.character(counts[g,])))
md$PDCD1 = as.vector(as.numeric(as.character(counts[p,])))

md$TM = as.factor(md$TM)
md = md[order(md$TM),]
pdf("outputs/HumanIntBlood_answerkey_UMAP.pdf", height=12, width=12)
g = ggplot(md, aes(UMAP_1, UMAP_2, color=TM))
g = g + geom_point(size = 1) + scale_color_manual(values=c("darkgreen", "gray60")) + theme_classic() 
print(g)
dev.off()

##f&l or c ot t

md$gate3 = 0

for (i in 1:nrow(md)){
  d = cd[cd$Sample == md$Sample[i],]
  cutoffs = as.numeric(as.character(d$cutoff_val))
  res = 0
  if (md$LTB[i] <= cutoffs[4]) { res = 1 }
  subres = 1
  if (md$GYPC[i] > cutoffs[3]) { subres = 0 }
  if (md$CCR7[i] > cutoffs[1]) { subres = 0 }
  if (md$FLT3LG[i] > cutoffs[2]) { subres = 0 }
  if (subres==1) { res = 1 }
  md$gate3[i] = res
}

md$gate3 = as.factor(md$gate3)
md = md[order(md$gate3),]
pdf("outputs/Fig6H_combo_final_UMAP.pdf", height=8, width=8)
g = ggplot(md, aes(UMAP_1, UMAP_2, color=gate3))
g = g + geom_point(size = 0.5) + scale_color_manual(values=c("gray60", "blue")) + theme_classic() 
print(g)
dev.off()

md$gateC = 0

for (i in 1:nrow(md)){
  d = cd[cd$Sample == md$Sample[i],]
  cutoffs = as.numeric(as.character(d$cutoff_val))
  res = 0
  if (md$LTB[i] <= cutoffs[4]) { res = 1 }
  md$gateC[i] = res
}

md$gateC = as.factor(md$gateC)
md = md[order(md$gateC),]
pdf("outputs/SuppFig4L_LTB_UMAP.pdf", height=8, width=8)
g = ggplot(md, aes(UMAP_1, UMAP_2, color=gateC))
g = g + geom_point(size = 0.5) + scale_color_manual(values=c("gray60", "blue")) + theme_classic() 
print(g)
dev.off()

md$gateC = 0

for (i in 1:nrow(md)){
  d = cd[cd$Sample == md$Sample[i],]
  cutoffs = as.numeric(as.character(d$cutoff_val))
  res = 0
  if (md$GYPC[i] <= cutoffs[3]) { res = 1 }
  md$gateC[i] = res
}

md$gateC = as.factor(md$gateC)
md = md[order(md$gateC),]
pdf("outputs/SuppFi4L_GYPC_UMAP.pdf", height=8, width=8)
g = ggplot(md, aes(UMAP_1, UMAP_2, color=gateC))
g = g + geom_point(size = 0.5) + scale_color_manual(values=c("gray60", "blue")) + theme_classic() 
print(g)
dev.off()

md$gateC = 0

for (i in 1:nrow(md)){
  d = cd[cd$Sample == md$Sample[i],]
  cutoffs = as.numeric(as.character(d$cutoff_val))
  res = 0
  if (md$FLT3LG[i] <= cutoffs[2]) { res = 1 }
  md$gateC[i] = res
}

md$gateC = as.factor(md$gateC)
md = md[order(md$gateC),]
pdf("outputs/SuppFig4L_FLT3LG_UMAP.pdf", height=8, width=8)
g = ggplot(md, aes(UMAP_1, UMAP_2, color=gateC))
g = g + geom_point(size = 0.5) + scale_color_manual(values=c("gray60", "blue")) + theme_classic() 
print(g)
dev.off()

md$gateC = 0

for (i in 1:nrow(md)){
  d = cd[cd$Sample == md$Sample[i],]
  cutoffs = as.numeric(as.character(d$cutoff_val))
  res = 0
  if (md$CCR7[i] <= cutoffs[1]) { res = 1 }
  md$gateC[i] = res
}

md$gateC = as.factor(md$gateC)
md = md[order(md$gateC),]ß
pdf("outputs/SuppFig4L_CCR7_UMAP.pdf", height=8, width=8)
g = ggplot(md, aes(UMAP_1, UMAP_2, color=gateC))
g = g + geom_point(size = 0.5) + scale_color_manual(values=c("gray60", "blue")) + theme_classic() 
print(g)
dev.off()

