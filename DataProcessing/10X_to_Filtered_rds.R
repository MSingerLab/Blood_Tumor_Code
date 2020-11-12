#!/usr/bin/env Rscript
library(Seurat)
library(dplyr)
library(data.table)
library(ggplot2)
library(cowplot)
library(viridis)
library(gridExtra)
library(RColorBrewer)
library(tibble)

# This script takes the output from the 10X Cell Ranger, filters cells,
# clusters data, and outputs Seurat/rds objects

# Define samples
human_samples <- c("K409_Blood", "K411_Blood", "K411_Blood_Longitudinal",
                   "K468_Blood", "K468_Blood_Longitudinal", "K484_Blood",
                   "K409_Tumor", "K409_LN", "K411_Tumor", "K468_Tumor",
                   "K484_Tumor")

mouse_samples <- c("M1_Blood", "M1_Tumor",
                   "M2_Blood", "M2_Tumor",
                   "M3_Blood", "M3_Tumor",
                   "M4_Blood", "M4_Tumor",
                   "M5_Blood", "M5_Tumor")

# Set working directory
setwd("DataProcessing")

# Read the cellranger/10X output. The data directory has the following files:
# barcodes.tsv
# features.tsv
# matrix.mtx
data <- Read10X(data.dir = "filtered_feature_bc_matrix")

# Set the sample name
# Example:
sample <- "K409_Blood"

# Create a Seurat object 
data2 <- CreateSeuratObject(data, min.cells = 3, min.features = 400, project = sample)

# Rename a column
idx = which('orig.ident' == colnames(data2@meta.data))
colnames(data2@meta.data)[idx] <- "Sample"


# Save unfiltered RDS
file_name <- paste0(sample, '_unfiltered.rds')
saveRDS(data2, file = file_name)


# If human sample, load human housekeeping genes
if(sample %in% human_samples){
  housekeeping_genes <- readLines('GeneLists/HK_Satija_human.txt')
} else {
  # Otherwise load mouse housekeeping genes
  housekeeping_genes <- readLines('GeneLists/HK_Satija_mouse.txt')
}


# Intersection of genes
s1 <- rownames(data2)
s2 <- housekeeping_genes
housekeeping_genes <- s1[ toupper(s1) %in% toupper(s2) ]
subset_counts_matrix <- data.frame(data2@assays$RNA@counts)[housekeeping_genes,]
subset_counts_matrix[subset_counts_matrix > 0] = 1 
summary(colSums(subset_counts_matrix))
hist(colSums(subset_counts_matrix),main="> .5 selected as cutoff for HK genes")
pdf("HK_genes_hist.pdf",width=7,height=7)
hist(colSums(subset_counts_matrix),main="> .5 selected as cutoff for HK genes")
dev.off()
counts <- colSums(subset_counts_matrix)
counts2 <- data.frame(counts)
counts3 <- counts2 %>% tibble::rownames_to_column() %>% filter(counts2 > length(housekeeping_genes)/2)


all.genes <- rownames(x = data2)

# If human sample, load human mitochondrial genes
if(sample %in% human_samples){
  mito_genes <- readLines('GeneLists/mito_genes_human.txt', warn = FALSE)
} else {
  # Otherwise load mouse mitochondrial genes
  mito_genes <- readLines('GeneLists/mito_genes_mouse.txt', warn = FALSE)
}


# Intersection of mito genes
# The number of counts to be filtered is determined on a sample-to-sample basis
# If human, remove cells if they express more than 500 mito genes
# If mouse, remove cells if the mito gene expression was higher than 2 standard deviations from the mean
s1 <- rownames(data2)
s2 <- mito_genes
mito_genes <- s1[ toupper(s1) %in% toupper(s2) ]
subset_counts_matrix <- data.frame(data2@assays$RNA@counts)[mito_genes,]
subset_counts_matrix[subset_counts_matrix > 0] =1 
summary(colSums(subset_counts_matrix))
hist(colSums(subset_counts_matrix),main="< 500/1157 genes (Broad list) selected as cutoff for mito genes")
pdf("mito_genes_hist.pdf",width=7,height=7)
hist(colSums(subset_counts_matrix),main="< 500/1157 genes (Broad list) selected as cutoff for mito genes")
dev.off()

counts4 <- colSums(subset_counts_matrix)
counts5 <- data.frame(counts4)
counts6 <- counts5 %>% tibble::rownames_to_column() %>% filter(counts4 < 500) # Change this number

cells <- intersect(counts6$rowname,counts3$rowname)




# Normalize and scale data
working_so <- NormalizeData(data2)
working_so <- subset(working_so, cells = gsub('.', '-', cells, fixed = TRUE))
working_so <- ScaleData(object = working_so, features = rownames(working_so))

# Find variable genes
working_so <- FindVariableFeatures(object = working_so, mean.function = ExpMean, dispersion.function = LogVMR)
hv.genes <- head(x = VariableFeatures(object = working_so), 1000)

# Run PCA, find neighbors, clustering, and UMAP
working_so <- RunPCA(object = working_so, pc.genes = hv.genes, pcs.print = 1:5,
                     genes.print = 5, pcs.compute = 50)
working_so <- FindNeighbors(object = working_so, dims = 1:30)
working_so <- FindClusters(object = working_so, resolution = 1.2)
working_so <- RunUMAP(object = working_so, reduction = "pca", dims = 1:15, n_neighbors = 15, min_dist = 0.3)

# Plot the UMAP
DimPlot(object = working_so, reduction = 'umap')

# Add the UMAP coordinates to the metadata
working_so@meta.data$UMAP_1 <- working_so@reductions$umap@cell.embeddings[, 1]
working_so@meta.data$UMAP_2 <- working_so@reductions$umap@cell.embeddings[, 2]

# Optional, save rds object after mito/HK filtering
file_name <- paste0(sample, '_postMitoHK.rds')
saveRDS(working_so, file = file_name)

#optionally read in existing Seurat object
#working_so <- readRDS("K409_Blood.rds")


# Determine which clusters to keep 
clusters <- levels(working_so@meta.data$seurat_clusters)

# Iterative function to see if a given cluster should be kept, results aggregate like a fold in Haskell

# Build function for human data
see_if_keep_cluster_human <- function(so, cluster, idents){
  copy_so <- so
  cso <- SubsetData(copy_so, ident.use = cluster, do.clean = TRUE, do.scale = TRUE)
  new_id <- NULL
  cd3e <- sum(GetAssayData(object = cso, slot = "data")["CD3E",]>0)/nrow(cso@meta.data)
  cd3d <- sum(GetAssayData(object = cso, slot = "data")["CD3D",]>0)/nrow(cso@meta.data)
  cd3g <- sum(GetAssayData(object = cso, slot = "data")["CD3G",]>0)/nrow(cso@meta.data)
  cd8b <- sum(GetAssayData(object = cso, slot = "data")["CD8B",]>0)/nrow(cso@meta.data)
  cd8a <- sum(GetAssayData(object = cso, slot = "data")["CD8A",]>0)/nrow(cso@meta.data)
  foxp3 <- sum(GetAssayData(object = cso, slot = "data")["FOXP3",]>0)/nrow(cso@meta.data)
  cd4 <- sum(GetAssayData(object = cso, slot = "data")["CD4",]>0)/nrow(cso@meta.data)
  mki67 <- sum(GetAssayData(object = cso, slot = "data")["MKI67",]>0)/nrow(cso@meta.data)
  if((cd3e > .3 && cd3d > .3) || (cd3e > .3 && cd3g > .3) || (cd3d > .3 && cd3g > .3)){
    if (cd8b > .3 && cd8a > .3 && foxp3 < .05 && cd4 < .05) {
      new_id <- cluster 
    }
    if (mki67 > .7 && (cd8a > .2 || cd8b > .2)) {
      new_id <- cluster
    }
  }
  idents <- c(idents,new_id)
  return(idents)
}


# Build function for mouse data
see_if_keep_cluster_mouse <- function(so, cluster, idents){
  copy_so <- so
  cso <- SubsetData(copy_so, ident.use = cluster, do.clean = TRUE, do.scale = TRUE)
  new_id <- NULL
  cd3e <- sum(GetAssayData(object = cso, slot = "data")["Cd3e",]>0)/nrow(cso@meta.data)
  cd3d <- sum(GetAssayData(object = cso, slot = "data")["Cd3d",]>0)/nrow(cso@meta.data)
  cd3g <- sum(GetAssayData(object = cso, slot = "data")["Cd3g",]>0)/nrow(cso@meta.data)
  cd8b <- sum(GetAssayData(object = cso, slot = "data")["Cd8b1",]>0)/nrow(cso@meta.data)
  cd8a <- sum(GetAssayData(object = cso, slot = "data")["Cd8a",]>0)/nrow(cso@meta.data)
  foxp3 <- sum(GetAssayData(object = cso, slot = "data")["Foxp3",]>0)/nrow(cso@meta.data)
  cd4 <- sum(GetAssayData(object = cso, slot = "data")["Cd4",]>0)/nrow(cso@meta.data)
  if (cd3e > .3 || cd3d > .3 || cd3g > .3) {
    if (cd8b > .3 && cd8a > .3 && foxp3 < .05) {
      new_id <- cluster 
    }
  }
  idents <- c(idents,new_id)
  return(idents)
}


# Emtpy identities
idents <- c()

# Keep certain clusters
if(sample %in% human_samples){
  # Keep clusters for human data
  for(i in clusters){
    idents <- see_if_keep_cluster_human(working_so, as.numeric(i),idents)
  }
} else {
  # Keep clusters for mouse data
  for(i in clusters){
    idents <- see_if_keep_cluster_mouse(working_so, as.numeric(i),idents)
  }
}


# Subset data after keeping certain clusters
working_so1 <- subset(working_so, idents = idents)
working_so1 <- ScaleData(object = working_so1, features = rownames(working_so1))
working_so1 <- FindVariableFeatures(object = working_so1, mean.function = ExpMean, dispersion.function = LogVMR)
hv.genes <- head(x = VariableFeatures(object = working_so1), 1000)
working_so1 <- RunPCA(object = working_so1, pc.genes = hv.genes, do.print = FALSE, pcs.print = 1:5, genes.print = 5, pcs.compute = 50)
working_so1 <- FindNeighbors(object = working_so1, dims = 1:30)
working_so1 <- FindClusters(object = working_so1, resolution = 1.2)
working_so1 <- RunUMAP(object = working_so1, reduction.use = "pca", dims = 1:15, n_neighbors = 15, min_dist = 0.3)

# Plot the UMAP
DimPlot(object = working_so1, reduction = 'umap')

working_so1@meta.data$UMAP_1 <- working_so1@reductions$umap@cell.embeddings[, 1]
working_so1@meta.data$UMAP_2 <- working_so1@reductions$umap@cell.embeddings[, 2]

#Generate UMAPs to check cluster selection 
t1 <- paste0("Clusters Kept=",paste(idents, collapse=', ' ))
p1 <- DimPlot(object = working_so, reduction = "umap", label = TRUE, pt.size = .4) + ggtitle(paste0(length(colnames(data2))," -> ",length(colnames(working_so))," cells (HK/mito filtering)")) + theme(plot.title = element_text(hjust = 0.5)) + coord_fixed() 
p2 <- DimPlot(object = working_so1, reduction = "umap", label = TRUE, pt.size = .4) + ggtitle(paste0(length(colnames(working_so))," -> ",length(colnames(working_so1))," cells (CD8 filtering)")) + theme(plot.title = element_text(hjust = 0.5)) + coord_fixed() + labs(subtitle = t1)
p3 <- FeaturePlot(working_so1,features=c("FOXP3","CD4","CD3G","CD8B","CD8A","CD3D","CD3E"), reduction = "umap",cols=c("gray","red"),pt.size = .4, coord.fixed = TRUE,combine=FALSE)
title <- paste0(sample,"_UMAP_original.pdf")
pdf(title,width=7,height=7)
p1
dev.off()
pdf(paste0(sample,"_UMAP_cd8s.pdf"),width=7,height=7)
p2
dev.off()
pdf(paste0(sample,"_UMAP_genes.pdf"),width=7,height=7)
p3
dev.off()


p4 <- FeaturePlot(working_so,features=c("FOXP3","CD4","CD3G","CD8B","CD8A","CD3D","CD3E"), reduction = "umap",cols=c("gray","red"),pt.size = .4, coord.fixed = TRUE,combine=FALSE)
title1 <- paste0(sample,"_UMAP_genes_noCD8_filtering.pdf")
pdf(title1,width=20,height=6)
p4
dev.off()


# Add TCR metadata from the clone pipeline
# These include all cells, including unfiltered cells
TCR_data <- read.csv('path/to/clone_pipeline_output.csv')

# Remove duplicates (identical rows)
TCR_data <- unique(TCR_data)

# Add rownames to TCR data
rownames(TCR_data) <- TCR_data$Barcode

# Keep only cells in TCR data that are in the filtered Seurat object
idxs <- match(colnames(working_so1), rownames(TCR_data))
TCR_data <- TCR_data[idxs, ]

# Check to see if all the filtered cells were kept
all(rownames(TCR_data) == colnames(working_so1))

# Add the TCR data to the Seurat object
meta <- working_so1@meta.data
meta <- cbind(meta, TCR_data[, c("Matching_pre_filter", "Frequency_pre_filter", "TCR")])
working_so1@meta.data <- meta

# Rename a column
idx = which('seurat_clusters' == colnames(working_so1@meta.data))
colnames(working_so1@meta.data)[idx] <- "Seurat_clusters"


# Save data
mtx_path <- paste0(sample, '_filtered_mtx.csv')
meta_path <- paste0(sample, '_filtered_meta.csv')
rds_path <- paste0(sample, '_filtered.rds')

write.csv(working_so1@meta.data, meta_path, quote = FALSE)
write.csv(working_so1@assays$RNA@counts, mtx_path, quote = FALSE)
saveRDS(working_so1, file = rds_path)
