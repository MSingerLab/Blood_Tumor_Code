library(Seurat)

# Define samples
samples <- c("M1_Blood", "M2_Blood", "M3_Blood",
             "M1_Tumor", "M2_Tumor", "M3_Tumor")

# Use this loop to load rds objects
for(i in 1:length(samples)){
  sample <- samples[i]
  
  # Load file name
  file_name <- paste0('path/to/data/', sample, '.rds')
  
  # Load rds object
  rds <- readRDS(file_name)
  
  # Assign rds object to sample
  assign(sample, rds)
}


# In each sample, replace 'matching' with 'sample_matching'
# For example, in the sample M1_Blood, the level 'matching' in the
# 'Matching' column will be replaced with 'M1_Blood_matching'
# This is done for convenience, since all the data will be integrated

seurat_objects <- c(M1_Blood, M2_Blood, M3_Blood,
                    M1_Tumor, M2_Tumor, M3_Tumor)

library(plyr)
for(i in 1:length(samples)){
  
  sample <- samples[i]
  rds <- seurat_objects[[i]]
  
  # Rename data
  rds$Matching <- revalue(rds$Matching, c("matching" = paste0(sample, "_matching")))
  
  # Assign rds object to sample
  assign(sample, rds)
  
}



# Integrate blood data
samples.list <- c(M1_Blood, M2_Blood, M3_Blood)

# First run SCTransform
for(i in 1:length(samples.list)){
  samples.list[[i]] <- SCTransform(samples.list[[i]], verbose = T)
}

# To avoid errors,
Sys.setenv('R_MAX_VSIZE'= 100000000000)
options(future.globals.maxSize = 3145728000)
# Note: 3000 * 1024^2 bytes to resolve error with 2GB list

# Keep top 3000 genes for integration
samples.features <- SelectIntegrationFeatures(object.list = samples.list, nfeatures = 3000)
samples.list <- PrepSCTIntegration(object.list = samples.list,
                                   anchor.features = samples.features,
                                   verbose = T)

# Find the integration anchors
samples.anchors <- FindIntegrationAnchors(object.list = samples.list,
                                          normalization.method = "SCT", 
                                          anchor.features = samples.features,
                                          verbose = T)

# Integrate the data
samples.integrated <- IntegrateData(anchorset = samples.anchors,
                                    normalization.method = "SCT", 
                                    verbose = T)

# Change the default assay
DefaultAssay(object = samples.integrated) <- "RNA"

# Scale data, find variable genes, run PCA
samples.integrated <- ScaleData(samples.integrated  , verbose = T)
#samples.integrated <- FindVariableFeatures(object = samples.integrated, mean.function = ExpMean,                                                 dispersion.function = LogVMR)
#samples.integrated <- FindVariableFeatures(samples.integrated)
samples.integrated <- RunPCA(samples.integrated, npcs = 30, verbose = T)

# Find nearest neighbors, run clustering, run UMAP
# Change resolution to get about 7 clusters
samples.integrated <- FindNeighbors(object = samples.integrated , dims = 1:16)
samples.integrated <- FindClusters(samples.integrated , resolution = 0.15, verbose=T)
samples.integrated <- RunUMAP(samples.integrated , reduction = "pca", dims = 1:16)

# Plot the UMAP
DimPlot(samples.integrated, reduction = "umap")

# Save RDS
saveRDS(samples.integrated, 'Path/to/Data/MouseIntegratedBlood.rds')







# Integrate tumor data
samples.list <- c(M1_Tumor, M2_Tumor, M3_Tumor)

# First run SCTransform
for(i in 1:length(samples.list)){
  samples.list[[i]] <- SCTransform(samples.list[[i]], verbose = T)
}

# To avoid errors,
Sys.setenv('R_MAX_VSIZE'= 100000000000)
options(future.globals.maxSize = 3145728000)
# Note: 3000 * 1024^2 bytes to resolve error with 2GB list

# Keep top 3000 genes for integration
samples.features <- SelectIntegrationFeatures(object.list = samples.list, nfeatures = 3000)
samples.list <- PrepSCTIntegration(object.list = samples.list,
                                   anchor.features = samples.features,
                                   verbose = T)

# Find the integration anchors
samples.anchors <- FindIntegrationAnchors(object.list = samples.list,
                                          normalization.method = "SCT", 
                                          anchor.features = samples.features,
                                          verbose = T)

# Integrate the data
samples.integrated <- IntegrateData(anchorset = samples.anchors,
                                    normalization.method = "SCT", 
                                    verbose = T)

# Change the default assay
DefaultAssay(object = samples.integrated) <- "RNA"

# Scale data, find variable genes, run PCA
samples.integrated <- ScaleData(samples.integrated  , verbose = T)
#samples.integrated <- FindVariableFeatures(object = samples.integrated, mean.function = ExpMean,                                                 dispersion.function = LogVMR)
#samples.integrated <- FindVariableFeatures(samples.integrated)
samples.integrated <- RunPCA(samples.integrated, npcs = 30, verbose = T)

# Find nearest neighbors, run clustering, run UMAP
# Change resolution to get about 7 clusters
samples.integrated <- FindNeighbors(object = samples.integrated , dims = 1:16)
samples.integrated <- FindClusters(samples.integrated , resolution = 0.15, verbose=T)
samples.integrated <- RunUMAP(samples.integrated , reduction = "pca", dims = 1:16)

# Plot the UMAP
DimPlot(samples.integrated, reduction = "umap")

# Save RDS
saveRDS(samples.integrated, 'Path/to/Data/MouseIntegratedTumor.rds')

