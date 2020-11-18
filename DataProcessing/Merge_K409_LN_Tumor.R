#####################################################################
## Merge K409 lymph node and tumor CD8 data into one Seurat object ##
#####################################################################

## load data
K409_LN <- readRDS("Data/K409_LN.rds")
K409_LN$tissue <- "LN"
K409_LN$sample <- "K409_LN"
K409_Tumor <- readRDS("K409_Tumor.rds")
K409_Tumor$tissue <- "Tumor"
K409_Tumor$sample <- "K409_Tumor"


## combine
K409_merged <- merge(K409_LN, K409_Tumor)

## SCTransform
K409_merged <- SCTransform(K409_merged, return.only.var.genes = FALSE)

## reduction
cat("Running reductions...\n")
K409_merged <- RunPCA(K409_merged)
K409_merged <- RunUMAP(K409_merged, reduction = "pca", dims = 1:30)

## clustering
cat("Clustering...\n")
K409_merged <- FindNeighbors(K409_merged)
K409_merged <- FindClusters(K409_merged, resolution = 0.3)

## save object
cat("Saving object...\n")
saveRDS(K409_merged, file = paste0(dir_out, "K409_merged.rds"))
