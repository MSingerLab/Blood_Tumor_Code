#########################
## Integrate data sets ##
#########################
##  Use Seurat v3 integration (SCTransform work flow)


#### configuration ####
library(Seurat)
library(ggplot2)
library(dplyr)
options(future.globals.maxSize = +Inf)

dir_data <- "../data/"
dir_out <- "../outputs/"
if (!dir.exists(dir_out)) {
  dir.create(dir_out, recursive = T)
  dir.create(paste0(dir_out, "Plots/"), recursive = F)
  dir.create(paste0(dir_out, "DEgenes/"), recursive = F)
}


#### integrate data #####
sink(paste0(dir_out, "R_terminal_output.txt"))

## load data
file_names <- list.files(dir_data, recursive = T, full.names = F)

# Keep human data, files that start with 'K'
file_names <- file_names[which(substr(file_names, 1, 1) == 'K')]

# Remove '.rds' at end of file names
sample_names <- sub(".rds", "", gsub(".*/", "", file_names))
so_list_all <- lapply(seq_along(file_names), function(i) {
  tmp <- readRDS(paste0(dir_data, file_names[i]))
  ## add meta data
  tmp$sample <- sample_names[i]
  tissue <- substr(sample_names[1], 6, 25)
  tmp
})
names(so_list_all) <- sample_names

## SCTransform
cat("\n\n\n===================  SCTransform  ===================\n")
so_list_all <- lapply(so_list_all, FUN = SCTransform, return.only.var.genes = FALSE)
sample_names <- names(so_list_all)


#### integrate data and clustering
## different data combinations
## blood_no_longitudinal: initial blood samples only; K409 blood, K411 blood, K468 blood, K484 blood
## tumor_ln_primary: all tumor samples, both LN met and primary; K409 LN met, K409 primary tumor; k311 LN met; K411 LN met; K468 axillary subcu. mass, K484 LN met
to_integrate_vec <- c("Blood", "Blood_wLongitudinal") 

for (to_integrate in to_integrate_vec) {
  if (to_integrate == "Blood") {
    so_list <- so_list_all[setdiff(grep("Blood", sample_names), grep("Longitudinal", sample_names))]
  } else if (to_integrate == "Blood_wLongitudinal") {
    so_list <- so_list_all 
  }
  cat("\n\n\n===================  Integrating ", to_integrate, "  ===================\n")
  print(names(so_list))
  
  ## Find anchors
  cat("Finding anchors...\n")
  features <- SelectIntegrationFeatures(object.list = so_list, nfeatures = Inf)
  so_list <- PrepSCTIntegration(object.list = so_list, 
                                anchor.features = features, 
                                verbose = TRUE)
  so_list <- lapply(X = so_list, FUN = RunPCA, verbose = TRUE, features = features)
  anchors <- FindIntegrationAnchors(object.list = so_list, 
                                    normalization.method = "SCT", 
                                    reduction = "rpca",
                                    anchor.features = features, 
                                    verbose = TRUE)
  
  ## Integrate
  cat("Integrating data...\n")
  so <- IntegrateData(anchorset = anchors, 
                      normalization.method = "SCT", 
                      verbose = TRUE)
  rm(so_list, anchors)
  DefaultAssay(so) <- "integrated"
  
  ## reduction
  cat("Running reductions...\n")
  so <- RunPCA(so)
  so <- RunUMAP(so, reduction = "pca", dims = 1:30)
  
  ## clustering
  cat("Clustering...\n")
  so <- FindNeighbors(so)
  so <- FindClusters(so, resolution = 0.3)
  
  ## save object
  cat("Saving object...\n")
  saveRDS(so, file = paste0(dir_out, to_integrate, "_integrated.rds"))
  
  ## visualization
  pdf(paste0(dir_out, "PLOTS/", to_integrate, "_cluster_overiew.pdf"))
  p <- so@meta.data %>% group_by(sample) %>% tally() %>%
    ggplot(aes(x = sample, y = n, fill = sample)) +
    geom_bar(stat = "identity") + geom_text(aes(label = n)) +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Sample Sizes (all samples)")
  print(p)
  p <- so@meta.data %>% group_by(seurat_clusters) %>% tally() %>%
    ggplot(aes(x = seurat_clusters, y = n, fill = seurat_clusters)) +
    geom_bar(stat = "identity") + geom_text(aes(label = n)) +
    theme_bw() + ggtitle("Cluster Sizes")
  print(p)
  p <- DimPlot(so, group.by = "sample") + ggtitle("color by samples")
  print(p)
  p <- DimPlot(so, group.by = "seurat_clusters") + ggtitle("color by clusters")
  print(p)
  dev.off()
  
  ## cluster differential expression
  cat("DE analyses...\n")
  DE_ranksum <- FindAllMarkers(so, only.pos = FALSE, max.cells.per.ident = 1000) ## downsample to speed up computation
  write.csv(DE_ranksum, file = paste0(dir_out, "DE/", to_integrate, "_cluster_deg_ranksum.csv"))
}

sink()
