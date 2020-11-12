
# Load libraries
library(Seurat)
library(DescTools)

# Load the function 'Keep_AB_Cells' which subsets the rds object
# to keep only cells with at least one alpha and one beta chain
source("Keep_AB_Cells.R")


# Create function that computes AUC for each marker
# Outputs csv files

compute_aucs <- function(sample, resol){
  # sample = sample name (ex: K409_Blood, K411_Blood, etc.)
  # resol = resolution
  
  # Load rds object, keep only AB cells
  rds <- Keep_AB_Cells(paste0(sample, '.rds'))
  
  # marker_data = count matrix (genes as rows, cells as columns)
  # md = metadata
  marker_data = rds@assays$RNA@counts
  md = rds@meta.data
  
  # Get all marker names
  orig_marker = rownames(marker_data)
  
  # Create new vector to include original and negative markers
  marker = vector(mode="character", length=(length(orig_marker)*2))
  
  for (i in 1:length(marker)){
    if (i%%2 !=0){
      marker[i] = orig_marker[(i+1)/2]
    }else{
      marker[i] = paste(orig_marker[(i/2)], "negation", sep="_")
    }
  }
  
  # Allocate memory to calculate auc for each marker
  auc = vector(mode="numeric", length=nrow(marker_data*2))
  
  # md = metadata
  # md$Barcode = cell barcodes
  md$Barcode = factor(md$Barcode, levels=as.character(colnames(marker_data)))
  md = md[order(md$Barcode),]
  md$marker = 0
  
  # For each marker,
  for (i in 1:(length(marker)/2)){
    
    # Save count data into 'data' variable (all cells, one gene)
    data = as.numeric(as.character(marker_data[i,]))
    
    # Create vector that binds a sequence that runs over the data
    # with a quantile vector
    dist = c(seq(min(data), max(data), length=(1/resol)), quantile(data, probs=seq(resol,(1-resol),by=resol)))
    
    # Initialize four vectors
    posx = vector(mode="numeric", length=length(dist))
    posy = vector(mode="numeric", length=length(dist))
    negx = vector(mode="numeric", length=length(dist))
    negy = vector(mode="numeric", length=length(dist))
    
    # Insert count data into metadata object
    md$marker = as.numeric(as.character(marker_data[i,]))
    
    # Tumor matching cells have TM == 1 and non-tumor matching cells have TM == 2
    md_match = md[md$TM==1,]
    md_nomatch = md[md$TM!=1,]
    
    for (j in 1:length(dist)){
      posy[j] = nrow(md_match[md_match$marker>=dist[j],])/nrow(md_match)
      posx[j] = 1 - (nrow(md_nomatch[md_nomatch$marker<dist[j],])/nrow(md_nomatch))
      negy[j] = nrow(md_match[md_match$marker<=dist[j],])/nrow(md_match)
      negx[j] = 1 - (nrow(md_nomatch[md_nomatch$marker>dist[j],])/nrow(md_nomatch))
    }
    posx = c(posx, 0, 1)
    posy = c(posy, 0, 1)
    negx = c(negx, 0, 1)
    negy = c(negy, 0, 1)
    
    # Use AUC function from DescTools to calculate the area under curve
    auc[(2*i)-1] = AUC(posx, posy)
    auc[2*i] = AUC(negx, negy)
                 
  }
  
  # Output is a dataframe with two columns:
  # First column is the marker name (including negation markers)
  # Second column is the AUC for that marker
  df = data.frame(cbind(marker, auc))
  return(df)
}



# Write short script to loop over each input
samples = c("K409_Blood", "K411_Blood", "K411_Blood_Longitudinal",
            "K468_Blood", "K468_Blood_Longitudinal", "K484_Blood")

for (i in 1:length(samples)){
  print(samples[i])
  input = paste0('data/', samples[i])
	df = compute_aucs(input, 0.05)
	df$auc = as.numeric(as.character(df$auc))
	df = df[order(-df$auc),]
	write.csv(df, paste0('outputs/', samples[i], "_aucs.csv"))
}

