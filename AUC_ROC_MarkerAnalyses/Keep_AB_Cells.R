Keep_AB_Cells <- function(rds_file){
  # A function that subsets the rds objects to keep only AB cells,
  # which are cells that have at least one alpha and one beta chain
  # in the TCR. This function works by removing any cells that do not
  # have 'alpha_chain_only', 'beta_chain_only', or 'notcr' in their
  # 'Matching' column/status.
  
  # Load Seurat library
  library(Seurat)
  
  # Define the 'not in' operator
  `%not_in%` = Negate(`%in%`)
  
  # Load the rds object
  rds <- readRDS(rds_file)
  
  # Keep cells that do not have a Matching status as 
  # 'alpha_chain_only', 'beta_chain_only', or 'notcr'
  rds <- subset(rds, subset = Matching %not_in% c('alpha_chain_only', 'beta_chain_only', 'notcr'))
  
  # For downstream analyses, label non-matching cells with TM = 2,
  # and label matching cells with TM = 1
  rds$TM = 2
  rds$TM[! rds$Matching=="not_matching"] = 1
  
  # Save as factor
  rds$TM <- factor(rds$TM)
  
  return(rds)
  
} 
