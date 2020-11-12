# This file defines matching status and counts clone size between
# blood and tumor samples using the filtered rds objects.
# The additional restriction is that matching status can only be
# applied to AB cells (or 'alpha/beta cells'), which are cells
# that have at least one alpha AND one beta chain in the TCR.
# All the TCRs are written in the form
#
#  alpha_chain|beta_chain
#
# For example, the TCR sequence CALSDRGGSNAKLTF|CASTNSAETLYF belongs
# to an AB cell, but the TCR sequence |CASSDGGAYAEQFF is missing an
# alpha chain, and would therefore not be considered an AB cell.
# Additionally, the clone size of all non-AB cells are defined to be 0.

# Load libraries
library(hash)
library(dplyr)
library(plyr)


############################################################
#                                                          #
#    Create function to enforce AB matching criteria       #
#                                                          #
############################################################


# Create function to enforce AB matching criteria
AB_Matching <- function(blood_rds, tumor_rds){
  
  # blood = blood Seurat/rds object
  # tumor = tumor Seurat/rds object
  
  # Get just the metadata
  blood_df <- blood_rds@meta.data
  tumor_df <- tumor_rds@meta.data
  
  # Turn TCR into character values instead of factors
  blood_df$TCR <- as.character(blood_df$TCR)
  tumor_df$TCR <- as.character(tumor_df$TCR)
  
  # Replace NA with 'notcr' in the 'TCR' column
  blood_df$TCR[ is.na(blood_df$TCR) ] <- 'notcr'
  tumor_df$TCR[ is.na(tumor_df$TCR) ] <- 'notcr'
  
  # Get a vector of all the unique TCRs for blood and tumor
  Unique_blood_TCR <- unique(blood_df$TCR)
  Unique_tumor_TCR <- unique(tumor_df$TCR)
  
  # From the vector, remove 'notcr'
  Unique_blood_TCR <- Unique_blood_TCR[ Unique_blood_TCR != 'notcr' ]
  Unique_tumor_TCR <- Unique_tumor_TCR[ Unique_tumor_TCR != 'notcr' ]

  
  # Count clone sizes
  # Non-AB cells will be dealt with later
  blood_sizes <- dplyr::count(blood_df, TCR, sort = TRUE)
  tumor_sizes <- dplyr::count(tumor_df, TCR, sort = TRUE)
  
  # Now get a vector of all the TCRs in blood and tumor
  all_blood_TCRs <- blood_df$TCR
  all_tumor_TCRs <- tumor_df$TCR
  
  # Run matching analysis on blood
  # Initialize empty vectors
  matching <- character(length = nrow(blood_df))
  clone_size <- rep(0, nrow(blood_df))
  
  for(i in 1:length(all_blood_TCRs)){
    # For each TCR
    TCR <- all_blood_TCRs[i]
    
    # If it has no sequenced TCR, call it 'notcr' and give it a clone size of 0
    if(TCR == 'notcr'){
      matching[i] <- 'notcr'
      clone_size[i] <- 0
    }
    
    # If it's missing an alpha chain, call it 'beta_chain_only', give it a
    # clone size of 0
    else if(substr(TCR, 1, 1) == "|"){                   # Check first character of TCR
      matching[i] <- "beta_chain_only"
      clone_size[i] <- 0
    }
    
    # If it's missing a beta chain, call it 'alpha_chain_only', give it a
    # clone size of 0
    else if(substr(TCR, nchar(TCR), nchar(TCR)) == "|"){ # Check last character of TCR
      matching[i] <- "alpha_chain_only"
      clone_size[i] <- 0
    }
    
    # If the TCR is in Unique_tumor_TCR, call it 'matching'
    else if(TCR %in% Unique_tumor_TCR){
      matching[i] <- "matching"
      clone_size[i] <- blood_sizes[blood_sizes$TCR == TCR, 'n']
    }
    
    # If the TCR is not in Unique_tumor_TCR, call it non-matching'
    else {
      matching[i] <- "not_matching"
      clone_size[i] <- blood_sizes[blood_sizes$TCR == TCR, 'n']
    }
    
  }
  
  # Add the results back to the blood dataframe
  blood_df$Matching <- matching
  blood_df$Clone.size <- clone_size
  
  
  
  # Run matching analysis on tumor
  # Initialize empty vectors
  matching <- character(length = nrow(tumor_df))
  clone_size <- rep(0, nrow(tumor_df))
  
  for(i in 1:length(all_tumor_TCRs)){
    # For each TCR
    TCR <- all_tumor_TCRs[i]
    
    # If it has no sequenced TCR, call it 'notcr' and give it a clone size of 0
    if(TCR == 'notcr'){
      matching[i] <- 'notcr'
      clone_size[i] <- 0
    }
    
    # If it's missing an alpha chain, call it 'beta_chain_only', give it a
    # clone size of 0
    else if(substr(TCR, 1, 1) == "|"){                   # Check first character of TCR
      matching[i] <- "beta_chain_only"
      clone_size[i] <- 0
    }
    
    # If it's missing a beta chain, call it 'alpha_chain_only', give it a
    # clone size of 0
    else if(substr(TCR, nchar(TCR), nchar(TCR)) == "|"){ # Check last character of TCR
      matching[i] <- "alpha_chain_only"
      clone_size[i] <- 0
    }
    
    # If the TCR is in Unique_blood_TCR, call it 'matching'
    else if(TCR %in% Unique_blood_TCR){
      matching[i] <- "matching"
      clone_size[i] <- tumor_sizes[tumor_sizes$TCR == TCR, 'n']
    }
    
    # If the TCR is not in Unique_blood_TCR, call it non-matching'
    else {
      matching[i] <- "not_matching"
      clone_size[i] <- tumor_sizes[tumor_sizes$TCR == TCR, 'n']
    }
    
  }
  
  # Add the results back to the tumor dataframe
  tumor_df$Matching <- matching
  tumor_df$Clone.size <- clone_size
  
  
  
  # Now that the matching is done, add the dataframes back to the Seurat objects
  blood_rds@meta.data <- blood_df
  tumor_rds@meta.data <- tumor_df
  
  return(list(blood_rds, tumor_rds))
  
}


#########################################################
#                                                       #
#     Run matching on filtered mouse rds objects        #
#                                                       #
#########################################################

# Run matching on filtered mouse rds objects
blood_samples <- c("M1_Blood", "M2_Blood", "M3_Blood", "M4_Blood", "M5_Blood")
tumor_samples <- c("M1_Tumor", "M2_Tumor", "M3_Tumor", "M4_Tumor", "M5_Tumor")

# Loop over each paired sample, run matching
for(i in 1:length(blood_samples)){
  
  # Load Seurat objects
  blood_rds <- readRDS(paste0('RDSobjects/', blood_samples[i], '.rds'))
  tumor_rds <- readRDS(paste0('RDSobjects/', tumor_samples[i], '.rds'))
  
  # Run matching
  output <- AB_Matching(blood_rds, tumor_rds)
  
  blood_rds <- output[[1]]
  tumor_rds <- output[[2]]
  
  # Save updated rds objects
  saveRDS(blood_rds, paste0('RDSobjects/', blood_samples[i], '.rds'))
  saveRDS(tumor_rds, paste0('RDSobjects/', tumor_samples[i], '.rds'))
  
}


#########################################################
#                                                       #
#     Run matching on K409 filtered rds objects         #
#                                                       #
#########################################################


# Matching overview for K409:
# Cells in tumor are matching if they match to blood
# Cells in lymph node (LN) are matching if they match to blood
# Cells in blood are matching if they match to either tumor or LN

# Method:
# Run matching between tumor and blood
# Run matching between LN and blood
# Merge tumor and LN, run the matching between
# the blood and merged tumor/LN data

K409_Blood <- readRDS('RDSobjects/K409_Blood.rds')
K409_Tumor <- readRDS('RDSobjects/K409_Tumor.rds')
K409_LN <- readRDS('RDSobjects/K409_LN.rds')

# Run matching between blood and tumor
output <- AB_Matching(K409_Blood, K409_Tumor)
K409_Blood <- output[[1]]
K409_Tumor <- output[[2]]

# Run matching between blood and LN
output <- AB_Matching(K409_Blood, K409_LN)
K409_Blood <- output[[1]]
K409_LN <- output[[2]]

# Merge tumor samples, run matching
tumor_ln <- merge(K409_Tumor, K409_LN)
output <- AB_Matching(K409_Blood, tumor_ln)

# Keep only blood data, since matching for tumor rds objects are done
K409_Blood <- output[[1]]

# Save objects
saveRDS(K409_Blood, 'RDSobjects/K409_Blood.rds')
saveRDS(K409_Tumor, 'RDSobjects/K409_Tumor.rds')
saveRDS(K409_LN, 'RDSobjects/K409_LN.rds')


#########################################################
#                                                       #
#     Run matching on K411 filtered rds objects         #
#                                                       #
#########################################################

# Matching overview for K411:
# Cells in blood are matching if they match to tumor
# Cells in blood longitudinal are matching if they match to tumor
# Cells in tumor are matching if they match to either blood sample

# Method:
# Run matching between blood and tumor
# Run matching between blood longitudinal and tumor
# Merge blood and blood longitudinal, run the matching between
# the tumor and merged blood/blood longitudinal data

K411_Blood <- readRDS('RDSobjects/K411_Blood.rds')
K411_Blood_Long <- readRDS('RDSobjects/K411_Blood_Longitudinal.rds')
K411_Tumor <- readRDS('RDSobjects/K411_Tumor.rds')

# Run matching between blood and tumor
output <- AB_Matching(K411_Blood, K411_Tumor)
K411_Blood <- output[[1]]
K411_Tumor <- output[[2]]

# Run matching between blood longitudinal and tumor
output <- AB_Matching(K411_Blood_Long, K411_Tumor)
K411_Blood_Long <- output[[1]]
K411_Tumor <- output[[2]]

# Merge blood samples, run matching
bloods <- merge(K411_Blood, K411_Blood_Long)
output <- AB_Matching(bloods, K411_Tumor)

# Keep only tumor data, since matching for blood rds objects are done
K411_Tumor <- output[[2]]

# Save RDS objects
saveRDS(K411_Blood, 'RDSobjects/K411_Blood.rds')
saveRDS(K411_Blood_Long, 'RDSobjects/K411_Blood_Longitudinal.rds')
saveRDS(K411_Tumor, 'RDSobjects/K411_Tumor.rds')


#########################################################
#                                                       #
#     Run matching on K468 filtered rds objects         #
#                                                       #
#########################################################

# Matching overview for K468:
# Cells in blood are matching if they match to tumor
# Cells in blood longitudinal are matching if they match to tumor
# Cells in tumor are matching if they match to either blood sample

# Method:
# Run matching between blood and tumor
# Run matching between blood longitudinal and tumor
# Merge blood and blood longitudinal, run the matching between
# the tumor and merged blood/blood longitudinal data

K468_Blood <- readRDS('RDSobjects/K468_Blood.rds')
K468_Blood_Long <- readRDS('RDSobjects/K468_Blood_Longitudinal.rds')
K468_Tumor <- readRDS('RDSobjects/K468_Tumor.rds')

# Run matching between blood and tumor
output <- AB_Matching(K468_Blood, K468_Tumor)
K468_Blood <- output[[1]]
K468_Tumor <- output[[2]]

# Run matching between blood longitudinal and tumor
output <- AB_Matching(K468_Blood_Long, K468_Tumor)
K468_Blood_Long <- output[[1]]
K468_Tumor <- output[[2]]

# Merge blood samples, run matching
bloods <- merge(K468_Blood, K468_Blood_Long)
output <- AB_Matching(bloods, K468_Tumor)

# Keep only tumor data, since matching for blood rds objects are done
K468_Tumor <- output[[2]]

# Save RDS objects
saveRDS(K468_Blood, 'RDSobjects/K468_Blood.rds')
saveRDS(K468_Blood_Long, 'RDSobjects/K468_Blood_Longitudinal.rds')
saveRDS(K468_Tumor, 'RDSobjects/K468_Tumor.rds')


#########################################################
#                                                       #
#     Run matching on K484 filtered rds objects         #
#                                                       #
#########################################################

# Matching overview for K468:
# Cells in blood are matching if they match to tumor
# Cells in tumor are matching if they match to blood

# Method:
# Run matching between blood and tumor

K484_Blood <- readRDS('RDSobjects/K484_Blood.rds')
K484_Tumor <- readRDS('RDSobjects/K484_Tumor.rds')

# Run matching
output <- AB_Matching(K484_Blood, K484_Tumor)
K484_Blood <- output[[1]]
K484_Tumor <- output[[2]]

# Save rds objects
saveRDS(K484_Blood, 'RDSobjects/K484_Blood.rds')
saveRDS(K484_Tumor, 'RDSobjects/K484_Tumor.rds')





