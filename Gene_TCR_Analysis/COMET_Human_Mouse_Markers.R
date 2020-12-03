
#--------------------------------------#
# Testing if COMET outputs overlap between mouse and human 
#--------------------------------------#

# This script uses the outputs from running COMET on each of the samples
# Both X and L are COMET parameters that are defined in the paper

# Link the path to the COMET outputs
dir_comet <- "COMET_Outputs/"

# List of all the mouse outputs
list_mouse <- c("M1_L_635_X_15", "M2_L_220_X_15", "M3_L_280_X_15",
                "M4_L_847_X_15","M5_L_921_X_15")


# List of all the human outputs
list_human <- c("K409_Blood_L_190_X_15",
                "K411_Blood_L_76_X_15",  "K411_Blood_Longitudinal_L_1494_X_15",
                "K468_Blood_L_782_X_15", "K468_Blood_Longitudinal_L_5477_X_15",
                "K484_Blood_L_482_X_15")

##
## Creating "Universe"
## human to mouse intersect with our default COMET lists to make "background universe". 
##
##

# Read in mapping file to convert between mouse and human gene names
table_h_to_m <- read.table("GeneLists/map_human_mouse.txt", stringsAsFactors = F)
names(table_h_to_m) <- c("human", "mouse")
table_h_to_m$mouse <- toupper(table_h_to_m$mouse)
table_h_to_m$human <- toupper(table_h_to_m$human)



# Read in default gene lists for mouse and human


# Get rows in human_mouse_map that have *mouse* cell surface genes. 
table_mouse_sc<-read.table("GeneLists/Default_genes.txt", stringsAsFactors = F)
table_mouse_sc<-toupper(table_mouse_sc$V1)
length(table_mouse_sc)
#[1] 796
rows_mapping_in_mouse_sc <- which(table_h_to_m$mouse %in% table_mouse_sc)
length(rows_mapping_in_mouse_sc)
# [1] 749

# Get rows in human_mouse_map that have *human* cell surface genes. 
table_human_sc <- read.table("GeneLists/human_Default_genes.txt", stringsAsFactors = F)
table_human_sc <- toupper(table_human_sc$V1)
length(table_human_sc)
# [1] 883

rows_mapping_in_human_sc <- which(table_h_to_m$human %in% table_human_sc)
length(rows_mapping_in_human_sc)
#[1] 883

# Take intersection. These are genes that are in the mapping and are cell surface for both species. 
rows_both_mouse_and_human <- intersect(rows_mapping_in_mouse_sc,rows_mapping_in_human_sc)
length(rows_both_mouse_and_human)
# [1] 749

which(!(rows_mapping_in_mouse_sc %in% rows_both_mouse_and_human))
# integer(0)


universe_table <- table_h_to_m[rows_both_mouse_and_human,]
dim(universe_table)
# [1] 749   2

## Add negations to table
tn <- universe_table
tn$mouse <- paste0(tn$mouse, "_negation")
tn$human <- paste0(tn$human, "_negation")

universe_table <- rbind(universe_table,tn)

## More conservative universe - all genes that are expressed in at least one sample (and their negations.)

##
## Read mouse, make table with all files
##

# Make table with human genes in universe. Filled with one. 
table_m <- data.frame(matrix(1, nrow = nrow(universe_table), ncol = length(list_mouse)))
names(table_m) <- list_mouse
row.names(table_m) <- universe_table$human
list_all_genes_expressed <- c()

## Update the locations with genes in specific file with their q-values. 
mf<-list_mouse[1]
for(mf in list_mouse){
  # This is the file inside of the COMET output
  cf <- paste0(dir_comet, mf, "/data/cluster_1_singleton_full_unranked.csv")
  ct <- read.csv(cf)
  row.names(ct) <- ct$gene_1

  ## Reduce to genes in universe. 
  rows_take_m <- which(row.names(ct) %in% universe_table$mouse)
  ct_in_universe <- ct[rows_take_m,]  
  if(length(which(!(row.names(ct_in_universe) %in% universe_table$mouse)))>0){stop("ERROR !!!")}
  
  ## Mapping to human genes
  cur_mapping_table <- universe_table[which(universe_table$mouse %in% row.names(ct_in_universe)),]
  row.names(cur_mapping_table) <- cur_mapping_table$mouse
  if (!(length(row.names(ct_in_universe)) == length(which(row.names(ct_in_universe) %in% row.names(cur_mapping_table))) )){stop("ERROR!!")}
  ct_in_universe <- merge(ct_in_universe, cur_mapping_table, by=0)
  row.names(ct_in_universe) <- ct_in_universe$human
  ct_in_universe$Row.names <- NULL
  list_all_genes_expressed <- unique(c(list_all_genes_expressed,row.names(ct_in_universe)))
  
  
  ## Update q-values where relevant with values from table. 
  col <- grep(mf, names(table_m))
  table_m[row.names(ct_in_universe),col] <- ct_in_universe$q_value
  
}



##
## Read human
##


# make table with human genes in universe. Filled with one. 
table_h <- data.frame(matrix(1, nrow = nrow(universe_table), ncol = length(list_human)))
names(table_h) <- list_human
row.names(table_h) <- universe_table$human

## update the locations with genes in specific file with their q-values. 
mf <- list_human[1]
for(mf in list_human){
  # This is the file inside of the COMET output
  cf <- paste0(dir_comet, mf, "/data/cluster_1_singleton_full_unranked.csv")
  ct <- read.csv(cf)
  row.names(ct) <- ct$gene_1
  
  ## Reduce to genes in universe. 
  rows_take_h <- which(row.names(ct) %in% universe_table$human)
  ct_in_universe< - ct[rows_take_h, ]  
  if(length(which(!(row.names(ct_in_universe) %in% universe_table$human)))>0){stop("ERROR !!!")}
  
  list_all_genes_expressed <- unique(c(list_all_genes_expressed, row.names(ct_in_universe)))
  
  
  ## update q-values where relevant with values from table. 
  col <- grep(mf, names(table_h))
  table_h[row.names(ct_in_universe), col] <- ct_in_universe$q_value
  
}


##
## Produce lists significant in mouse and in human
##

vec_thresh_pval <- c(0.01, 0.05)
thresh_samples_sig <- c(1, 2, 3)
universe_values <- c(nrow(universe_table), length(list_all_genes_expressed))
universe_values <- nrow(universe_table)

## Make table with results and print out. 

table_results <- data.frame(matrix(0, nrow = length(vec_thresh_pval)*length(thresh_samples_sig)*length(universe_values),ncol=10 ))
names(table_results) <- c("pvalue_threshold", "min_number_of_samples_significant",
                          "pvalue_significance", "number_genes_overlaping",
                          "number_genes_mouse", "number_genes_human",
                          "number_universe", "list_intersecting_markers",
                          "list_mouse_only", "list_human_only")


counter <- 1;

for(n_u in universe_values){
  for(pv in vec_thresh_pval){
    for(ts in thresh_samples_sig){
      num_samples_sig_mouse <- apply(table_m,1,function(x){length(which(x<=pv))})
      list_sig_genes_mouse <- names(num_samples_sig_mouse[which(num_samples_sig_mouse>=ts)])
      length(list_sig_genes_mouse)
      
      num_samples_sig_human <- apply(table_h,1,function(x){length(which(x<=pv))})
      list_sig_genes_human <- names(num_samples_sig_human[which(num_samples_sig_human>=ts)])
      length(list_sig_genes_human)
      
      list_intersect <- intersect(list_sig_genes_mouse,list_sig_genes_human)
      list_m_only <- list_sig_genes_mouse[which(!(list_sig_genes_mouse %in% list_sig_genes_human))]
      list_h_only <- list_sig_genes_human[which(!(list_sig_genes_human %in% list_sig_genes_mouse))]
      
      which(list_m_only %in% list_h_only)
      which(list_h_only %in% list_m_only)
      
      n_i <- length(list_intersect)
      n_r <- length(list_sig_genes_mouse)
      n_b <- length(list_sig_genes_human)
      
      pval <- phyper(n_i-1, n_r, (n_u-n_r), n_b, lower.tail = F)
      
      print(paste0("thresh_pval: ", pv, " thresh_samples_sig: ", ts, " universe: ", n_u))
      print(paste0("n_i: ", n_i, " n_m: ", n_r, " n_h: ", n_b, " n_u: ", n_u))
      print(paste0("pvalue: ", pval))
      print(list_intersect)
      print("\n")
      
      str_int_genes <- paste(list_intersect, collapse = ", ")
      str_m_only_genes <- paste(list_m_only, collapse = ", ")
      str_h_only_genes <- paste(list_h_only, collapse = ", ")
      
      table_results[counter, ] <- c(pv,ts,pval,n_i,n_r,n_b,n_u,str_int_genes,str_m_only_genes,str_h_only_genes)
      counter <- counter + 1
    }
  }
  print("\n\n")
}



write.csv(table_results, file = "TableS11_overlap_mouse_human.csv")






