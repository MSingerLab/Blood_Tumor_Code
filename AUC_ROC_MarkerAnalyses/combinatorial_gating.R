
# Load libraries
library(Seurat)

# Load the function 'Keep_AB_Cells' which subsets the rds object
# to keep only cells with at least one alpha and one beta chain
source("Keep_AB_Cells.R")

# Define some functions
fill_in_structure <- function(gatenum, genes, gate_structures){
  v =  gate_structures[[gatenum]]
  pos = 1
  i = 1
  while (pos <= length(v)){
    if (v[pos]=="") { 
      v[pos] = genes[i] 
      i = i + 1
    }
    pos = pos + 1
  }
  statement = ""
  for (i in 1:length(v)){
    statement = paste(statement, v[i], sep="")
  }
  return(statement)
}

set_up_gates <- function(markers) {
  gate_structures <- list()
  one_gates <- list()
  one_gates[[1]] <- function(v) { v[1] }
  gate_structures[[1]] <- c("")
  
  two_gates <- list()
  two_gates[[1]] <- function(v) { v[1] & v[2] }
  two_gates[[2]] <- function(v) { v[1] | v[2] }
  gate_structures[[2]] <-  c("", " & ", "")
  gate_structures[[3]] <- c("", " | ", "")
  
  three_gates <- list()
  three_gates[[1]] <- function(v) { (v[1] & v[2]) & v[3] }
  three_gates[[2]] <- function(v) { (v[1] | v[2]) & v[3] }
  three_gates[[3]] <- function(v) { (v[1] & v[2]) | v[3] }
  three_gates[[4]] <- function(v) { (v[1] | v[2]) | v[3] }
  gate_structures[[4]] <- c("(", "", " & ", "", ") & ", "")
  gate_structures[[5]] <- c("(", "", " | ", "", ") & ", "")
  gate_structures[[6]] <- c("(", "", " & ", "", ") | ", "")
  gate_structures[[7]] <- c("(", "", " | ", "", ") | ", "")
  
  gates <- c(one_gates, two_gates, three_gates)
  num_genes <- c(1, rep(2, 2), rep(3, 4))
  
  if (length(markers)==4) { 
    four_gates <- list()
    four_gates[[1]] <- function(v) { (v[1] & v[2]) & (v[3] & v[4]) }
    four_gates[[2]] <- function(v) { (v[1] | v[2]) & (v[3] | v[4]) }
    four_gates[[3]] <- function(v) { (v[1] & v[2]) & (v[3] | v[4]) }
    four_gates[[4]] <- function(v) { (v[1] | v[2]) & (v[3] & v[4]) }
    four_gates[[5]] <- function(v) { (v[1] & v[2]) | (v[3] & v[4]) }
    four_gates[[6]] <- function(v) { (v[1] | v[2]) | (v[3] | v[4]) }
    four_gates[[7]] <- function(v) { (v[1] & v[2]) | (v[3] | v[4]) }
    four_gates[[8]] <- function(v) { (v[1] | v[2]) | (v[3] & v[4]) }
    gate_structures[[8]] <- c("(", "", " & ", "", ") & (", "", " & ", "", ")")
    gate_structures[[9]] <- c("(", "", " | ", "", ") & (", "", " | ", "", ")")
    gate_structures[[10]] <- c("(", "", " & ", "", ") & (", "", " | ", "", ")")
    gate_structures[[11]] <- c("(", "", " | ", "", ") & (", "", " & ", "", ")")
    gate_structures[[12]] <- c("(", "", " & ", "", ") | (", "", " & ", "", ")")
    gate_structures[[13]] <- c("(", "", " | ", "", ") | (", "", " | ", "", ")")
    gate_structures[[14]] <- c("(", "", " & ", "", ") | (", "", " | ", "", ")")
    gate_structures[[15]] <- c("(", "", " | ", "", ") | (", "", " & ", "", ")")
    c = length(four_gates) + 1
    shift = sum(length(one_gates), length(two_gates))
    for (i in 1:length(three_gates)){
      th_g = three_gates[[i]]
      four_gates[[c]] <- function(v) { th_g(v[1:3]) & v[4]} 
      gate_structures[[c + length(gates)]] <- c("(", gate_structures[[i + shift]], ") & ", "")
      c = c + 1
      four_gates[[c]] <- function(v) { th_g(v[1:3]) | v[4]} 
      gate_structures[[c + length(gates)]] <- c("(", gate_structures[[i + shift]], ") | ", "")
      c = c + 1
    }
    gates <- c(gates, four_gates)
    num_genes <- c(num_genes, rep(4, length(four_gates)))
  }
  
  gate_names <- vector(mode="numeric")
  one_options <- data.frame(markers)
  num_gates <- length(markers)
  which_gate <- rep(1, length(markers))
  n = 1
  for (i in 1:length(markers)){
    gate_names[n] <- as.character(markers[i])
    n = n + 1
  }
  two_options <- combinations(length(markers), 2, markers)
  num_gates <- num_gates + (nrow(two_options)*length(two_gates))
  which_gate <- c(which_gate, sort(rep(c(2,3), nrow(two_options))))
  
  for (i in 2:3){
    for (j in 1:nrow(two_options)){
      gate_names[n] <- as.character(fill_in_structure(i, two_options[j,], gate_structures))
      n = n + 1
    }
  }
  
  third = vector(mode="character")
  if (length(markers)==4){
    three_options <- rbind(two_options, two_options)
    for (j in 1:nrow(two_options)){
      leftover <- markers[!(markers %in% two_options[j,])]
      third[j] <- leftover[1]
      third[(j+nrow(two_options))] <- leftover[2]
    }
  } else {
    three_options <- two_options
    for (j in 1:nrow(two_options)){
      leftover <- markers[!(markers %in% two_options[j,])]
      third[j] <- leftover
    }
  }
  
  three_options <- data.frame(cbind(three_options, third))
  num_gates <- num_gates + (nrow(three_options)*length(three_gates))
  which_gate <- c(which_gate, sort(rep(c(4,5,6,7), nrow(three_options))))
  three_options <- sapply(three_options, as.character)
  for (i in 4:7){
    for (j in 1:nrow(three_options)){
      gate_names[n] <- as.character(fill_in_structure(i, three_options[j,], gate_structures))
      n <- n + 1
    }
  }
  
  if (length(markers)==4){
    from_two <- two_options
    third <- vector(mode="character")
    fourth <- vector(mode="character")
    for (j in 1:nrow(two_options)){
      leftover <- markers[!(markers %in% two_options[j,])]
      third[j] <- leftover[1]
      fourth[j] <- leftover[2]
    }
    from_two <- cbind(from_two, third, fourth)
    from_three <- three_options
    fourth <- vector(mode="character")
    for (j in 1:nrow(three_options)){
      fourth[j] <- markers[!(markers %in% three_options[j,])]
    }
    from_three <- cbind(from_three, fourth)
    four_options <- data.frame(rbind(from_two, from_three))
    four_options <- sapply(four_options, as.character)
    num_gates <- num_gates + (nrow(four_options)*length(four_gates))
    which_gate <- c(which_gate, sort(rep(seq(8,23,by=1), nrow(four_options))))
    for (i in 8:23){
      for (j in 1:nrow(four_options)){
        gate_names[n] <- as.character(fill_in_structure(i, as.character(four_options[j,]), gate_structures))
        n <- n + 1
      }
    }
    return(list(gates, gate_structures, gate_names, num_gates, which_gate, num_genes, options, one_options, two_options, three_options, four_options))
  } else {
    return(list(gates, gate_structures, gate_names, num_gates, which_gate, num_genes, options, one_options, two_options, three_options))
  }
}

get_sensitivity <- function(did_pass, d) { 
  d$pass = did_pass
  dm = d[d$TM==1,]
  return(nrow(dm[dm$pass==1,])/nrow(dm)) 
}

get_specificity <- function(did_pass, d) { 
  d$pass = did_pass
  dnm = d[d$TM!=1,]
  return(nrow(dnm[dnm$pass==0,])/nrow(dnm)) 
}

get_pass_popsize <- function(did_pass, d) { 
  d$pass = did_pass
  d_pass = d[d$pass==1,]
  return(nrow(d_pass))
}

addPareto <- function(fdf){
  fdf$pareto = "not Pareto-optimal"
  for (i in 1:nrow(fdf)){
    better = fdf[fdf$sensitivity > fdf$sensitivity[i],]
    better = better[better$specificity > fdf$specificity[i],]
    if ((nrow(better)==0) & (fdf$y[i] >= fdf$x[i])) { 
      fdf$pareto[i] = "Pareto-optimal" 
    }
  }
  fdf$pareto = factor(fdf$pareto, levels=c("not Pareto-optimal", "Pareto-optimal"))
  return(fdf)
}

preset_cutoffs <- function(genes, cutoffs, s){
  df = data.frame(matrix(ncol=length(cutoffs), nrow=1))
  for (i in 1:ncol(cutoffs)){
    df[1,i] = cutoffs[s,i]
  }
  colnames(df) = genes
  return(df)
}

compute_combinations <- function(id, markers, sim, neg){
  
  rds = Keep_AB_Cells(paste0('data/', id, '.rds'))
  data = rds@meta.data
  
  # If using human data,
  if (substr(id, 1, 1) == "K"){
    counts = rds@assays$RNA@counts
    for (i in 1:length(markers)){
      ind = which(as.character(rownames(counts))==gsub("_.*$", "", markers[i]))
      data = data.frame(cbind(data, as.numeric(as.character(counts[ind,]))))
      colnames(data)[ncol(data)] = markers[i]
    }
  }
  
  
  res = data.frame(matrix(ncol=(num_gates*2), nrow=nrow(sim)))
  data_match = data[data$TM==1,]
  data_nomatch = data[data$TM!=1,]
  
  for (r in 1:nrow(res)){
    pps = vector(mode="numeric")
    ##each 2 columns in res corresponds to a gate
    j = 1
    while(j < num_gates){
      g = gates[[which_gate[j]]]
      ng = num_genes[which_gate[j]]
      if (ng==1) { 
        options = data.frame(markers)
      }
      if (ng==2) { options = two_options }
      if (ng==3) { options = three_options }
      if (ng ==4) { options = four_options }
      for (k in 1:nrow(options)){
        gene_order = as.vector(options[k,])
        workspace = matrix(nrow=nrow(data), ncol=ng)
        for(q in 1:ng){
          g_code = which(colnames(sim)==gene_order[q])
          col = which(colnames(data)==gene_order[q])
          cutoff = sim[r, g_code]
          if (neg[g_code]){
            workspace[,q] = ifelse(data[,col] <= cutoff, 1, 0)
          } else {
            workspace[,q] = ifelse(data[,col] >= cutoff, 1, 0)
          }
        }
        pass_gate = vector(mode="logical")
        for (q in 1:nrow(workspace)){
          pass_gate[q] = g(workspace[q,])
        }
        res[r,((j*2)-1)] = get_sensitivity(pass_gate, data)
        res[r,(j*2)] = get_specificity(pass_gate, data)
        pps[j] = get_pass_popsize(pass_gate, data)
        j = j + 1
      }
    }
  }
  
  dfr = data.frame(matrix(nrow=(nrow(sim)*(num_gates)), ncol=3))
  colnames(dfr) = c("gate", "sensitivity", "specificity")
  for (i in 1:num_gates){
    e = i*2 - 1
    p = i*2
    dfr[(1+(i-1)*nrow(sim)):(i*nrow(sim)),1] = i 
    dfr[(1+(i-1)*nrow(sim)):(i*nrow(sim)),2:3] = res[,e:p] 
  }
  
  dfr$x = 1 - dfr$specificity
  dfr$y = dfr$sensitivity
  gate_name = gsub("_negation", "", gate_names)
  dfr$gate = gate_name
  df = dfr
  df = addPareto(df)
  df$logic_num = which_gate
  legend = data.frame(cbind(logic_num = seq(1,length(gates),1), num_genes))
  df = left_join(df, legend, by="logic_num")
  df$num_genes = as.factor(df$num_genes)
  corners = data.frame(matrix(nrow=2, ncol=ncol(df)))
  colnames(corners) = colnames(df)
  corners$x = c(0, 1)
  corners$y = c(0, 1)
  corners$pareto = rep("Pareto-optimal", 2)
  corners$num_genes = c(1,1)
  df = data.frame(rbind(df, corners))
  po = df[df$pareto=="Pareto-optimal",]
  df$pareto = as.character(df$pareto)
  
  for (i in 1:nrow(df)){
    if(!is.na(df$gate[i])){
      if (df$gate[i] == "(CCR7 | GYPC) | (LTB | FLT3LG)") { df$pareto[i] = "sens" }
      if (df$gate[i] == "((CCR7 & GYPC) & FLT3LG) | LTB") { df$pareto[i] = "balance" }
      if (df$gate[i] == "(CCR7 & LTB) & (GYPC & FLT3LG)") { df$pareto[i] = "spec" }
    }
  }
  
  df$num_genes = as.numeric(as.character(df$num_genes))
  df = df[order(df$num_genes),]
  df$num_genes = as.factor(df$num_genes)
  df$pareto = factor(df$pareto, levels=c("not Pareto-optimal", "Pareto-optimal", "balance", "spec", "sens"))
  df = df[order(df$pareto),]
  
  pdf(paste0('outputs/', paste(id, paste(markers, collapse='_'), sep="_"), "_pareto_plot.pdf"), height=6, width=6)
  g = ggplot()
  g = g + geom_line(aes(po$x, po$y), linetype="dashed", color="deepskyblue3")
  g = g + geom_point(aes(df$x, df$y, color=df$pareto, shape=df$num_genes, size=df$pareto), show.legend = FALSE) + theme_bw()  + scale_color_manual(values = c("dimgray", "deepskyblue3", "darkslateblue", "darkmagenta", "darkgreen")) + scale_shape_manual(values=c(16,15,17,18)) + scale_size_manual(values=c(2, 2, 6, 6, 6))
  g = g + xlab("1 - Specificity") + ylab("Sensitivity") + scale_x_continuous(limits=c(0,1), breaks = seq(0, 1, .2)) + scale_y_continuous(limits = c(0,1), breaks = seq(0, 1, .2)) ##+ xlim(c(0,1)) + ylim(c(0,1))
  g = g + geom_abline(slope=1, intercept = 0, size=0.25) 
  print(g)
  dev.off()
  
  write.csv(df, paste0('outputs/', id, "_combinatorial_sens_spec.csv"))
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

library(gtools)
library(DescTools)
library(ggplot2)
library(ggrepel)

markers = c("NKG2D_TotalSeqC", "CX3CR1_TotalSeqC", "CD39_TotalSeqC")
neg = rep(FALSE, 3)
gate_info = set_up_gates(markers)
gates <- gate_info[[1]]
gate_structures <- gate_info[[2]]
gate_names <- gate_info[[3]]
num_gates <- gate_info[[4]]
which_gate <- gate_info[[5]]
num_genes <- gate_info[[6]]
options <- gate_info[[7]]
one_options <- gate_info[[8]]
two_options <- gate_info[[9]]
three_options <- gate_info[[10]]
if (length(markers)==4){
  four_options <- gate_info[[11]]
}

samples = c("M4_Blood", "M5_Blood")
cutoffs = set_up_cutoffs(markers, samples)

for (s in 1:length(samples)){
  id = samples[s]
  sim = preset_cutoffs(markers, cutoffs, s)
  colnames(sim) = markers
  compute_combinations(id, markers, sim, neg)
}

rm(list = setdiff(ls(), lsf.str()))

markers = c("LTB_negation", "CCR7_negation", "GYPC_negation", "FLT3LG_negation")
neg = c(TRUE, TRUE, TRUE, TRUE)
gate_info = set_up_gates(markers)
gates <- gate_info[[1]]
gate_structures <- gate_info[[2]]
gate_names <- gate_info[[3]]
num_gates <- gate_info[[4]]
which_gate <- gate_info[[5]]
num_genes <- gate_info[[6]]
options <- gate_info[[7]]
one_options <- gate_info[[8]]
two_options <- gate_info[[9]]
three_options <- gate_info[[10]]
if (length(markers)==4){
  four_options <- gate_info[[11]]
}

samples = c('K409_Blood', 'K411_Blood', 'K411_Blood_Longitudinal', 'K468_Blood', 'K468_Blood_Longitudinal', 'K484_Blood')
cutoffs = set_up_cutoffs(markers, samples)

for (s in 1:length(samples)){
  id = samples[s]
  sim = preset_cutoffs(markers, cutoffs, s)
  colnames(sim) = markers
  compute_combinations(id, markers, sim, neg)
}
