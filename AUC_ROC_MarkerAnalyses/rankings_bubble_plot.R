source("kaitlyncolors.R")

# Define some functions
processCometFile <- function(id){
  d = read.csv(paste(id, "singletons_HS_full.csv", sep="_"))
  d$sample = id
  return(d)
}

processAUCFile <- function(id){
  d = read.csv(paste(id, "aucs.csv", sep="_")) ##
  d$sample = id
  return(d)
}

getSigQs <- function(id){
  d = read.csv(paste(id, "_singletons_HS_full.csv", sep=""))
  d = d[d$q_value<0.05,]
  return(as.character(d$gene_1))
}

# Sample names
samples = c("K411_Blood","K411_Blood_Longitudinal", "K468_Blood",
            "K468_Blood_Longitudinal", "K409_Blood",  "K484_Blood")

# Load all the data by changing directories
setwd('data/')
all_comet_results <- processCometFile(samples[1])
all_auc_results <- processAUCFile(samples[1])
all_auc_results = all_auc_results[order(-all_auc_results$auc),]
marker = as.character(all_auc_results$marker)
sig_qs = getSigQs(samples[1])

for (i in 2:length(samples)){
  all_comet_results = data.frame(rbind(all_comet_results, processCometFile(samples[i])))
  sig_qs = c(sig_qs, getSigQs(samples[i]))
  a = processAUCFile(samples[i])
  all_auc_results = data.frame(rbind(all_auc_results, a))
}

# Return to main directory
setwd('..')

# Start working on analyses
t = data.frame(table(sig_qs))
t = t[order(-t$Freq),]
t = data.frame(cbind(1:length(unique(sig_qs)), t)) 
colnames(t)[ncol(t)-1] = "marker"
near_consensus = as.character(t$marker[t$Freq >=4])
exclude = c("CD8A", "CCL4", "CCL5", "MIF_negation")
near_consensus = near_consensus[!(near_consensus %in% exclude)]

consensus_results = all_auc_results[all_auc_results$marker %in% near_consensus,]
ranking_df = data.frame(cbind(near_consensus))
colnames(ranking_df)[1] = "marker"
ranking_df$marker = as.character(ranking_df$marker)
for (i in 1:length(samples)){
  a = consensus_results[consensus_results$sample==samples[i],]
  a = a[order(-a$auc),]
  a$rank = seq(1, nrow(a), 1)
  a = a[,c("marker","rank")]
  colnames(a) = c("marker", paste(samples[i], "rank", sep="_"))
  a$marker = as.character(a$marker)
  ranking_df = left_join(ranking_df, a, by="marker")
}

rt = ranking_df[,2:ncol(ranking_df)]
rownames(rt) = ranking_df$marker

continuous_rank = vector(mode="numeric")
for (i in 1:nrow(rt)){
  sum = 0
  for (j in 1:ncol(rt)){
    sum = sum + as.numeric(as.character(rt[i,j]))
  }
  continuous_rank[i] = sum/ncol(rt)
}

rt = data.frame(cbind(continuous_rank, rt))
#rt$continuous_rank = continuous_rank
rt = rt[order(rt$continuous_rank),]
rt = rt[rownames(rt) %in% near_consensus,]
write.csv(rt, "outputs/TableS10_human_comet_results.csv")

colnames(all_comet_results)[1] = "marker"
colnames(all_auc_results)[4] = "sample"

all_c = left_join(all_comet_results, all_auc_results, by=c("marker", "sample"))

df_c = all_c

df_c$logq = -log10(df_c$q_value)
df_c$loglogq = 0
for (i in 1:nrow(df_c)){
  df_c$loglogq[i] = min(10, log2(df_c$logq[i]+1))
}

markerstc = c("GYPC_negation", "LTB_negation", "CCR7_negation", "PDCD1", "FLT3LG_negation")

df_c_t = df_c
df_c_t$lbl = as.character(df_c_t$marker)
df_c_t$auc = as.numeric(as.character(df_c_t$auc))
df_c_t = df_c_t[df_c_t$auc>0.5,]
df_c_t$lbl[!(df_c_t$marker %in% markerstc)] = " "
df_c_t = df_c_t[order(df_c_t$lbl),]
pdf("outputs/Fig6D_bubble_plot.pdf", height=4, width=9)
g = ggplot(df_c_t, aes(sample, loglogq, color=lbl, size=auc)) + scale_color_manual(values=c("gray30", generatePalette(length(markerstc))))
g + geom_jitter(alpha=0.5) + theme_classic() + geom_hline(yintercept = log2(-log10(0.05) + 1), linetype="dashed") + geom_hline(yintercept = log2(-log10(10^-5) + 1), linetype="dashed")
dev.off()

HS = (read.csv("data/my_human_surface_genes.csv"))$x
all_auc_results = all_auc_results[!(grepl("negation", all_auc_results$marker)),]

all_auc_results$marker = as.character(all_auc_results$marker)
aucs_surface = all_auc_results[all_auc_results$marker %in% HS,]

aucs_inverse = aucs_surface
aucs_inverse$auc = 1 - aucs_inverse$auc

aucs_averaged = aucs_surface %>% group_by(marker) %>% summarize(mean_auc = mean(auc))
marks = aucs_averaged[aucs_averaged$marker %in% c("CCR7_negation", "GYPC_negation", "FLT3LG_negation", "LTB_negation"),]

# Plot the histogram
pdf("outputs/SuppFig4K_auc_histogram.pdf", height=4, width=6)
g = ggplot()
g = g + geom_histogram(aes(x=aucs_averaged$mean_auc), bins=40) + theme_classic() + xlim(c(0,1)) + xlab("AUC")
g = g + geom_vline(xintercept=(0.763), linetype="dashed", color=dark_blues[1]) 
g = g + geom_vline(xintercept=(0.696), linetype="dashed", color=dark_purples[1]) 
g = g + geom_vline(xintercept=(0.644), linetype="dashed", color=hot_pinks[1]) 
g = g + geom_vline(xintercept=(0.660), linetype="dashed", color=greens[1]) 
g + ylab("number of human surface markers")
dev.off()



