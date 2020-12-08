source("../mycolors.R")

processCometFile <- function(id){
  d = read.csv(paste(id, "singletons_HS_full.csv", sep="_"))
  d$sample = id
  return(d)
}

processAUCFile <- function(id){
  d = read.csv(paste(id, "aucs_CHECK.csv", sep="_")) ##
  d$sample = id
  return(d)
}

getSigQs <- function(id){
  d = read.csv(paste(id, "_singletons_HS_full.csv", sep=""))
  d = d[d$q_value<0.05,]
  return(as.character(d$gene_1))
}

samples = c("k411Ac1","k411Bc1", "k468Ac1", "k468Bc1", "k409c1",  "k484c1_1124")

all_comet_results <- processCometFile(samples[1])
all_auc_results <- processAUCFile(samples[1])
all_auc_results = all_auc_results[order(-all_auc_results$auc),]
marker = as.character(all_auc_results$marker)
sig_qs = getSigQs(samples[1])

for (i in 2:length(samples)){
  all_comet_results = data.frame(rbind(all_comet_results, processCometFile(samples[i])))
  res = getSigQs(samples[i])
  if ("LAG3" %in% res){
    print(samples[i])
    print(all_comet_reuslts)
  }
  sig_qs = c(sig_qs, res)
  a = processAUCFile(samples[i])
  all_auc_results = data.frame(rbind(all_auc_results, a))
}

t = data.frame(table(sig_qs))
t = t[order(-t$Freq),]
t = data.frame(cbind(1:length(unique(sig_qs)), t)) 
colnames(t)[ncol(t)-1] = "marker"
near_consensus = as.character(t$marker[t$Freq >=4])
exclude = c("CD8A", "CCL4", "CCL5", "MIF_negation")
near_consensus = near_consensus[!(near_consensus %in% exclude)]

print(all_comet_results[all_comet_results$gene_1=="LAG3",])

consensus_results = all_auc_results[all_auc_results$marker %in% near_consensus,]

nrp = consensus_results %>% group_by(marker) %>% dplyr::summarise(K409_auc = auc[sample=="k409c1"], K411_auc = auc[sample=="k411Ac1"], K411longitudinal_auc = auc[sample=="k411Ac1"], K468_auc = auc[sample=="k468Ac1"], K468longitudinal_auc = auc[sample=="k468Bc1"], K484_auc = auc[sample=="k484c1_1124"], mean_auc = mean(auc))
nrp = nrp[order(-nrp$mean_auc),]

colnames(all_comet_results)[1] = "marker"
colnames(all_auc_results)[4] = "sample"

all_c = left_join(all_comet_results, all_auc_results, by=c("marker", "sample"))

df_c = all_c

df_c$logq = -log10(df_c$q_value)
df_c$loglogq = 0
for (i in 1:nrow(df_c)){
  df_c$loglogq[i] = min(10, log2(df_c$logq[i]+1))
}

markerstc = c("GYPC_negation", "LTB_negation", "CCR7_negation", "PDCD1", "FLT3LG_negation", "LGALS1", "KLRD1")##, "LGALS1")

df_c_t = df_c
df_c_t$lbl = as.character(df_c_t$marker)
df_c_t$auc = as.numeric(as.character(df_c_t$auc))
df_c_t = df_c_t[df_c_t$auc>0.5,]
df_c_t$lbl[!(df_c_t$marker %in% markerstc)] = " "
df_c_t = df_c_t[order(df_c_t$lbl),]
colors = c("#9a7fb9", "#52ab7d", "#dd78ae", "#f9cda4", "#f27071", "#7274b5", "#fae379")#db6fa9", "#7072b5", "#977cb7", "#fae16e", "#5eb186", "#f16c6e", "#f8c597")
#pdf("Fig6D_bubble_1127.pdf", height=4, width=9)

g = ggplot(df_c_t, aes(sample, loglogq, color=lbl, size=auc)) + scale_color_manual(values=c("gray30", colors))
g = g + geom_jitter(alpha=0.5) + theme_classic(base_size=15) + geom_hline(yintercept = log2(-log10(0.05) + 1), linetype="dashed") + geom_hline(yintercept = log2(-log10(10^-5) + 1), linetype="dashed")
g + xlab("") + ylab("loglogq") + scale_x_discrete(labels=c("k409c1" = "K409", "k411Ac1" = "K411",
                              "k411Bc1" = "K411\nlongitudinal", "k468Ac1" = "K468", "k468Bc1" = "K468\nlongitudinal", "k484c1_1124" = "K484"))
#dev.off()

HS = (read.csv("my_human_surface_genes.csv"))$x
all_auc_results = all_auc_results[!(grepl("negation", all_auc_results$marker)),]

all_auc_results$marker = as.character(all_auc_results$marker)
aucs_surface = all_auc_results[all_auc_results$marker %in% HS,]

aucs_inverse = aucs_surface
aucs_inverse$auc = 1 - aucs_inverse$auc

aucs_averaged = aucs_surface %>% group_by(marker) %>% summarize(mean_auc = mean(auc))
marks = aucs_averaged[aucs_averaged$marker %in% c("CCR7_negation", "GYPC_negation", "FLT3LG_negation", "LTB_negation"),]

pdf("figS3k_auc_histogram.pdf", height=4, width=6)
g = ggplot()
g = g + geom_histogram(aes(x=aucs_averaged$mean_auc), bins=40) + theme_classic() + xlim(c(0,1)) + xlab("AUC")
g = g + geom_vline(xintercept=(0.763), linetype="dashed", color=dark_blues[1]) 
g = g + geom_vline(xintercept=(0.696), linetype="dashed", color=dark_purples[1]) 
g = g + geom_vline(xintercept=(0.644), linetype="dashed", color=hot_pinks[1]) 
g = g + geom_vline(xintercept=(0.660), linetype="dashed", color=greens[1]) 
g + ylab("number of human surface markers")
dev.off()



