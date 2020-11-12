
# Load libraries
library(DescTools)
library(ggplot2)
library(ggpubr)
library(dplyr)


processAUCFile <- function(id){
  d = read.csv(paste(id, "aucs.csv", sep="_"))  ##
  d$sample = id
  
  # Keep all positive (the non-negation) markers
  d = d[!(grepl("negation", d$marker)),]
  d$auc = as.numeric(as.character(d$auc))
  d = d[order(-d$auc),]
  write.csv(d, paste(id, "aucs_ranked.csv", sep="_"))
  return(d)
}

setwd('data/')
HS = (read.csv("my_human_surface_genes.csv"))$x

K409 = processAUCFile("K409_Blood")
K411_first = processAUCFile("K411_Blood")
K411_second = processAUCFile("K411_Blood_Longitudinal")
K468_first = processAUCFile("K468_Blood")
K468_second = processAUCFile("K468_Blood_Longitudinal")
K484 = processAUCFile("K484_Blood")
setwd('../')

samples = list(K409, K484, K468_first, K468_second, K411_first, K411_second)
names = c("K409", "K484", "K468_first", "K468_second", "K411_first", "K411_second")

mdf = data.frame()
xsample = vector(mode="character")
ysample = vector(mode="character")
for (i in 1:(length(samples)-1)){
  for (j in (i+1):length(samples)){
    if (i!=j){
      si = samples[[i]]
      sj = samples[[j]]
      mini = inner_join(si, sj, by="marker")
      xsample = c(xsample, rep(names[i], nrow(mini)))
      ysample = c(ysample, rep(names[j], nrow(mini)))
      if (substr(names[i], 1, 4)==substr(names[j], 1, 4)){
        mini$type = "within"
      } else {
        mini$type = "across"
      }
      mdf = data.frame(rbind(mdf, mini))
    }
  }
}

mdf$surface = 0
mdf$surface[mdf$marker %in% HS] = 1
mdf$surface = factor(mdf$surface, levels=c(0,1))

set.seed = 27
for (i in 1:nrow(mdf)){
  r = runif(1, 0, 1)
  if (r > 0.5){
    y = mdf$auc.y[i]
    mdf$auc.y[i] = mdf$auc.x[i]
    mdf$auc.x[i] = y
    ys = ysample[i]
    ysample[i] = xsample[i]
    xsample[i] = ys
  }
}

mdf$ysample = ysample
mdf$xsample = xsample


cyto = mdf[mdf$surface==0,]
surf = mdf[mdf$surface==1,]
cytow = cyto[cyto$type=="within",]
cytoa = cyto[cyto$type=="across",]
surfw = surf[surf$type=="within",]
surfa = surf[surf$type=="across",]

pdf("outputs/Fig5C_auc_scatter.pdf", height=8, width=12)
g = ggplot()
g = g + geom_point(aes(surfa$auc.x, surfa$auc.y), pch=1, color="black", size=2, stroke=1)
g = g + geom_point(aes(surfa$auc.x, surfa$auc.y), color="aquamarine4", alpha=0.7)
g = g + geom_point(aes(cytoa$auc.x, cytoa$auc.y), color="aquamarine4", alpha=0.7) 
g = g + geom_point(aes(surfw$auc.x, surfw$auc.y), pch=1, color="black", size=2, stroke=1)
g = g + geom_point(aes(surfw$auc.x, surfw$auc.y), color="darkmagenta", alpha=0.7) 
g = g + geom_point(aes(cytow$auc.x, cytow$auc.y), color="darkmagenta", alpha=0.7) + xlim(c(0.08,1)) + ylim(c(0.08,0.9))
g = g + theme_classic() + geom_abline(intercept = 0, slope=1, linetype="dashed") + stat_cor() + xlab("AUC on first sample") + ylab("AUC on second sample")
g + stat_cor(aes(x=c(surfa$auc.x, cytoa$auc.x), y=c(surfa$auc.y, cytoa$auc.y)), color="aquamarine4", label.x=0.1, label.y=0.85) + stat_cor(aes(x=c(surfw$auc.x, cytow$auc.x), y=c(surfw$auc.y, cytow$auc.y)), color="darkmagenta", label.x=0.1, label.y=0.82)
dev.off()

pdf("outputs/auc_scatter_onlysurface.pdf", height=8, width=12)
g = ggplot(surf, aes(x=auc.x, y=auc.y))
g + geom_point() + stat_cor() + theme_classic()
dev.off()





dev.off()