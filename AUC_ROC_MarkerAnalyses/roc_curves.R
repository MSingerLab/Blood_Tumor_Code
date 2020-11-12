
# Load libraries
library(Seurat)

# Load the function 'Keep_AB_Cells' which subsets the rds object
# to keep only cells with at least one alpha and one beta chain
source("Keep_AB_Cells.R")
source('kaitlyncolors.R')

##first, the mouse ROC curve for exhuastion markers (Fig 1i)

so = Keep_AB_Cells('data/MouseIntegratedBlood.rds')
counts = so@assays$RNA
md = so@meta.data


resol = 0.05
markerstc = c("Btla", "Ctla4", "Havcr2", "Lag3", "Cd160", "Pdcd1", "Tigit")
mouse_aucs = vector(mode="numeric")
neg = rep(FALSE, 7)
df = data.frame(matrix(ncol=3))
colnames(df) = c("marker", "x", "y")
for (j in 1:length(markerstc)){
  ind = which(as.character(rownames(counts))==gsub("_.*$", "", markerstc[j]))
  mdt = md
  mdt$TM = as.factor(mdt$TM)
  mdt$marker = as.numeric(as.character(counts[ind,]))
  mdt_match = mdt[mdt$TM==1,]
  mdt_nomatch = mdt[mdt$TM!=1,]
  gates = c(seq(min(mdt$marker), max(mdt$marker), length=(1/resol)), quantile(mdt$marker, probs=seq(resol,(1-resol),by=resol)))
  x = vector(mode="numeric")
  y = vector(mode="numeric")
  for (k in 1:length(gates)){
    if (!neg[j]){
      x[k] = 1 - (nrow(mdt_nomatch[mdt_nomatch$marker < gates[k],])/nrow(mdt_nomatch))
      y[k] = nrow(mdt_match[mdt_match$marker >= gates[k],])/nrow(mdt_match)
    } else {
      x[k] = 1 - (nrow(mdt_nomatch[mdt_nomatch$marker > gates[k],])/nrow(mdt_nomatch))
      y[k] = nrow(mdt_match[mdt_match$marker <= gates[k],])/nrow(mdt_match)
    }
  }
  x = c(x, 0, 1)
  y = c(y, 0, 1)
  mouse_aucs = c(mouse_aucs, AUC(x,y))
  print(paste(markerstc[j], AUC(x,y), sep=": "))
  marker = rep(markerstc[j], length(x))
  minidf = data.frame(cbind(marker, x, y))
  df = data.frame(rbind(df, minidf))
}
# Remoave rows with NAs that were created during initialization of data.frame
df = na.omit(df)
df$marker = factor(df$marker, levels=markerstc)
df$x = as.numeric(as.character(df$x))
df$y = as.numeric(as.character(df$y))
pdf(paste("outputs/MouseIntBlood", paste(markerstc, collapse='_'), "ROC_curves_check.pdf", sep="_"), height=4, width=4)
g = ggplot(df, aes(x, y, color=marker))
g = g + geom_line(size=0.5, show.legend = TRUE) + theme_classic() + scale_color_manual(values = generatePalette(length(markerstc))) ##scale_color_manual(values=c(generatePalette(length(markerstc))))#+ scale_color_manual(values = c(generatePalette((length(markerstc)-1), 27), "black"))
g = g + xlab("1 - Specificity") + ylab("Sensitivity") ##+ xlim(c(0,1)) + ylim(c(0,1))
print(g + geom_abline(slope=1, intercept=0, linetype="dashed"))
dev.off()

write.csv(data.frame(cbind(markerstc, mouse_aucs)), "outputs/Mouse_exhuastion_aucs_check.csv")

##next, ROC curves for good markers in mouse (Fig2h)
markers = c("NKG2D_TotalSeqC", "CD39_TotalSeqC", "CX3CR1_TotalSeqC")
neg = rep(FALSE, 3)
samples = c("M4_Blood", "M5_Blood")
aucs = vector(mode="numeric")
resol <- 0.05
for (i in 1:length(samples)){
  rds = Keep_AB_Cells(paste0('data/', samples[i], '.rds'))
  md = rds@meta.data
  df = data.frame(matrix(ncol=3))
  colnames(df) = c("marker", "x", "y")
  for (j in 1:length(markers)){
    ind = which(colnames(md)==markers[j])
    mdt = md
    mdt$TM = as.factor(mdt$TM)
    mdt$marker = mdt[,ind]
    mdt_match = mdt[mdt$TM==1,]
    mdt_nomatch = mdt[mdt$TM!=1,]
    gates = c(seq(min(mdt$marker), max(mdt$marker), length=(1/resol)), quantile(mdt$marker, probs=seq(resol,(1-resol),by=resol)))
    x = vector(mode="numeric")
    y = vector(mode="numeric")
    for (k in 1:length(gates)){
      if (!neg[j]){
        x[k] = 1 - (nrow(mdt_nomatch[mdt_nomatch$marker < gates[k],])/nrow(mdt_nomatch))
        y[k] = nrow(mdt_match[mdt_match$marker >= gates[k],])/nrow(mdt_match)
      } else {
        x[k] = 1 - (nrow(mdt_nomatch[mdt_nomatch$marker > gates[k],])/nrow(mdt_nomatch))
        y[k] = nrow(mdt_match[mdt_match$marker <= gates[k],])/nrow(mdt_match)
      }
    }
    x = c(x, 0, 1)
    y = c(y, 0, 1)
    aucs = c(aucs, AUC(x,y))
    print(paste(markers[j], AUC(x,y), sep=": "))
    marker = rep(markers[j], length(x))
    minidf = data.frame(cbind(marker, x, y))
    df = data.frame(rbind(df, minidf))
  }
  # Remoave rows with NAs that were created during initialization of data.frame
  df = na.omit(df)
  df$marker = factor(df$marker, levels=markers)
  df$x = as.numeric(as.character(df$x))
  df$y = as.numeric(as.character(df$y))
  pdf(paste0('outputs/', samples[i], "_good_markers_ROC_curves.pdf"), height=4, width=4)
  g = ggplot(df, aes(x, y, color=marker))
  g = g + geom_line(size=0.5, show.legend = TRUE) + theme_classic() + scale_color_manual(values = generatePalette(length(markers))) ##scale_color_manual(values=c(generatePalette(length(markerstc))))#+ scale_color_manual(values = c(generatePalette((length(markerstc)-1), 27), "black"))
  g = g + xlab("1 - Specificity") + ylab("Sensitivity") ##+ xlim(c(0,1)) + ylim(c(0,1))
  print(g + geom_abline(slope=1, intercept=0, linetype="dashed"))
  dev.off()
}

##next, exhuastion ROC curves for human samples (Fig5b-c)
markers = c("BTLA", "CTLA4", "HAVCR2", "LAG3", "CD160", "PDCD1", "TIGIT")
neg = rep(FALSE, 7)
samples = c("K409_Blood", "K411_Blood", "K411_Blood_Longitudinal", "K468_Blood", "K468_Blood_Longitudinal", "K484_Blood")
aucs = vector(mode="numeric")
resol <- 0.05
for (i in 1:length(samples)){
  rds <- Keep_AB_Cells(paste0('data/', samples[i], '.rds'))
  md = rds@meta.data
  counts = rds@assays$RNA@counts
  df = data.frame(matrix(ncol=3))
  colnames(df) = c("marker", "x", "y")
  for (j in 1:length(markers)){
    ind = which(rownames(counts)==markers[j])
    mdt = md
    mdt$TM = as.factor(mdt$TM)
    mdt$marker = as.numeric(as.character(counts[ind,]))
    mdt_match = mdt[mdt$TM==1,]
    mdt_nomatch = mdt[mdt$TM!=1,]
    gates = c(seq(min(mdt$marker), max(mdt$marker), length=(1/resol)), quantile(mdt$marker, probs=seq(resol,(1-resol),by=resol)))
    x = vector(mode="numeric")
    y = vector(mode="numeric")
    for (k in 1:length(gates)){
      if (!neg[j]){
        x[k] = 1 - (nrow(mdt_nomatch[mdt_nomatch$marker < gates[k],])/nrow(mdt_nomatch))
        y[k] = nrow(mdt_match[mdt_match$marker >= gates[k],])/nrow(mdt_match)
      } else {
        x[k] = 1 - (nrow(mdt_nomatch[mdt_nomatch$marker > gates[k],])/nrow(mdt_nomatch))
        y[k] = nrow(mdt_match[mdt_match$marker <= gates[k],])/nrow(mdt_match)
      }
    }
    x = c(x, 0, 1)
    y = c(y, 0, 1)
    aucs = c(aucs, AUC(x,y))
    print(paste(markers[j], AUC(x,y), sep=": "))
    marker = rep(markers[j], length(x))
    minidf = data.frame(cbind(marker, x, y))
    df = data.frame(rbind(df, minidf))
  }
  # Remoave rows with NAs that were created during initialization of data.frame
  df = na.omit(df)
  df$marker = factor(df$marker, levels=markers)
  df$x = as.numeric(as.character(df$x))
  df$y = as.numeric(as.character(df$y))
  pdf(paste0('outputs/', samples[i], "_exhaustion_markers_ROC_curves.pdf"), height=4, width=4)
  g = ggplot(df, aes(x, y, color=marker))
  g = g + geom_line(size=0.5, show.legend = TRUE) + theme_classic() + scale_color_manual(values = generatePalette(length(markers))) ##scale_color_manual(values=c(generatePalette(length(markerstc))))#+ scale_color_manual(values = c(generatePalette((length(markerstc)-1), 27), "black"))
  g = g + xlab("1 - Specificity") + ylab("Sensitivity") ##+ xlim(c(0,1)) + ylim(c(0,1))
  print(g + geom_abline(slope=1, intercept=0, linetype="dashed"))
  dev.off()
}


##finally, ROC curves for the good markers in human samples (Fig5e)
markers = c("CCR7", "FLT3LG", "GYPC", "LTB")
neg = rep(TRUE, 4)
samples = c("K409_Blood", "K411_Blood", "K411_Blood_Longitudinal", "K468_Blood", "K468_Blood_Longitudinal", "K484_Blood")
aucs = vector(mode="numeric")
resol <- 0.05
for (i in 1:length(samples)){
  rds <- Keep_AB_Cells(paste0('data/', samples[i], '.rds'))
  md = rds@meta.data
  counts = rds@assays$RNA@counts
  df = data.frame(matrix(ncol=3))
  colnames(df) = c("marker", "x", "y")
  for (j in 1:length(markers)){
    ind = which(rownames(counts)==markers[j])
    mdt = md
    mdt$TM = as.factor(mdt$TM)
    mdt$marker = as.numeric(as.character(counts[ind,]))
    mdt_match = mdt[mdt$TM==1,]
    mdt_nomatch = mdt[mdt$TM!=1,]
    gates = c(seq(min(mdt$marker), max(mdt$marker), length=(1/resol)), quantile(mdt$marker, probs=seq(resol,(1-resol),by=resol)))
    x = vector(mode="numeric")
    y = vector(mode="numeric")
    for (k in 1:length(gates)){
      if (!neg[j]){
        x[k] = 1 - (nrow(mdt_nomatch[mdt_nomatch$marker < gates[k],])/nrow(mdt_nomatch))
        y[k] = nrow(mdt_match[mdt_match$marker >= gates[k],])/nrow(mdt_match)
      } else {
        x[k] = 1 - (nrow(mdt_nomatch[mdt_nomatch$marker > gates[k],])/nrow(mdt_nomatch))
        y[k] = nrow(mdt_match[mdt_match$marker <= gates[k],])/nrow(mdt_match)
      }
    }
    x = c(x, 0, 1)
    y = c(y, 0, 1)
    aucs = c(aucs, AUC(x,y))
    print(paste(markers[j], AUC(x,y), sep=": "))
    marker = rep(markers[j], length(x))
    minidf = data.frame(cbind(marker, x, y))
    df = data.frame(rbind(df, minidf))
  }
  # Remoave rows with NAs that were created during initialization of data.frame
  df = na.omit(df)
  df$marker = factor(df$marker, levels=markers)
  df$x = as.numeric(as.character(df$x))
  df$y = as.numeric(as.character(df$y))
  pdf(paste0('outputs/', samples[i], "_good_markers_ROC_curves.pdf"), height=4, width=4)
  g = ggplot(df, aes(x, y, color=marker))
  g = g + geom_line(size=0.5, show.legend = TRUE) + theme_classic() + scale_color_manual(values = generatePalette(length(markers))) ##scale_color_manual(values=c(generatePalette(length(markerstc))))#+ scale_color_manual(values = c(generatePalette((length(markerstc)-1), 27), "black"))
  g = g + xlab("1 - Specificity") + ylab("Sensitivity") ##+ xlim(c(0,1)) + ylim(c(0,1))
  print(g + geom_abline(slope=1, intercept=0, linetype="dashed"))
  dev.off()
}
