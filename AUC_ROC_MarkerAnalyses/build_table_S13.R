
samples = c("HumanIntegratedBlood_wLongitudinal", "K409_Blood",
            "K411_Blood", "K411_Blood_Longitudinal",
            "K468_Blood", "K468_Blood_Longitudinal",
            "K484_Blood")

interval = c(0.025, 0.975)

gate = vector(mode="character")
df = data.frame(matrix(nrow=7, ncol=0))

for (i in 1:length(samples)){
  boot = read.csv(paste0('outputs/', samples[i], "_specsens_bootstrapped.csv"))
  c = 2
  point_est_sens = vector(mode="numeric")
  ci_sens = vector(mode="character")
  p_sens = vector(mode="numeric")
  point_est_spec = vector(mode="numeric")
  ci_spec = vector(mode="character")
  p_spec = vector(mode="numeric")
  for (j in 1:((ncol(boot)-1)/4)){
    gate[j] = gsub("_.*$", "", colnames(boot)[c])
    point_est_sens[j] = round(mean(boot[,c]), 3)
    ci_sens[j] = paste(round(quantile(boot[,c], probs=interval)[1], 3), round(quantile(boot[,c], probs=interval)[2],3), sep="-")
    p_sens[j] = round(max(0.0001, nrow(boot[boot[,(c+2)] > point_est_sens[j],])/nrow(boot)), 4)
    c = c + 1
    point_est_spec[j] = round(mean(boot[,c]), 3)
    ci_spec[j] = paste(round(quantile(boot[,c], probs=interval)[1], 3), round(quantile(boot[,c], probs=interval)[2],3), sep="-")
    p_spec[j] = round(max(0.0001, nrow(boot[boot[,(c+2)] > point_est_spec[j],])/nrow(boot)), 4)
    c = c + 3
  }
  df = data.frame(cbind(df, point_est_sens, ci_sens, p_sens, point_est_spec, ci_spec, p_spec))
  colnames(df)[(ncol(df)-5):ncol(df)] = paste(samples[i], colnames(df)[(ncol(df)-5):ncol(df)], sep="_")
}
df = data.frame(cbind(gate, df))
df$gate = factor(df$gate, levels=c("CCR7", "FLT3LG", "GYPC", "LTB", "combo", "ands", "ors"))
df = df[order(df$gate),]
write.csv(df, "outputs/spec_sens_statistics.csv")









gate = vector(mode="character")
df = data.frame(matrix(nrow=4, ncol=0))

for (i in 1:length(samples)){
  boot = read.csv(paste0('outputs/', samples[i], "_aucs_bootstrapped.csv"))
  c = 2
  point_est = vector(mode="numeric")
  ci = vector(mode="character")
  p = vector(mode="numeric")
  for (j in 1:((ncol(boot)-1)/2)){
    gate[j] = gsub("_.*$", "", colnames(boot)[c])
    point_est[j] = round(mean(boot[,c]), 3)
    ci[j] = paste(round(quantile(boot[,c], probs=interval)[1], 3), round(quantile(boot[,c], probs=interval)[2],3), sep="-")
    p[j] = round(max(0.0001, nrow(boot[boot[,(c+1)] > point_est_sens[j],])/nrow(boot)), 4)
    c = c + 2
  }
  df = data.frame(cbind(df, point_est, ci, p))
  colnames(df)[(ncol(df)-2):ncol(df)] = paste(samples[i], colnames(df)[(ncol(df)-2):ncol(df)], sep="_")
}
df = data.frame(cbind(gate, df))
df$gate = factor(df$gate, levels=c("CCR7", "FLT3LG", "GYPC", "LTB"))
df = df[order(df$gate),]
write.csv(df, "outputs/auc_statistics.csv")

