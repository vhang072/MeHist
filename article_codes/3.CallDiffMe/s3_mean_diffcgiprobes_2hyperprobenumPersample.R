cancerTypes <- c("LAML",
                 "ACC",
                 "BLCA",
                 "BRCA",
                 "CESC",
                 "CHOL",
                 "COAD",
                 "ESCA",
                 "GBM",
                 "HNSC",
                 "KICH",
                 "KIRC",
                 "KIRP",
                 "LGG",
                 "LIHC",
                 "LUAD",
                 "LUSC",
                 "PAAD",
                 "PCPG",
                 "PRAD",
                 "READ",
                 "SARC",
                 "SKCM",
                 "STAD",
                 "THCA",
                 "THYM",
                 "UCEC")
mean_rmna <- function(x){
  return(mean(x,na.rm = T)) 
}

sum_rmna <- function(x){
  return(sum(x,na.rm = T)) 
}

used_lines <- 1:length(cancerTypes)
cor_mean <- NULL
for (i in used_lines){
  print(cancerTypes[i])
  info0 <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_TCGA/_DiffCGIprobe_",cancerTypes[i],".txt"),
                     sep = "\t",header = T)
  info <- na.omit(info0)
  data <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Primary/_cgiProbe_TCGA-",cancerTypes[i],".txt"),
                     sep = "\t",header = T)
  rownames(data) <- data$probeid
  
  data_colrange <- 5:ncol(data)
  hyper <- info[(info$Mean_T > info$Mean_N) & (info$AberrantFrac_T_0.2 > 0.25) & (info$Mean_N < 0.3),]
  hypo <- info[(info$Mean_T < info$Mean_N) & (info$AberrantFrac_T_0.2 > 0.25) & (info$Mean_N > 0.7),]
  hyper_wide <- info[(info$Mean_T > info$Mean_N) & (info$AberrantFrac_T_0.2 > 0.05) & (info$Mean_N < 0.3),]
  hypo_wide <- info[(info$Mean_T < info$Mean_N) & (info$AberrantFrac_T_0.2 > 0.05) & (info$Mean_N > 0.7),]
  
  hyper_probe <- data[data$probeid %in% hyper$probeid,]
  hypo_probe <- data[data$probeid %in% hypo$probeid,]
  hyperwide_probe <- data[hyper_wide$probeid,data_colrange]
  hypowide_probe <- data[hypo_wide$probeid,data_colrange]
  
  
  hyper_mean <- apply(hyper_probe[,data_colrange],2,mean_rmna) - mean(hyper$Mean_N,na.rm = T)
  hypo_mean <- mean(hypo$Mean_N,na.rm = T) - apply(hypo_probe[,data_colrange],2,mean_rmna)
  tmp <-  matrix(rep(hyper_wide$Mean_N,ncol(hyperwide_probe)),nrow = nrow(hyperwide_probe),byrow = F)
  hyper_val <- hyperwide_probe - tmp
  hyper_num <- apply((hyper_val > 0.2),2,sum_rmna)
  tmp <-  matrix(rep(hypo_wide$Mean_N,ncol(hypowide_probe)),nrow = nrow(hypowide_probe),byrow = F)
  hypo_val <-  tmp - hypowide_probe
  hypo_num <- apply((hypo_val > 0.2),2,sum_rmna)
  
  #hyper corr
  corrs <- NULL
  for (idx in 1:nrow(hyperwide_probe)) {
    tmp.cor <- cor.test(as.numeric(hyperwide_probe[idx,]),hyper_mean,method ="pearson",exact = F)
    tmp.cor2 <- cor.test(as.numeric(hyperwide_probe[idx,]),hyper_mean,method ="spearman",exact = F)
    corrs <- rbind(corrs,c(tmp.cor$estimate,tmp.cor$p.value,tmp.cor2$estimate,tmp.cor2$p.value))
  }
  
  corrs <- as.data.frame(corrs)
  colnames(corrs) <- c("rho.Pearson","p.Pearson","rho.Spearman","p.Spearman")
  rownames(corrs) <- rownames(hyperwide_probe)
  write.table(corrs,paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_TCGA/_CGIprobeHyperCors_",cancerTypes[i],".txt"),
              sep = "\t",row.names = T,col.names = T, quote = F)
  
  #hypo corr
  corrs2 <- NULL
  for (idx in 1:nrow(hypowide_probe)) {
    tmp.cor <- cor.test(as.numeric(hypowide_probe[idx,]),hypo_mean,method ="pearson",exact = F)
    tmp.cor2 <- cor.test(as.numeric(hypowide_probe[idx,]),hypo_mean,method ="spearman",exact = F)
    corrs2 <- rbind(corrs2,c(tmp.cor$estimate,tmp.cor$p.value,tmp.cor2$estimate,tmp.cor2$p.value))
  }
  corrs2 <- as.data.frame(corrs2)
  colnames(corrs2) <- c("rho.Pearson","p.Pearson","rho.Spearman","p.Spearman")
  rownames(corrs2) <- rownames(hypowide_probe)
  write.table(corrs2,paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_TCGA/_CGIprobeHypoCors_",cancerTypes[i],".txt"),
              sep = "\t",row.names = T,col.names = T, quote = F)
  
  x1 <- cor((hyper_mean),t(hyper_probe[,data_colrange]))
  hist(x1,breaks = seq(-1,1,0.05))
  x2 <- cor((hypo_mean),t(hypo_probe[,data_colrange]))
  hist(x2,breaks = seq(-1,1,0.05))
  x1[is.na(x1)] <- 0
  x2[is.na(x2)] <- 0
  
  cor_mean <- rbind(cor_mean,c(mean(x1),mean(x2)))
  mean_vals <- rbind(hyper_mean,hypo_mean,hyper_num,hypo_num)
  mean_vals <- cbind(c("hyper_mean","hypo_mean","hyperprobe_num","hypoprobe_num"),mean_vals)
  mean_vals <- as.data.frame(mean_vals)
  colnames(mean_vals) <- c("mean_type",colnames(data)[data_colrange])
  write.table(mean_vals,paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_TCGA/_DiffCGIprobeMean_",cancerTypes[i],".txt"),
              sep = "\t",row.names = F,col.names = T, quote = F)
}
cor_mean <- cbind(cancerTypes,cor_mean)
cor_mean <- as.data.frame(cor_mean)
colnames(cor_mean) <- c("cancer_types","hyper_cor_mean","hypo_cor_mean")

write.table(cor_mean,paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_TCGA/__DiffCGIprobeCorSummary.txt"),
            sep = "\t",row.names = F,col.names = T, quote = F)
