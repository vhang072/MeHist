cancerTypes1 <- c(
                 "KIRC",
                 "KIRP",
                 "PCPG",
                 "KICH",
                 "THYM",
                 "LAML",
                 "THCA")
cancerTypes2 <- c(
                  "PRAD",
                  "PAAD",
                  "BRCA",
                  "LUAD",
                  "STAD",
                  "SARC",
                  "ACC")

cancerTypes3 <- c("CESC",
                  "LIHC",
                  "HNSC",
                  "BLCA",
                  "CHOL",
                  "UCEC",
                  "GBM",
                  "LUSC",
                  "LGG",
                  "COAD",
                  "ESCA",
                  "READ",
                  "SKCM")

library(stringr)
library(AdaptGauss)
library(dplyr)
library(survival)
library(survminer)
library(gtsummary)


cdr <- read.table("Y:/4.basic_data/TCGA_PancanAtlas/clinical/TCGA-CDR.txt",sep = "\t",header = T)
rn <- NULL 
for (x in 1:nrow(cdr)){
  rn <- c(rn, str_replace_all(cdr$bcr_patient_barcode[x],'-',"."))
}
rownames(cdr) <- rn


merged_mepersample <- NULL
median_num <- NULL
pdf("I:/Projects/P7.Methy-PanCancer/pics/tcga_survival_purity/typessurvival_puritylowepimut.pdf",height = 6,width = 6)
pancancer_out <- NULL

#type1
for (idx in 1:length(cancerTypes1)){
  tmp_out <- NULL
  data <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_TCGA/ABSOLUTE_MeanCGIprobe/_SigPur.",cancerTypes1[idx],".CGIMean_ABSOLUTE.txt"),
                     sep = "\t",header = T)
  coln <- NULL
  colx <- NULL
  for (x in 1:nrow(data)){
    if (substr(rownames(data)[x],1,12) %in% coln == F){
      coln <- c(coln, substr(rownames(data)[x],1,12))
      colx <- c(colx, x)
    }
  }
  data_f <- data[colx,]
  rownames(data_f) <- coln
  overlap_sps <- intersect(rownames(cdr),rownames(data_f))
  
  data_u <- data_f[overlap_sps,]
  cdr_u <- cdr[overlap_sps,]
  
  cdr_data <- cbind(cdr_u[,26:33],data_u[,c("Hyper.Num","Hypo.Num","Purity","label")])
  
  
  #m1 <- median(cdr_data$Hyper.Num)
  #std1 <- sqrt(var(cdr_data$Hyper.Num))
  #m2 <- median(cdr_data$Hypo.Num)
  #std2 <- sqrt(var(cdr_data$Hypo.Num))  
  #out1 <- EMGauss(cdr_data$Hyper.Num,K = 2,Means = c(0,100),SDs = c(20,20),MaxNumberofIterations =20)
  #thred1 <- quantile(cdr_data$Hyper.Num,out1$Weights[1])
  #out2 <- EMGauss(cdr_data$Hypo.Num,K = 2,Means = c(0,100),SDs = c(20,20),MaxNumberofIterations =20)
  #thred2 <- quantile(cdr_data$Hypo.Num,out2$Weights[1])
  #out3 <- EMGauss(cdr_data$Hypo.Num *3 + cdr_data$Hyper.Num,K = 2,Means = c(0,100),SDs = c(20,20),MaxNumberofIterations =20)
  #thred3 <- quantile(cdr_data$Hypo.Num *3 + cdr_data$Hyper.Num,out3$Weights[1])
  
  cdr_data$me_label <- "HighEpiMut"
  thrednew <- quantile(cdr_data$Hypo.Num *2 + cdr_data$Hyper.Num, 0.7)
  cdr_data$me_label[cdr_data$Hypo.Num *2 + cdr_data$Hyper.Num <= thrednew] <- "LowEpiMut"
  #cdr_data$me_label[(cdr_data$Hyper.Num <= (m1+std1)) & (cdr_data$Hypo.Num <= (m2+std2))] <- "LowEpiMut"
  cdr_data$me_label <- factor(cdr_data$me_label,levels = c("LowEpiMut","HighEpiMut"))
  write.table(cdr_data,paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Survival_purity/typesDichotomy_",cancerTypes1[idx],".txt"),
               sep = "\t",quote = F)
  
  #os
  cdr_data_os <- cdr_data[(is.na(cdr_data$OS.time) == F),]#& (cdr_data$label == "Included")
  tnum <- sort(table(cdr_data_os$me_label),decreasing = T)
  if ((tnum[1]>=5) & (tnum[2]>=5)) {
    coxfit <- coxph(Surv(OS.time,OS) ~ me_label, data=cdr_data_os)
    x.out <- summary(coxfit)
    tmp_out <- c(tmp_out,x.out$conf.int[1],x.out$sctest[3])
    
    fit <- survfit(Surv(OS.time,OS) ~ me_label, data=cdr_data_os)
    surv_pvalue(fit)
    p <- ggsurvplot(fit,
                    pval = TRUE, conf.int = F,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
                    
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    palette = c("#0076BB", "#E69F00", "#56B4E9", "#800000"),
                    pval.method = T,title = paste0(cancerTypes1[idx],"-","OS"))
    print(p)
  } else {
    tmp_out <- c(tmp_out,NA,NA)
  }
  #dss
  cdr_data_dss <- cdr_data[(is.na(cdr_data$DSS.time) == F),] #& (cdr_data$label == "Included")
  tnum <- sort(table(cdr_data_dss$me_label),decreasing = T)
  if ((tnum[1]>=5) & (tnum[2]>=5)) {
    coxfit <- coxph(Surv(DSS.time,DSS) ~ me_label, data=cdr_data_dss)
    x.out <- summary(coxfit)
    tmp_out <- c(tmp_out,x.out$conf.int[1],x.out$sctest[3])
    
    fit <- survfit(Surv(DSS.time,DSS) ~ me_label, data=cdr_data_dss)
    p <- ggsurvplot(fit,
                    pval = TRUE, conf.int = F,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
                    
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    palette = c("#0076BB", "#E69F00", "#56B4E9", "#800000"),
                    pval.method = T,title = paste0(cancerTypes1[idx],"-","DSS"))
    print(p)
  } else {
    tmp_out <- c(tmp_out, NA,NA)
  }
  #dfi
  cdr_data_dfi <- cdr_data[(is.na(cdr_data$DFI.time) == F) ,] # & (cdr_data$label == "Included")
  tnum <- sort(table(cdr_data_dfi$me_label),decreasing = T)
  if ((tnum[1]>=5) & (tnum[2]>=5)) {
    coxfit <- coxph(Surv(DFI.time,DFI) ~ me_label, data=cdr_data_dfi)
    x.out <- summary(coxfit)
    tmp_out <- c(tmp_out,x.out$conf.int[1],x.out$sctest[3])
    
    fit <- survfit(Surv(DFI.time,DFI) ~ me_label, data=cdr_data_dfi)
    p <- ggsurvplot(fit,
                    pval = TRUE, conf.int = F,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
                    
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    palette = c("#0076BB", "#E69F00", "#56B4E9", "#800000"),
                    pval.method = T,title = paste0(cancerTypes1[idx],"-","DFI"))
    print(p)
  }  else {
    tmp_out <- c(tmp_out, NA,NA)
  }
  #pfi
  cdr_data_pfi <- cdr_data[(is.na(cdr_data$PFI.time) == F) ,] #& (cdr_data$label == "Included")
  tnum <- sort(table(cdr_data_pfi$me_label),decreasing = T)
  if ((tnum[1]>=5) & (tnum[2]>=5)) {
    coxfit <- coxph(Surv(PFI.time,PFI) ~ me_label, data=cdr_data_pfi)
    x.out <- summary(coxfit)
    tmp_out <- c(tmp_out,x.out$conf.int[1], x.out$sctest[3])
    
    fit <- survfit(Surv(PFI.time,PFI) ~ me_label, data=cdr_data_pfi)
    p <- ggsurvplot(fit,
                    pval = TRUE, conf.int = F,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
                    
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    palette = c("#0076BB", "#E69F00", "#56B4E9", "#800000"),
                    pval.method = T,title = paste0(cancerTypes1[idx],"-","PFI"))
    print(p)
  }  else {
    tmp_out <- c(tmp_out,NA,NA)
  }
  pancancer_out <- rbind(pancancer_out,c(tmp_out,max(nrow(cdr_data_os),
                                                     nrow(cdr_data_dss),
                                                     nrow(cdr_data_dfi),
                                                     nrow(cdr_data_pfi))))
  
}
#type2
for (idx in 1:length(cancerTypes2)){
  tmp_out <- NULL
  data <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_TCGA/ABSOLUTE_MeanCGIprobe/_SigPur.",cancerTypes2[idx],".CGIMean_ABSOLUTE.txt"),
                     sep = "\t",header = T)
  coln <- NULL
  colx <- NULL
  for (x in 1:nrow(data)){
    if (substr(rownames(data)[x],1,12) %in% coln == F){
      coln <- c(coln, substr(rownames(data)[x],1,12))
      colx <- c(colx, x)
    }
  }
  data_f <- data[colx,]
  rownames(data_f) <- coln
  overlap_sps <- intersect(rownames(cdr),rownames(data_f))
  
  data_u <- data_f[overlap_sps,]
  cdr_u <- cdr[overlap_sps,]
  
  cdr_data <- cbind(cdr_u[,26:33],data_u[,c("Hyper.Num","Hypo.Num","Purity","label")])
  cdr_data <- cdr_data[cdr_data$Purity >= min(min(cdr_data$Purity[cdr_data$label == "Included"]),0.5),]
  
  #m1 <- median(cdr_data$Hyper.Num)
  #std1 <- sqrt(var(cdr_data$Hyper.Num))
  #m2 <- median(cdr_data$Hypo.Num)
  #std2 <- sqrt(var(cdr_data$Hypo.Num))  
  out3 <- EMGauss(cdr_data$Hypo.Num *2 + cdr_data$Hyper.Num,K = 2,Means = c(0,100),SDs = c(20,20),MaxNumberofIterations =20)
  thred3 <- quantile(cdr_data$Hypo.Num *2 + cdr_data$Hyper.Num,out3$Weights[1])
  
  cdr_data$me_label <- "HighEpiMut"
  thrednew <- quantile(cdr_data$Hypo.Num *4 + cdr_data$Hyper.Num, 0.3)
  cdr_data$me_label[cdr_data$Hypo.Num *4 + cdr_data$Hyper.Num <= thrednew] <- "LowEpiMut"
  
  #cdr_data$me_label <- "HighEpiMut"
  #cdr_data$me_label[cdr_data$Hypo.Num *2 + cdr_data$Hyper.Num <= 8000] <- "LowEpiMut"
  #cdr_data$me_label[(cdr_data$Hyper.Num <= m1) & (cdr_data$Hypo.Num <= m2)] <- "LowEpiMut"
  # for (xx in 1:nrow(cdr_data)){
  #   if (cdr_data$Hyper.Num[xx] <= (m1 -m1/m2*cdr_data$Hypo.Num[xx])){
  #     cdr_data$me_label[xx] <- "LowEpiMut"
  #   }
  # }
  
  cdr_data$me_label <- factor(cdr_data$me_label,levels = c("LowEpiMut","HighEpiMut"))
  write.table(cdr_data,paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Survival_purity/typesDichotomy_",cancerTypes2[idx],".txt"),
              sep = "\t",quote = F)
  
  #os
  cdr_data_os <- cdr_data[(is.na(cdr_data$OS.time) == F),]#& (cdr_data$label == "Included")
  tnum <- sort(table(cdr_data_os$me_label),decreasing = T)
  if ((tnum[1]>=5) & (tnum[2]>=5)) {
    coxfit <- coxph(Surv(OS.time,OS) ~ me_label, data=cdr_data_os)
    x.out <- summary(coxfit)
    tmp_out <- c(tmp_out,x.out$conf.int[1],x.out$sctest[3])
    
    fit <- survfit(Surv(OS.time,OS) ~ me_label, data=cdr_data_os)
    surv_pvalue(fit)
    p <- ggsurvplot(fit,
                    pval = TRUE, conf.int = F,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
                    
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    palette = c("#0076BB", "#E69F00", "#56B4E9", "#800000"),
                    pval.method = T,title = paste0(cancerTypes2[idx],"-","OS"))
    print(p)
  } else {
    tmp_out <- c(tmp_out,NA,NA)
  }
  #dss
  cdr_data_dss <- cdr_data[(is.na(cdr_data$DSS.time) == F),] #& (cdr_data$label == "Included")
  tnum <- sort(table(cdr_data_dss$me_label),decreasing = T)
  if ((tnum[1]>=5) & (tnum[2]>=5)) {
    coxfit <- coxph(Surv(DSS.time,DSS) ~ me_label, data=cdr_data_dss)
    x.out <- summary(coxfit)
    tmp_out <- c(tmp_out,x.out$conf.int[1],x.out$sctest[3])
    
    fit <- survfit(Surv(DSS.time,DSS) ~ me_label, data=cdr_data_dss)
    p <- ggsurvplot(fit,
                    pval = TRUE, conf.int = F,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
                    
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    palette = c("#0076BB", "#E69F00", "#56B4E9", "#800000"),
                    pval.method = T,title = paste0(cancerTypes2[idx],"-","DSS"))
    print(p)
  } else {
    tmp_out <- c(tmp_out, NA,NA)
  }
  #dfi
  cdr_data_dfi <- cdr_data[(is.na(cdr_data$DFI.time) == F) ,] # & (cdr_data$label == "Included")
  tnum <- sort(table(cdr_data_dfi$me_label),decreasing = T)
  if ((tnum[1]>=5) & (tnum[2]>=5)) {
    coxfit <- coxph(Surv(DFI.time,DFI) ~ me_label, data=cdr_data_dfi)
    x.out <- summary(coxfit)
    tmp_out <- c(tmp_out,x.out$conf.int[1],x.out$sctest[3])
    
    fit <- survfit(Surv(DFI.time,DFI) ~ me_label, data=cdr_data_dfi)
    p <- ggsurvplot(fit,
                    pval = TRUE, conf.int = F,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
                    
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    palette = c("#0076BB", "#E69F00", "#56B4E9", "#800000"),
                    pval.method = T,title = paste0(cancerTypes2[idx],"-","DFI"))
    print(p)
  }  else {
    tmp_out <- c(tmp_out, NA,NA)
  }
  #pfi
  cdr_data_pfi <- cdr_data[(is.na(cdr_data$PFI.time) == F) ,] #& (cdr_data$label == "Included")
  tnum <- sort(table(cdr_data_pfi$me_label),decreasing = T)
  if ((tnum[1]>=5) & (tnum[2]>=5)) {
    coxfit <- coxph(Surv(PFI.time,PFI) ~ me_label, data=cdr_data_pfi)
    x.out <- summary(coxfit)
    tmp_out <- c(tmp_out,x.out$conf.int[1], x.out$sctest[3])
    
    fit <- survfit(Surv(PFI.time,PFI) ~ me_label, data=cdr_data_pfi)
    p <- ggsurvplot(fit,
                    pval = TRUE, conf.int = F,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
                    
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    palette = c("#0076BB", "#E69F00", "#56B4E9", "#800000"),
                    pval.method = T,title = paste0(cancerTypes2[idx],"-","PFI"))
    print(p)
  }  else {
    tmp_out <- c(tmp_out,NA,NA)
  }
  pancancer_out <- rbind(pancancer_out,c(tmp_out,max(nrow(cdr_data_os),
                                                     nrow(cdr_data_dss),
                                                     nrow(cdr_data_dfi),
                                                     nrow(cdr_data_pfi))))
}
  



#type3
for (idx in 1:length(cancerTypes3)){
  tmp_out <- NULL
  data <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_TCGA/ABSOLUTE_MeanCGIprobe/_SigPur.",cancerTypes3[idx],".CGIMean_ABSOLUTE.txt"),
                     sep = "\t",header = T)
  coln <- NULL
  colx <- NULL
  for (x in 1:nrow(data)){
    if (substr(rownames(data)[x],1,12) %in% coln == F){
      coln <- c(coln, substr(rownames(data)[x],1,12))
      colx <- c(colx, x)
    }
  }
  data_f <- data[colx,]
  rownames(data_f) <- coln
  overlap_sps <- intersect(rownames(cdr),rownames(data_f))
  
  data_u <- data_f[overlap_sps,]
  cdr_u <- cdr[overlap_sps,]
  
  cdr_data <- cbind(cdr_u[,26:33],data_u[,c("Hyper.Num","Hypo.Num","Purity","label")])
  cdr_data <- cdr_data[cdr_data$Purity >= min(min(cdr_data$Purity[cdr_data$label == "Included"]),0.5),]
  
  #m1 <- median(cdr_data$Hyper.Num)
  out1 <- EMGauss(cdr_data$Hyper.Num,K = 2,Means = c(0,100),SDs = c(20,20),MaxNumberofIterations =20)
  thred1 <- quantile(cdr_data$Hyper.Num,out1$Weights[1])
  out2 <- EMGauss(cdr_data$Hypo.Num,K = 2,Means = c(0,100),SDs = c(20,20),MaxNumberofIterations =20)
  thred2 <- quantile(cdr_data$Hypo.Num,out2$Weights[1])
  #km_mat <- cbind(cdr_data$Hyper.Num,cdr_data$Hypo.Num*max(cdr_data$Hyper.Num)/max(cdr_data$Hypo.Num))
  #kmout <- kmeans(km_mat,centers = 2,iter.max = 20)
  #h1 <- mean(cdr_data$Hyper.Num[kmout$cluster == 1])
  #h2 <- mean(cdr_data$Hyper.Num[kmout$cluster == 2])
  lmodel <- lm(Hyper.Num ~ Hypo.Num,data = cdr_data)

  cdr_data$me_label <- "Cluster1"
  if (cancerTypes3[idx] == "LIHC"){
    cdr_data$me_label[cdr_data$Hypo.Num > thred2] <- "Cluster2"
  } else {
    cdr_data$me_label[cdr_data$Hyper.Num > thred1] <- "Cluster2"
  }
  
  # for (nn in 1:nrow(cdr_data)){
  #   if (cdr_data$Hyper.Num[nn] >= lmodel$coefficients[1] + lmodel$coefficients[2]*cdr_data$Hypo.Num[nn]){
  #     cdr_data$me_label[nn] <- "Cluster2"
  #   }
  # }
  
  cdr_data$me_label <- factor(cdr_data$me_label,levels = c("Cluster1","Cluster2"))
  write.table(cdr_data,paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Survival_purity/typesDichotomy_",cancerTypes3[idx],".txt"),
              sep = "\t",quote = F)
  
  #os
  cdr_data_os <- cdr_data[(is.na(cdr_data$OS.time) == F),]#& (cdr_data$label == "Included")
  tnum <- sort(table(cdr_data_os$me_label),decreasing = T)
  if ((tnum[1]>=5) & (tnum[2]>=5)) {
    coxfit <- coxph(Surv(OS.time,OS) ~ me_label, data=cdr_data_os)
    x.out <- summary(coxfit)
    tmp_out <- c(tmp_out,x.out$conf.int[1],x.out$sctest[3])
    
    fit <- survfit(Surv(OS.time,OS) ~ me_label, data=cdr_data_os)
    surv_pvalue(fit)
    p <- ggsurvplot(fit,
                    pval = TRUE, conf.int = F,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
                    
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    palette = c("#0076BB", "#E69F00", "#56B4E9", "#800000"),
                    pval.method = T,title = paste0(cancerTypes3[idx],"-","OS"))
    print(p)
  } else {
    tmp_out <- c(tmp_out,NA,NA)
  }
  #dss
  cdr_data_dss <- cdr_data[(is.na(cdr_data$DSS.time) == F),] #& (cdr_data$label == "Included")
  tnum <- sort(table(cdr_data_dss$me_label),decreasing = T)
  if ((tnum[1]>=5) & (tnum[2]>=5)) {
    coxfit <- coxph(Surv(DSS.time,DSS) ~ me_label, data=cdr_data_dss)
    x.out <- summary(coxfit)
    tmp_out <- c(tmp_out,x.out$conf.int[1],x.out$sctest[3])
    
    fit <- survfit(Surv(DSS.time,DSS) ~ me_label, data=cdr_data_dss)
    p <- ggsurvplot(fit,
                    pval = TRUE, conf.int = F,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
                    
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    palette = c("#0076BB", "#E69F00", "#56B4E9", "#800000"),
                    pval.method = T,title = paste0(cancerTypes3[idx],"-","DSS"))
    print(p)
  } else {
    tmp_out <- c(tmp_out, NA,NA)
  }
  #dfi
  cdr_data_dfi <- cdr_data[(is.na(cdr_data$DFI.time) == F) ,] # & (cdr_data$label == "Included")
  tnum <- sort(table(cdr_data_dfi$me_label),decreasing = T)
  if ((tnum[1]>=5) & (tnum[2]>=5)) {
    coxfit <- coxph(Surv(DFI.time,DFI) ~ me_label, data=cdr_data_dfi)
    x.out <- summary(coxfit)
    tmp_out <- c(tmp_out,x.out$conf.int[1],x.out$sctest[3])
    
    fit <- survfit(Surv(DFI.time,DFI) ~ me_label, data=cdr_data_dfi)
    p <- ggsurvplot(fit,
                    pval = TRUE, conf.int = F,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
                    
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    palette = c("#0076BB", "#E69F00", "#56B4E9", "#800000"),
                    pval.method = T,title = paste0(cancerTypes3[idx],"-","DFI"))
    print(p)
  }  else {
    tmp_out <- c(tmp_out, NA,NA)
  }
  #pfi
  cdr_data_pfi <- cdr_data[(is.na(cdr_data$PFI.time) == F) ,] #& (cdr_data$label == "Included")
  tnum <- sort(table(cdr_data_pfi$me_label),decreasing = T)
  if ((tnum[1]>=5) & (tnum[2]>=5)) {
    coxfit <- coxph(Surv(PFI.time,PFI) ~ me_label, data=cdr_data_pfi)
    x.out <- summary(coxfit)
    tmp_out <- c(tmp_out,x.out$conf.int[1], x.out$sctest[3])
    
    fit <- survfit(Surv(PFI.time,PFI) ~ me_label, data=cdr_data_pfi)
    p <- ggsurvplot(fit,
                    pval = TRUE, conf.int = F,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
                    
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    palette = c("#0076BB", "#E69F00", "#56B4E9", "#800000"),
                    pval.method = T,title = paste0(cancerTypes3[idx],"-","PFI"))
    print(p)
  }  else {
    tmp_out <- c(tmp_out,NA,NA)
  }
  pancancer_out <- rbind(pancancer_out,c(tmp_out,max(nrow(cdr_data_os),
                                                     nrow(cdr_data_dss),
                                                     nrow(cdr_data_dfi),
                                                     nrow(cdr_data_pfi))))
}

dev.off()
rownames(pancancer_out) <- c(cancerTypes1,cancerTypes2,cancerTypes3)
colnames(pancancer_out) <- c("OS.HR","OS.Pvalue",
                             "DSS.HR","DSS.Pvalue",
                             "DFI.HR","DFI.Pvalue",
                             "PFI.HR","PFI.Pvalue",
                             "SampleNumInSurvival")
pancancer_out <- as.data.frame(pancancer_out)
write.table(pancancer_out, "Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Survival_purity/typesPvalue_table.txt",
            sep = "\t",quote = F)



## distributions
pdf("I:/Projects/P7.Methy-PanCancer/pics/tcga_survival_purity/typesdistribution.pdf",height = 4,width = 6)
for (idx in 1: length(cancerTypes1)){
  disdata <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Survival_purity/typesDichotomy_",cancerTypes1[idx],".txt"),
                         sep = "\t",header = T)
  disdata$me_label <- factor(disdata$me_label,levels = c("LowEpiMut","HighEpiMut"))
  g <- ggplot(disdata, aes(x = Hypo.Num,y = Hyper.Num, fill = me_label)) + geom_point(shape = 21,colour = "black",size = 5) + 
    theme_bw() + scale_fill_manual(values=c("#0076BB", "#E69F00", "#56B4E9", "#800000")) + ggtitle(cancerTypes1[idx])
  print(g)
}

for (idx in 1: length(cancerTypes2)){
  disdata <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Survival_purity/typesDichotomy_",cancerTypes2[idx],".txt"),
                        sep = "\t",header = T)
  disdata$me_label <- factor(disdata$me_label,levels = c("LowEpiMut","HighEpiMut"))
  g <- ggplot(disdata, aes(x = Hypo.Num,y = Hyper.Num, fill = me_label)) + geom_point(shape = 21,colour = "black",size = 5) + 
    theme_bw() + scale_fill_manual(values=c("#0076BB", "#E69F00", "#56B4E9", "#800000")) + ggtitle(cancerTypes2[idx])
  print(g)
}

for (idx in 1: length(cancerTypes3)){
  disdata <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Survival_purity/typesDichotomy_",cancerTypes3[idx],".txt"),
                        sep = "\t",header = T)
  disdata$me_label <- factor(disdata$me_label,levels = c("Cluster1","Cluster2"))
  g <- ggplot(disdata, aes(x = Hypo.Num,y = Hyper.Num, fill = me_label)) + geom_point(shape = 21,colour = "black",size = 5) + 
    theme_bw() + scale_fill_manual(values=c("#0076BB", "#E69F00", "#56B4E9", "#800000")) + ggtitle(cancerTypes3[idx])
  print(g)
}

dev.off()








  
  #out1 <- EMGauss(cdr_data$Hyper.Num[cdr_data$label == "Included"],K = 2,Means = c(0,100),SDs = c(20,20),MaxNumberofIterations =20)
  #thred1 <- quantile(cdr_data$Hyper.Num,out1$Weights[1])
 





