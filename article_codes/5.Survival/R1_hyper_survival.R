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

library(stringr)
library(AdaptGauss)
library(dplyr)


cdr <- read.table("Y:/4.basic_data/TCGA_PancanAtlas/clinical/TCGA-CDR.txt",sep = "\t",header = T)
rn <- NULL 
for (x in 1:nrow(cdr)){
  rn <- c(rn, str_replace_all(cdr$bcr_patient_barcode[x],'-',"."))
}
rownames(cdr) <- rn


merged_mepersample <- NULL
median_num <- NULL
pdf("I:/Projects/P7.Methy-PanCancer/pics/hypersurvival.pdf",height = 6,width = 6)
pancancer_out <- NULL
for (idx in 1:length(cancerTypes)){
  tmp_out <- NULL
  data0 <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_TCGA/_DiffCGIprobeMean_",cancerTypes[idx],".txt"),
                      sep = "\t",header = T)
  data <- data0[,2:ncol(data0)]
  tmp.merged <- data.frame(idz = colnames(data),mean_hyper = as.numeric(unname(data[1,])),mean_hypo = as.numeric(unname(data[2,])),
                           num_hyper = as.numeric(unname(data[3,])),num_hypo = as.numeric(unname(data[4,])),cancertype = rep(cancerTypes[idx],ncol(data)))
  merged_mepersample <- rbind(merged_mepersample,tmp.merged)
  median_num <- c(median_num,median(as.numeric(unname(data[3,]))))
  colnm <- NULL
  colx <- NULL
  for (x in 1:ncol(data)){
    tmp.name <- substr(colnames(data)[x],1,12)
    if (tmp.name %in% colnm == F){
      colnm <- c(colnm,tmp.name)
      colx <- c(colx,x)
    }
  }
  data_f <- data[,colx]
  colnames(data_f) <- colnm

  overlap_sps <- intersect(rownames(cdr),colnames(data_f))
  
  data_u <- data_f[,overlap_sps]
  rownames(data_u) <- c("hyper_mean","hypo_mean","hyper_num","hypo_num")
  cdr_u <- cdr[overlap_sps,]
  
  cdr_data <- cbind(cdr_u[,26:33],t(data_u))
  
  out <- EMGauss(cdr_data$hyper_num,K = 2,Means = c(0,100),SDs = c(20,20),MaxNumberofIterations =20)
  
  thred <- quantile(cdr_data$hyper_num,out$Weights[1])
  
  cdr_data$hyper_label <- factor(ifelse(cdr_data$hyper_num > thred, "Highly_Hyper","Control"),
                                     levels = c("Control","Highly_Hyper"))
  write.table(cdr_data,paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Survival_hyper/HyperDichotomy_",cancerTypes[idx],".txt"),
              sep = "\t",quote = F)
  
  
  
  library(survival)
  library(survminer)
  library(gtsummary)
  #os
  cdr_data_os <- cdr_data[is.na(cdr_data$OS.time) == F,]
  if ((sum(cdr_data_os$hyper_label ==  "Highly_Hyper")>=5) & (sum(cdr_data_os$hyper_label ==  "Control")>=5)) {
    coxfit <- coxph(Surv(OS.time,OS) ~ hyper_label, data=cdr_data_os)
    x.out <- summary(coxfit)
    tmp_out <- c(tmp_out,x.out$conf.int[1],x.out$sctest[3])
    
    fit <- survfit(Surv(OS.time,OS) ~ hyper_label, data=cdr_data_os)
    surv_pvalue(fit)
    p <- ggsurvplot(fit,
               pval = TRUE, conf.int = F,
               risk.table = TRUE, # Add risk table
               risk.table.col = "strata", # Change risk table color by groups
               linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
               
               ggtheme = theme_bw(), # Change ggplot2 theme
               palette = c( "#2E9FDF","#E7B800"),
               pval.method = T,title = paste0(cancerTypes[idx],"-","OS"))
    print(p)
  } else {
    tmp_out <- c(tmp_out,NA,NA)
  }
  #dss
  cdr_data_dss <- cdr_data[is.na(cdr_data$DSS.time) == F,]
  if ((sum(cdr_data_dss$hyper_label ==  "Highly_Hyper")>=5) & (sum(cdr_data_dss$hyper_label ==  "Control")>=5)) {
    coxfit <- coxph(Surv(DSS.time,DSS) ~ hyper_label, data=cdr_data_dss)
    x.out <- summary(coxfit)
    tmp_out <- c(tmp_out,x.out$conf.int[1],x.out$sctest[3])
    
    fit <- survfit(Surv(DSS.time,DSS) ~ hyper_label, data=cdr_data_dss)
    p <- ggsurvplot(fit,
                    pval = TRUE, conf.int = F,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
                    
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    palette = c( "#2E9FDF","#E7B800"),
                    pval.method = T,title = paste0(cancerTypes[idx],"-","DSS"))
    print(p)
  } else {
    tmp_out <- c(tmp_out, NA,NA)
  }
  #dfi
  cdr_data_dfi <- cdr_data[is.na(cdr_data$DFI.time) == F,]
  if ((sum(cdr_data_dfi$hyper_label ==  "Highly_Hyper")>=5) & (sum(cdr_data_dfi$hyper_label ==  "Control")>=5)) {
    coxfit <- coxph(Surv(DFI.time,DFI) ~ hyper_label, data=cdr_data_dfi)
    x.out <- summary(coxfit)
    tmp_out <- c(tmp_out,x.out$conf.int[1],x.out$sctest[3])
    
    fit <- survfit(Surv(DFI.time,DFI) ~ hyper_label, data=cdr_data_dfi)
    p <- ggsurvplot(fit,
                    pval = TRUE, conf.int = F,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
                    
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    palette = c( "#2E9FDF","#E7B800"),
                    pval.method = T,title = paste0(cancerTypes[idx],"-","DFI"))
    print(p)
  }  else {
    tmp_out <- c(tmp_out, NA, NA)
  }
  #pfi
  cdr_data_pfi <- cdr_data[is.na(cdr_data$PFI.time) == F,]
  if ((sum(cdr_data_pfi$hyper_label ==  "Highly_Hyper")>=5) & (sum(cdr_data_pfi$hyper_label ==  "Control")>=5)) {
    coxfit <- coxph(Surv(PFI.time,PFI) ~ hyper_label, data=cdr_data_pfi)
    x.out <- summary(coxfit)
    tmp_out <- c(tmp_out,x.out$conf.int[1],x.out$sctest[3])
    
    fit <- survfit(Surv(PFI.time,PFI) ~ hyper_label, data=cdr_data_pfi)
    p <- ggsurvplot(fit,
                    pval = TRUE, conf.int = F,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
                    
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    palette = c( "#2E9FDF","#E7B800"),
                    pval.method = T,title = paste0(cancerTypes[idx],"-","PFI"))
    print(p)
  }  else {
    tmp_out <- c(tmp_out,NA,NA)
  }
  pancancer_out <- rbind(pancancer_out,tmp_out)
}
dev.off()
rownames(pancancer_out) <- cancerTypes
colnames(pancancer_out) <- c("OS.HR","OS.Pvalue","DSS.HR","DSS.Pvalue","DFI.HR","DFI.Pvalue","PFI.HR","PFI.Pvalue")
pancancer_out <- as.data.frame(pancancer_out)
write.table(pancancer_out, "Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Survival_hyper/Hyper_survivalHRPvalue_table.txt",
            sep = "\t",quote = F)


## distributions

lv <- cancerTypes[order(median_num,decreasing = T)]
merged_mepersample$cancertype <- factor(merged_mepersample$cancertype,levels = lv)

library(ggplot2)
pdf("I:/Projects/P7.Methy-PanCancer/pics/hyper_hypo_num.pdf",height = 4,width = 16)
p <- ggplot(merged_mepersample,aes(x = cancertype,y = num_hyper)) + 
  geom_violin(draw_quantiles = 0.5,linewidth = 0.1,width = 0.8,color = "black",scale = "width",adjust = 0.8) + 
  geom_dotplot(color = "#ff8000",binaxis='y', stackdir='center', dotsize=0.2,binwidth = 200) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
print(p)

p <- ggplot(merged_mepersample,aes(x = cancertype,y = num_hypo)) + 
  geom_violin(draw_quantiles = 0.5,linewidth = 0.1,width = 0.8,color = "black",scale = "width",adjust = 0.8) + 
  geom_dotplot(color = "#6a5acd",binaxis='y', stackdir='center', dotsize=0.2,binwidth = 50) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1))
print(p)
dev.off()



###pan-cancer survival
rownm <- NULL
rowx <- NULL
for (x in 1:nrow(merged_mepersample)){
  tmp.name <- substr(merged_mepersample$idz[x],1,12)
  if (tmp.name %in% rownm == F){
    rownm <- c(rownm,tmp.name)
    rowx <- c(rowx,x)
  }
}
merged_mepersample_f <- merged_mepersample[rowx,]
rownames(merged_mepersample_f) <- rownm

overlap_sps <- intersect(rownames(cdr),rownames(merged_mepersample_f))

merged_mepersample_u <- merged_mepersample_f[overlap_sps,]
cdr_u <- cdr[overlap_sps,]

cdr_data_pan <- cbind(cdr_u[,26:33],merged_mepersample_u)

pdf("I:/Projects/P7.Methy-PanCancer/pics/pancancer_survival.pdf",height = 6,width = 6)
hist(cdr_data_pan$num_hyper,breaks = 50)
hist(cdr_data_pan$num_hypo,breaks = 50)
out <- EMGauss(cdr_data_pan$num_hyper,K = 2,Means = c(0,100),SDs = c(20,20),MaxNumberofIterations =20)

thred <- quantile(cdr_data_pan$num_hyper,out$Weights[1])

cdr_data_pan$hyper_label <- factor(ifelse(cdr_data_pan$num_hyper > thred, "Highly_Hyper","Control"),
                               levels = c("Control","Highly_Hyper"))

out2 <- EMGauss(cdr_data_pan$num_hypo,K = 2,Means = c(0,100),SDs = c(20,20),MaxNumberofIterations =20)

thred2 <- quantile(cdr_data_pan$num_hypo,out2$Weights[1])

cdr_data_pan$hypo_label <- factor(ifelse(cdr_data_pan$num_hypo > thred2, "Highly_Hypo","Control"),
                                   levels = c("Control","Highly_Hypo"))
#os

cdr_data_os <- cdr_data_pan[is.na(cdr_data_pan$OS.time) == F,]
if ((sum(cdr_data_os$hyper_label ==  "Highly_Hyper")>=5) & (sum(cdr_data_os$hyper_label ==  "Control")>=5)) {
  fit <- survfit(Surv(OS.time,OS) ~ hyper_label, data=cdr_data_os)
  p <- ggsurvplot(fit,
                  pval = TRUE, conf.int = F,
                  risk.table = TRUE, # Add risk table
                  risk.table.col = "strata", # Change risk table color by groups
                  linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
                  
                  ggtheme = theme_bw(), # Change ggplot2 theme
                  palette = c( "#c0c0c0","#ff8000"),
                  pval.method = T,title = paste0("Pancancer-hyper","-","OS"),
                  xlim = c(0,6000))
  coxfit <- coxph(Surv(OS.time,OS) ~ hyper_label, data=cdr_data_os)
  x.out <- summary(coxfit)
  print(x.out$conf.int[1])
  print(x.out$sctest[3])
  
  print(p)
}

if ((sum(cdr_data_os$hypo_label ==  "Highly_Hypo")>=5) & (sum(cdr_data_os$hypo_label ==  "Control")>=5)) {
  fit <- survfit(Surv(OS.time,OS) ~ hypo_label, data=cdr_data_os)
  p <- ggsurvplot(fit,
                  pval = TRUE, conf.int = F,
                  risk.table = TRUE, # Add risk table
                  risk.table.col = "strata", # Change risk table color by groups
                  linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
                  
                  ggtheme = theme_bw(), # Change ggplot2 theme
                  palette = c( "#c0c0c0","#6a5acd"),
                  pval.method = T,title = paste0("Pancancer-hypo","-","OS"),
                  xlim = c(0,6000))
  coxfit <- coxph(Surv(OS.time,OS) ~ hypo_label, data=cdr_data_os)
  x.out <- summary(coxfit)
  print(x.out$conf.int[1])
  print(x.out$sctest[3])
  print(p)
}
dev.off()