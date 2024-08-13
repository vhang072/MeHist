cancerTypes <- c(
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
                 "LAML",
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
library(ggplot2)


pdf("I:/Projects/P7.Methy-PanCancer/pics/hyeprhypodichotomyShow.pdf",height = 4,width = 6)
for (idx in 1: length(cancerTypes)){
  x1 <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Survival_hyper/HyperDichotomy_",cancerTypes[idx],".txt"),
                   sep = "\t",header = T)
  x2 <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Survival_hypo/HyperDichotomy_",cancerTypes[idx],".txt"),
                   sep = "\t",header = T)  
  catval <- data.frame(idz = rownames(x1),hyper.num = x1$hyper_num,hypo.num = x1$hypo_num)
  catval$label <- "Control"
  catval$label[(x1$hyper_label == "Highly_Hyper") & (x2$hypo_label == "Control")] <- "Hyper"
  catval$label[(x1$hyper_label == "Control") & (x2$hypo_label == "Highly_Hypo")] <- "Hypo"
  catval$label[(x1$hyper_label == "Highly_Hyper") & (x2$hypo_label == "Highly_Hypo")] <- "Intermediate"
  catval$label <- factor(catval$label,levels = c("Control","Hyper","Hypo","Intermediate"))
  g <- ggplot(catval, aes(x = hypo.num,y = hyper.num, fill = label)) + geom_point(shape = 21,colour = "black",size = 3) + 
    theme_bw() + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9", "#800000")) + ggtitle(cancerTypes[idx])
  print(g)
}

dev.off()



####survivals
pancancer_out <- NULL
library(survival)
library(survminer)
library(gtsummary)
pdf("I:/Projects/P7.Methy-PanCancer/pics/hyeprhypotogetherSurvival.pdf",height = 6,width = 6)
for (idx in 1: length(cancerTypes)){
  tmp_out <- NULL
  x1 <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Survival_hyper/HyperDichotomy_",cancerTypes[idx],".txt"),
                   sep = "\t",header = T)
  x2 <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Survival_hypo/HyperDichotomy_",cancerTypes[idx],".txt"),
                   sep = "\t",header = T)  
  catval <- data.frame(idz = rownames(x1),hyper.num = x1$hyper_num,hypo.num = x1$hypo_num)
  catval$label <- "Control"
  catval$label[(x1$hyper_label == "Highly_Hyper") & (x2$hypo_label == "Control")] <- "Hyper"
  catval$label[(x1$hyper_label == "Control") & (x2$hypo_label == "Highly_Hypo")] <- "Hypo"
  catval$label[(x1$hyper_label == "Highly_Hyper") & (x2$hypo_label == "Highly_Hypo")] <- "Intermediate"
  catval$label <- factor(catval$label,levels = c("Control","Hyper","Hypo","Intermediate"))
  
  toge <- cbind(x1[,1:8], catval)
  
  #os
  cdr_data_os <- toge[is.na(toge$OS.time) == F,]
  ns <- c(sum(cdr_data_os$label ==  "Hyper"), sum(cdr_data_os$label ==  "Hypo"),
          sum(cdr_data_os$label ==  "Intermediate"),sum(cdr_data_os$label ==  "Control"))
  if ((sort(ns,decreasing = T)[1] >=5) & (sort(ns,decreasing = T)[2] >=5)) {
    coxfit <- coxph(Surv(OS.time,OS) ~ label, data=cdr_data_os)
    x.out <- summary(coxfit)
    tmp_out <- c(tmp_out,x.out$conf.int[1],x.out$sctest[3])
    
    fit <- survfit(Surv(OS.time,OS) ~ label, data=cdr_data_os)
    p <- ggsurvplot(fit,
                    pval = TRUE, conf.int = F,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
                    
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    palette = c("#999999", "#E69F00", "#56B4E9", "#800000"),
                    pval.method = T,title = paste0(cancerTypes[idx],"-","OS"))
    print(p)
  } else {
    tmp_out <- c(tmp_out, NA, NA)
  }
  
  #dss
  cdr_data_dss <- toge[is.na(toge$DSS.time) == F,]
  ns <- c(sum(cdr_data_dss$label ==  "Hyper"), sum(cdr_data_dss$label ==  "Hypo"),
          sum(cdr_data_dss$label ==  "Intermediate"),sum(cdr_data_dss$label ==  "Control"))
  if ((sort(ns,decreasing = T)[1] >=5) & (sort(ns,decreasing = T)[2] >=5)) {
    coxfit <- coxph(Surv(DSS.time,DSS) ~ label, data=cdr_data_dss)
    x.out <- summary(coxfit)
    tmp_out <- c(tmp_out,x.out$conf.int[1],x.out$sctest[3])
    
    fit <- survfit(Surv(DSS.time,DSS) ~ label, data=cdr_data_dss)
    p <- ggsurvplot(fit,
                    pval = TRUE, conf.int = F,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
                    
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    palette = c( "#999999", "#E69F00", "#56B4E9", "#800000"),
                    pval.method = T,title = paste0(cancerTypes[idx],"-","DSS"))
    print(p)
  } else {
    tmp_out <- c(tmp_out, NA, NA)
  }  
  
  #dfi
  cdr_data_dfi <- toge[is.na(toge$DFI.time) == F,]
  ns <- c(sum(cdr_data_dfi$label ==  "Hyper"), sum(cdr_data_dfi$label ==  "Hypo"),
          sum(cdr_data_dfi$label ==  "Intermediate"),sum(cdr_data_dfi$label ==  "Control"))
  if ((sort(ns,decreasing = T)[1] >=5) & (sort(ns,decreasing = T)[2] >=5)) {
    coxfit <- coxph(Surv(DFI.time,DFI) ~ label, data=cdr_data_dfi)
    x.out <- summary(coxfit)
    tmp_out <- c(tmp_out,x.out$conf.int[1],x.out$sctest[3])
    
    fit <- survfit(Surv(DFI.time,DFI) ~ label, data=cdr_data_dfi)
    p <- ggsurvplot(fit,
                    pval = TRUE, conf.int = F,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
                    
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    palette = c( "#999999", "#E69F00", "#56B4E9", "#800000"),
                    pval.method = T,title = paste0(cancerTypes[idx],"-","DFI"))
    print(p)
  } else {
    tmp_out <- c(tmp_out, NA, NA)
  }    
  
  #pfi
  cdr_data_pfi <- toge[is.na(toge$PFI.time) == F,]
  ns <- c(sum(cdr_data_pfi$label ==  "Hyper"), sum(cdr_data_pfi$label ==  "Hypo"),
          sum(cdr_data_pfi$label ==  "Intermediate"),sum(cdr_data_pfi$label ==  "Control"))
  if ((sort(ns,decreasing = T)[1] >=5) & (sort(ns,decreasing = T)[2] >=5)) {
    coxfit <- coxph(Surv(PFI.time,PFI) ~ label, data=cdr_data_pfi)
    x.out <- summary(coxfit)
    
    tmp_out <- c(tmp_out,x.out$conf.int[1],x.out$sctest[3])
    
    fit <- survfit(Surv(PFI.time,PFI) ~ label, data=cdr_data_pfi)
    p <- ggsurvplot(fit,
                    pval = TRUE, conf.int = F,
                    risk.table = TRUE, # Add risk table
                    risk.table.col = "strata", # Change risk table color by groups
                    linetype = "solid", # Change line type by groups surv.median.line = "hv", # Specify median survival
                    
                    ggtheme = theme_bw(), # Change ggplot2 theme
                    palette = c( "#999999", "#E69F00", "#56B4E9", "#800000"),
                    pval.method = T,title = paste0(cancerTypes[idx],"-","PFI"))
    print(p)
  } else {
    tmp_out <- c(tmp_out, NA, NA)
  }  
  
  pancancer_out <- rbind(pancancer_out,tmp_out)
}
dev.off()

rownames(pancancer_out) <- cancerTypes
colnames(pancancer_out) <- c("OS.HR","OS.Pvalue","DSS.HR","DSS.Pvalue","DFI.HR","DFI.Pvalue","PFI.HR","PFI.Pvalue")
pancancer_out <- as.data.frame(pancancer_out)
write.table(pancancer_out, "Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Survival_together/HyperHypoInter_survivalHRPvalue_table.txt",
            sep = "\t",quote = F)



