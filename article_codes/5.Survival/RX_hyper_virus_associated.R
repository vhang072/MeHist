#cancerTypes <- c("ESCA","GBM","LAML","LGG","LUSC","STAD","UCEC")
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

viral <- read.table("Y:/4.basic_data/TCGA_PancanAtlas/Clinical_TCGAperType/viral_associated/pancancer_viralExpression.txt",
                    sep = "\t",header = T)
rownames(viral) <- str_replace_all(viral$TCGA_patientID,"-",".")
viral_f <- viral[,57:ncol(viral)]

out_vis <- NULL
for (idx in 1:length(cancerTypes)){
  subclu <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Survival_hyper/HyperDichotomy_",cancerTypes[idx],".txt"),
                       sep = "\t",header = T)
  overlaps <- intersect(rownames(subclu),rownames(viral_f))
  cls <- subclu[overlaps,]
  vis <- viral_f[overlaps,]
  if (length(overlaps)>10){
    for (x in 1:ncol(vis)){
      x1 <- ifelse(vis[,x]>1,1,0)
      x2 <- cls[,"hyper_label"]
      if (sum(x1)>=3){
        t <- table(factor(x1),x2)
        test1 <- fisher.test(t)
        out_vis <- rbind(out_vis,c(cancerTypes[idx],colnames(vis)[x],test1$p.value,t[1,1],t[1,2],t[2,1],t[2,2]))
      }
    }
  }
}  
out_vis <- as.data.frame(out_vis)
  colnames(out_vis) <- c("Cancer","Virus","Pvalue","T1.1","T1.2","T2.1","T2.2")
  
  write.table(out_vis,"Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/virus_hyper/virus_hyper.txt",sep = "\t",
              row.names = F,quote = F)

  
  ##draw
  ct<- "STAD"
  gz <- "EBV"
  subclu <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Survival_hyper/HyperDichotomy_",ct,".txt"),
                       sep = "\t",header = T)
  viral <- read.table("Y:/4.basic_data/TCGA_PancanAtlas/Clinical_TCGAperType/viral_associated/pancancer_viralExpression.txt",
                      sep = "\t",header = T)
  rownames(viral) <- str_replace_all(viral$TCGA_patientID,"-",".")
  viral_f <- viral[,57:ncol(viral)]
  
  
  overlaps <- intersect(rownames(subclu),rownames(viral_f))
  cluf <- subclu[overlaps,]
  virf <- viral_f[overlaps,gz]
  virf[virf>=1] <- 1
  virf[virf< 1] <- 0
  cluf$vir <- factor(as.numeric(virf))
  
  library(ggplot2)
  pdf(paste0("I:/Projects/P7.Methy-PanCancer/pics/mutation_virus/",ct,"-",gz,".pdf"),height = 4,width = 5)
  g <- ggplot(cluf,aes(x = hypo_num,y = hyper_num,fill = vir)) + geom_point(shape = 21,color = "black",size = 2) +
    theme_bw() + scale_fill_manual(values=c("#999999", "#800000")) + ggtitle(paste0(ct,"-",gz))
  print(g)
  dev.off()
  