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

cor_hypo <- NULL
pval_hypo <- NULL
for (x in 1:length(cancerTypes)) {
  data1 <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_TCGA/_DiffCGIprobeMean_",cancerTypes[x],".txt"),
                      sep = "\t",header = T)
  data2 <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_L1/Diff_TCGA/_DiffL1probeMean_",cancerTypes[x],".txt"),
                      sep = "\t",header = T)
  data3 <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_soloWCGW/Diff_TCGA/_DiffsoloWCGWprobeMean_",cancerTypes[x],".txt"),
                      sep = "\t",header = T)
  data4 <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_opensea/Diff_TCGA/_DiffopenseaprobeMean_",cancerTypes[x],".txt"),
                      sep = "\t",header = T)
  samples <- intersect(intersect(intersect(colnames(data1)[2:ncol(data1)],colnames(data2)[2:ncol(data2)]),
                      colnames(data3)[2:ncol(data3)]),
            colnames(data4)[2:ncol(data4)])
  y <- rbind(data1[2,samples],data2[2,samples],data4[2,samples],data3[2,samples])

  tmp.cor <- NULL
  tmp.pval <- NULL
  for (idx in 2:nrow(y)){
    corr1 <- cor.test(t(y[1,]),t(y[idx,]))
    tmp.cor <- c(tmp.cor,-1*corr1$estimate)
    tmp.pval <- c(tmp.pval,corr1$p.value)
  }
  cor_hypo <- rbind(cor_hypo,tmp.cor)
  pval_hypo <- rbind(pval_hypo,tmp.pval)
}

items <- c("L1","OpenSea","SoloWCGW")

rownames(cor_hypo) <- cancerTypes
colnames(cor_hypo) <- paste0("Cor.",items)

rownames(pval_hypo) <- cancerTypes
colnames(pval_hypo) <- paste0("P.",items)

tmp.order <- order(as.numeric(cor_hypo[,1]),decreasing = T)
order.cancer <- rownames(cor_hypo)[tmp.order]

show <- NULL
for (x in 1:nrow(cor_hypo)){
  for (y in 1:ncol(cor_hypo)){
    show <- rbind(show,c(rownames(cor_hypo)[x],colnames(cor_hypo)[y],cor_hypo[x,y],pval_hypo[x,y]))
  }
}
show <- as.data.frame(show)
colnames(show) <- c("Cancer","Cor.Item","Cor.Pearson","Pvalue")
show$Cor.Pearson <- as.numeric(show$Cor.Pearson)
show$Pvalue <- as.numeric(show$Pvalue)
show$Cancer <- factor(show$Cancer,levels = order.cancer)
show$Cor.Item <- factor(show$Cor.Item,levels = paste0("Cor.",items))

show$sig <- sign(show$Cor.Pearson)*(-log10(show$Pvalue))
show$sig[show$sig>30] = 30
show$sig[show$sig< -30] = -30
show$Cor.Item <- factor(show$Cor.Item,levels = c("Cor.SoloWCGW","Cor.OpenSea","Cor.L1"))

library(ggplot2)
library(scales)
pdf("I:/Projects/P7.Methy-PanCancer/pics/hypo_cors/TCGA_hypoCor.pdf",height = 2.1,width = 9)
g <- ggplot(show,aes(x=Cancer,y=Cor.Item,fill = sig,size = abs(Cor.Pearson))) + geom_point(shape = 21)+
  scale_color_gradientn(colors = paletteer_d("dichromat::DarkRedtoBlue_12")  ,limits = c(-30,30),aesthetics = "fill") +
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_size_continuous(limits = c(0,1),range = c(0, 6))
print(g)
dev.off()
