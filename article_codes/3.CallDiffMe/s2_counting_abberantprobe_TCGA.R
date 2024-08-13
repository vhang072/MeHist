cancerTypes <- c("LAML",
                 "ACC",
                 "BLCA",
                 "BRCA",
                 "CESC",
                 "CHOL",
                 "COAD",
                 "DLBC",
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
                 "OV",
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
probelabel <- NULL

pdf(file = "Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/_PDF_pic/normal_mean_hist.pdf",height = 5,width = 7)

for (i in 1:length(cancerTypes)){
  data <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_TCGA/_DiffCGIprobe_",cancerTypes[i],".txt"),sep = "\t",
                     header = T)
  data$label <- "Control"
  data$label[data$Mean_T > data$Mean_N & data$AberrantFrac_T_0.1 > 0.5] <- "Global_Hyper"
  data$label[data$Mean_T < data$Mean_N & data$AberrantFrac_T_0.1 > 0.5] <- "Global_Hypo"
  probelabel <- cbind(probelabel,data$label)
  p1 <- hist(data$Mean_N[data$label == "Global_Hyper"],breaks = seq(0,1,0.05),main = paste0(cancerTypes[i],"_Global_Hyper"))
  print(p1)
  p2 <- hist(data$Mean_N[data$label == "Global_Hypo"],breaks = seq(0,1,0.05),main = paste0(cancerTypes[i],"_Global_Hypo"))
  print(p2) 
  hyper <- data[data$Mean_T > data$Mean_N,]
  hypo <- data[data$Mean_T < data$Mean_N,]
  p3 <- hist(hyper$AberrantFrac_T_0.1,breaks = seq(0,1,0.05),main = paste0(cancerTypes[i],"_Hyper_Frc"))
  print(p3)
  p4 <- hist(hypo$AberrantFrac_T_0.1,breaks = seq(0,1,0.05),main = paste0(cancerTypes[i],"_Hypo_Frac"))
  print(p4) 
}
dev.off()

probelabel <- as.data.frame(probelabel)
colnames(probelabel) <- cancerTypes
probelabel <- cbind(data[,1:4],probelabel)



diffnum <- NULL
for (i in 5:ncol(probelabel)){
  diffnum <- rbind(diffnum,c(sum(probelabel[,i] == "Global_Hyper"),sum(probelabel[,i] == "Global_Hypo")))
}
rownames(diffnum) <- cancerTypes
colnames(diffnum) <- c("hyper","hypo")
diffnum <- as.data.frame(diffnum)
diffnum$hyperfrac <- diffnum$hyper/max(diffnum$hyper)
diffnum$hypofrac <- diffnum$hypo/max(diffnum$hypo)
diffnum$total <- diffnum$hyperfrac + diffnum$hypofrac


plot(sort(diffnum$total))
