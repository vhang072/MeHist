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


#pdf(file = "Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/_PDF_pic/normal_mean_hist_cellline.pdf",height = 5,width = 7)

path <- c("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_TCGA/_DiffCGIprobe_",
          "Y:/4.basic_data/TCGA_PancanAtlas/methylation_L1/Diff_TCGA/_DiffL1probe_",
          "Y:/4.basic_data/TCGA_PancanAtlas/methylation_soloWCGW/Diff_TCGA/_DiffsoloWCGWprobe_",
          "Y:/4.basic_data/TCGA_PancanAtlas/methylation_opensea/Diff_TCGA/_DiffOpenSeaprobe_")



diff_probes <- NULL
for (x in 1:length(cancerTypes)){
  tmp <- NULL
  for (y in 1:length(path)){
    data <- read.table(paste0(path[y],cancerTypes[x],".txt"),sep = "\t",
                       header = T)
    data$label <- "Control"
    data$label[data$Mean_T > data$Mean_N & data$AberrantFrac_T_0.2 > 0.25 & data$Mean_N < 0.25] <- "Global_Hyper"
    data$label[data$Mean_T < data$Mean_N & data$AberrantFrac_T_0.2 > 0.25 & data$Mean_N > 0.75] <- "Global_Hypo"
    
    tmp <- c(tmp,sum(data$label == "Global_Hyper"),
             sum(data$label == "Global_Hypo"),
             nrow(data))
  }
  diff_probes <- rbind(diff_probes,tmp)
}

diff_probes <- cbind(cancerTypes,diff_probes)
diff_probes <- as.data.frame(diff_probes)
colnames(diff_probes) <- c("cancerTypes","hyper_cgi","hypo_cgi","totalprobe_cgi",
                           "hyper_L1","hypo_L1","totalprobe_L1",
                           "hyper_soloWCGW","hypo_soloWCGW","totalprobe_soloWCGW",
                           "hyper_opensea","hypo_opensea","totalprobe_opensea")
for (idx in 2:ncol(diff_probes)){
  diff_probes[,idx] <- as.numeric(diff_probes[,idx])
}
lv <- diff_probes$cancerTypes[order(diff_probes$hyper_cgi,decreasing = T)]
diff_probes$cancerTypes <- factor(diff_probes$cancerTypes,levels = lv)
library(ggplot2)
pdf("Y:/4.basic_data/TCGA_PancanAtlas/methylation_integrated/TCGA/summary_cgidiff.pdf",height = 5,width = 8)
ggplot(diff_probes,aes(x = cancerTypes,y = hyper_cgi)) + geom_col(fill = "#FF8000")+
  geom_col(aes(x = cancerTypes,y = -1*hypo_cgi),fill = "#6A5ACD") + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

write.table(diff_probes,"Y:/4.basic_data/TCGA_PancanAtlas/methylation_integrated/TCGA/_DiffProbeNumSummary.txt",
            sep = "\t",row.names = F,col.names = T, quote = F)
