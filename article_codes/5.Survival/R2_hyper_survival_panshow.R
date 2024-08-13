surv <- read.table("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Survival_hyper/Hyper_survivalHRPvalue_table.txt",
                   sep = "\t",header = T)

lvs <- c("CESC","LGG","COAD","CHOL","ESCA","HNSC","READ","UCEC","GBM","BRCA",
         "PRAD","BLCA","LUAD","LUSC","LIHC","LAML","PAAD","SKCM","ACC","STAD",
         "SARC","KIRC","KIRP","PCPG","KICH","THCA","THYM")
show  <- NULL
for (idx in 1:nrow(surv)){
  if (is.na(surv[idx,1]) == F){
    show <- rbind(show, c(rownames(surv)[idx],"OS",surv[idx,1],surv[idx,2]))
  }
  if (is.na(surv[idx,3]) == F){
    show <- rbind(show, c(rownames(surv)[idx],"DSS",surv[idx,3],surv[idx,4]))
  }
  if (is.na(surv[idx,5]) == F){
    show <- rbind(show, c(rownames(surv)[idx],"DFI",surv[idx,5],surv[idx,6]))
  }
  if (is.na(surv[idx,7]) == F){
    show <- rbind(show, c(rownames(surv)[idx],"PFI",surv[idx,7],surv[idx,8]))
  }
}

show <- as.data.frame(show)
colnames(show) <- c("Cancer","Type","HR","Pvalue")
show$Cancer <- factor(show$Cancer,levels = lvs)
show$Type <- factor(show$Type,levels = rev(c("OS","DSS","DFI","PFI")))
show$HR <-as.numeric(show$HR)
show$Pvalue <-as.numeric(show$Pvalue)
show$log2HR <- log2(show$HR)

library(ggplot2)
library(scales)
library(paletteer)
show$sig <- sign(show$log2HR)*(-log10(show$Pvalue))
show$sig[show$sig>6] = 6
show$sig[show$sig< -6] = -6
show$label <- 0
show$label[abs(show$sig) >= abs(log10(0.05))] <- 1
show2 <- show[show$label == 1,]

pdf("I:/Projects/P7.Methy-PanCancer/pics/TCGA_hypersurvival_integrated_points.pdf",height = 2.5,width = 15)
g <- ggplot(show,aes(y=Type,x=Cancer,fill = sig,size = 2*abs(log2HR))) + geom_point(shape = 21)+ 
  geom_point(mapping = aes(y = Type, x= Cancer),data = show2,shape = 3,size = 2, show.legend = F )+
  scale_color_gradientn(colors = paletteer_d("dichromat::BluetoOrange_12")  ,limits = c(-6,6),aesthetics = "fill") + theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(g)
dev.off()

