# cancerTypes <- c("aero_digestive_tract_head and neck",
#                  "aero_digestive_tract_oesophagus",
#                  "blood_B_cell_lymphoma",
#                  "blood_acute_myeloid_leukaemia",
#                  "blood_B_cell_leukemia",
#                  "bone_ewings_sarcoma",
#                  "bone_osteosarcoma",
#                  "breast_breast",
#                  "digestive_system_large_intestine",
#                  "digestive_system_stomach",
#                  "digestive_system_liver",
#                  "kidney_kidney",
#                  "lung_lung_NSCLC_adenocarcinoma",
#                  "lung_lung_small_cell_carcinoma",
#                  "lung_mesothelioma",
#                  "lung_lung_NSCLC_squamous_cell_carcinoma",
#                  "nervous_system_glioma",
#                  "nervous_system_neuroblastoma",
#                  "pancreas_pancreas",
#                  "skin_melanoma",
#                  "thyroid_thyroid",
#                  "urogenital_system_ovary",
#                  "urogenital_system_bladder",
#                  "urogenital_system_cervix",
#                  "urogenital_system_endometrium",
#                  "urogenital_system_prostate")


#more than 20 cell lines for one cancer type
cancerTypes <- c("aero_digestive_tract_head and neck",
                 "aero_digestive_tract_oesophagus",
                 "blood_B_cell_lymphoma",
                 "blood_acute_myeloid_leukaemia",
                 "bone_ewings_sarcoma",
                 "breast_breast",
                 "digestive_system_large_intestine",
                 "digestive_system_stomach",
                 "kidney_kidney",
                 "lung_lung_NSCLC_adenocarcinoma",
                 "lung_lung_small_cell_carcinoma",
                 "lung_mesothelioma",
                 "nervous_system_glioma",
                 "nervous_system_neuroblastoma",
                 "pancreas_pancreas",
                 "skin_melanoma",
                 "urogenital_system_ovary",
                 "urogenital_system_bladder")

cor_hypo <- NULL
pval_hypo <- NULL
for (x in 1:length(cancerTypes)) {
  data1 <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_CCLE/_DiffCGIprobeMean_",cancerTypes[x],".txt"),
                      sep = "\t",header = T)
  data2 <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_L1/Diff_CCLE/_DiffL1probeMean_",cancerTypes[x],".txt"),
                      sep = "\t",header = T)
  data3 <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_soloWCGW/Diff_CCLE/_DiffPMDsoloWCGWProbeMean_",cancerTypes[x],".txt"),
                      sep = "\t",header = T)
  data4 <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_opensea/Diff_CCLE/_DiffOpenseaProbeMean_",cancerTypes[x],".txt"),
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
show$Pvalue <- as.numeric(show$Pvalue) + 1e-50
show$Cancer <- factor(show$Cancer,levels = order.cancer)
show$Cor.Item <- factor(show$Cor.Item,levels = paste0("Cor.",items))

show$sig <- sign(show$Cor.Pearson)*(-log10(show$Pvalue))
show$sig[show$sig>20] = 20
show$sig[show$sig< -20] = -20
show$Cor.Item <- factor(show$Cor.Item,levels = c("Cor.SoloWCGW","Cor.OpenSea","Cor.L1"))

library(ggplot2)
library(scales)
pdf("I:/Projects/P7.Methy-PanCancer/pics/hypo_cors/CCLE_hypoCor.pdf",height = 3,width = 7)
g <- ggplot(show,aes(x=Cancer,y=Cor.Item,fill = sig,size = abs(Cor.Pearson))) + geom_point(shape = 21)+
  scale_color_gradientn(colors = paletteer_d("dichromat::DarkRedtoBlue_12")  ,limits = c(-20,20),aesthetics = "fill") +
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_size_continuous(limits = c(0,1),range = c(0, 6))
print(g)
dev.off()
