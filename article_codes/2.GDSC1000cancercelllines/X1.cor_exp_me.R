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
  
#hypoxia signature
genelist1 <- c("ALDOA","MIF","TUBB6","P4HA1","SLC2A1","PGAM1","ENO1","LDHA","CDKN3","TPI1","NDRG1","VEGFA","ACOT7","ADM")
#tumor proliferation signature
genelist2 <- c("MKI67","NDC80","NUF2","PTTG1","RRM2","BIRC5","CCNB1","CEP55","UBE2C","CDC20","TYMS")
#cellcycle
library(stringr)
genelist2a <- read.table("I:/Projects/P7.Methy-PanCancer/genome_reference/CellCycle_reactome_R-HSA-164170/Cell_cycle_proteins_[R-HSA-1640170].tsv",
                         sep = "\t",header = T)
cycle <- matrix(unlist(str_split(genelist2a$MoleculeName," ")),ncol = 2,byrow = T)
genelist2b <- cycle[,2]

#methylation enzymes
genelist3 <- c("DNMT1","UHRF2","DNMT3A","DNMT3B","DNMT3L","TET1","TET2","TET3")


cors <- NULL
pvals <- NULL

cors2 <- NULL #num
pvals2 <- NULL #num
for (x in 1:length(cancerTypes)) {
  data1 <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_CCLE/_DiffCGIprobeMean_",cancerTypes[x],".txt"),
                      sep = "\t",header = T)
  data2 <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/RNA_CancerCellLines_rmGeneTitles_anotheranno/Exp_",cancerTypes[x],".txt"),
                      sep = "\t",header = T)
  colna <- NULL
  for (y in 2:ncol(data1)){
    colna <- c(colna,str_split(colnames(data1)[y],"_AVG.")[[1]][1])
  }
  colnames(data1)[2:ncol(data1)] <- colna
  samples <- intersect(colnames(data1)[2:ncol(data1)],colnames(data2)[2:ncol(data2)])
  
  
  d1 <- data1[,samples]
  d2 <- data2[is.na(data2$GENE_SYMBOLS)== F,samples]
  rownames(d2) <- data2$GENE_SYMBOLS[is.na(data2$GENE_SYMBOLS)==F]
  
  glist1 <- intersect(genelist1,rownames(d2))
  glist1_mat <- d2[glist1,]
  val1_mat <-  (glist1_mat - matrix(rep(apply(glist1_mat,1,mean),ncol(glist1_mat)),nrow(glist1_mat),byrow = F))/
    sqrt(matrix(rep(apply(glist1_mat,1,var),ncol(glist1_mat)),nrow(glist1_mat),byrow = F))
  val1 <- apply(val1_mat,2,mean)
  
  glist2 <- intersect(genelist2,rownames(d2))
  glist2_mat <- d2[glist2,]
  val2_mat <-  (glist2_mat - matrix(rep(apply(glist2_mat,1,mean),ncol(glist2_mat)),nrow(glist2_mat),byrow = F))/
    sqrt(matrix(rep(apply(glist2_mat,1,var),ncol(glist2_mat)),nrow(glist2_mat),byrow = F))
  val2 <- apply(val2_mat,2,mean)
 
  glist2a <- intersect(genelist2b,rownames(d2))
  glist2a_mat <- d2[glist2a,]
  val2a_mat <-  (glist2a_mat - matrix(rep(apply(glist2a_mat,1,mean),ncol(glist2a_mat)),nrow(glist2a_mat),byrow = F))/
    sqrt(matrix(rep(apply(glist2a_mat,1,var),ncol(glist2a_mat)),nrow(glist2a_mat),byrow = F))
  val2a <- apply(val2a_mat,2,mean) 
  
  glist3 <- intersect(genelist3,rownames(d2))
  glist3_mat <- d2[glist3,]
  val3_mat <-  glist3_mat - matrix(rep(apply(glist3_mat,1,mean),ncol(glist3_mat)),nrow(glist3_mat),byrow = F)
  
  
  causes_val <- rbind(d1,val1,val2,val2a,val3_mat)
  varnames <- c("Hyper","Hypo","Hyper.num","Hypo.num","Hypoxia","Proliferation","CellCycle",glist3)
  rownames(causes_val) <- varnames
  
  #mean value
  data <- t(causes_val)
  t.cor <- NULL
  t.pval <- NULL
  for (idx in c(2,5:ncol(data))){
    corr1 <- cor.test(data[,1],data[,idx])
    t.cor <- c(t.cor,corr1$estimate)
    t.pval <- c(t.pval,corr1$p.value)
  }
  
  cors <- rbind(cors,t.cor)
  pvals <- rbind(pvals,t.pval)
  
  #num
  t.cor <- NULL
  t.pval <- NULL
  for (idx in c(4,5:ncol(data))){
    corr1 <- cor.test(data[,3],data[,idx])
    t.cor <- c(t.cor,corr1$estimate)
    t.pval <- c(t.pval,corr1$p.value)
  }
  
  cors2 <- rbind(cors2,t.cor)
  pvals2 <- rbind(pvals2,t.pval)
}


#mean value
rownames(pvals) <- cancerTypes
colnames(pvals) <- paste0("P.",varnames[c(2,5:length(varnames))])

rownames(cors) <- cancerTypes
colnames(cors) <- paste0("Cor.",varnames[c(2,5:length(varnames))])

tmp.order <- order(as.numeric(cors[,1]),decreasing = T)
order.cancer <- rownames(cors)[tmp.order]

show <- NULL
for (x in 1:nrow(cors)){
  for (y in 1:ncol(cors)){
    show <- rbind(show,c(rownames(cors)[x],colnames(cors)[y],cors[x,y],pvals[x,y]))
  }
}
show <- as.data.frame(show)
colnames(show) <- c("Cancer","Cor.Item","Cor.Pearson","Pvalue")
show$Cor.Pearson <- as.numeric(show$Cor.Pearson)
show$Pvalue <- as.numeric(show$Pvalue)
show$Cancer <- factor(show$Cancer,levels = order.cancer)
show$Cor.Item <- factor(show$Cor.Item,levels = paste0("Cor.",varnames[c(2,5:length(varnames))]))

library(ggplot2)
library(scales)
show$sig <- sign(show$Cor.Pearson)*(-log10(show$Pvalue))
show$sig[show$sig>6] = 6
show$sig[show$sig< -6] = -6


library(paletteer)
pdf("I:/Projects/P7.Methy-PanCancer/pics/ccle_persample/ccle_coritems_mean_addUHRF1.pdf",height = 5,width = 6.5)
ggplot(show,aes(y=Cancer,x=Cor.Item,fill = sig,size = abs(Cor.Pearson))) + geom_point(shape = 21)+
  scale_color_gradientn(colors = paletteer_d("dichromat::DarkRedtoBlue_12")  ,limits = c(-6,6),aesthetics = "fill") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

  
#num
rownames(pvals2) <- cancerTypes
colnames(pvals2) <- paste0("P.",varnames[4:length(varnames)])

rownames(cors2) <- cancerTypes
colnames(cors2) <- paste0("Cor.",varnames[4:length(varnames)])

pvals2 <- as.data.frame(pvals2)
cors2 <- as.data.frame(cors2)
write.table(cors2, paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_CCLE/","___SigPur.",
                             "corrs",".CGIMean_",".txt"),sep = "\t",row.names = T, col.names = T,quote = F)
write.table(pvals2, paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_CCLE/","___SigPur.",
                              "pvals",".CGIMean_",".txt"),sep = "\t",row.names = T, col.names = T,quote = F)

tmp.order <- order(as.numeric(cors2[,1]),decreasing = T)
order.cancer <- rownames(cors2)[tmp.order]


show2 <- NULL
for (x in 1:nrow(cors2)){
  for (y in 1:ncol(cors2)){
    show2 <- rbind(show2,c(rownames(cors2)[x],colnames(cors2)[y],cors2[x,y],pvals2[x,y]))
  }
}
show2 <- as.data.frame(show2)
colnames(show2) <- c("Cancer","Cor.Item","Cor.Pearson","Pvalue")
show2$Cor.Pearson <- as.numeric(show2$Cor.Pearson)
show2$Pvalue <- as.numeric(show2$Pvalue)
show2$Cancer <- factor(show2$Cancer,levels = order.cancer)
show2$Cor.Item <- factor(show2$Cor.Item,levels = paste0("Cor.",varnames[4:length(varnames)]))
write.table(show2,"Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_CCLE/ccle_hyper_corritemss.txt",
            sep = '\t',quote = F,row.names = F)
library(ggplot2)
library(scales)
show2$sig <- sign(show2$Cor.Pearson)*(-log10(show2$Pvalue))
show2$sig[show2$sig>6] = 6
show2$sig[show2$sig< -6] = -6


library(paletteer)
pdf("I:/Projects/P7.Methy-PanCancer/pics/ccle_persample/ccle_coritems_num_addUHRF1.pdf",height = 5,width = 6.5)
ggplot(show2,aes(y=Cancer,x=Cor.Item,fill = sig,size = abs(Cor.Pearson))) + geom_point(shape = 21)+
  scale_color_gradientn(colors = paletteer_d("dichromat::DarkRedtoBlue_12")  ,limits = c(-6,6),aesthetics = "fill") +
  theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
