cancerTypes <- c("aero_digestive_tract_head and neck",
                 "aero_digestive_tract_oesophagus",
                 "blood_B_cell_lymphoma",
                 "blood_acute_myeloid_leukaemia",
                 "blood_B_cell_leukemia",
                 "bone_ewings_sarcoma",
                 "bone_osteosarcoma",
                 "breast_breast",
                 "digestive_system_large_intestine",
                 "digestive_system_stomach",
                 "digestive_system_liver",
                 "kidney_kidney",
                 "lung_lung_NSCLC_adenocarcinoma",
                 "lung_lung_small_cell_carcinoma",
                 "lung_mesothelioma",
                 "lung_lung_NSCLC_squamous_cell_carcinoma",
                 "nervous_system_glioma",
                 "nervous_system_neuroblastoma",
                 "pancreas_pancreas",
                 "skin_melanoma",
                 "thyroid_thyroid",
                 "urogenital_system_ovary",
                 "urogenital_system_bladder",
                 "urogenital_system_cervix",
                 "urogenital_system_endometrium",
                 "urogenital_system_prostate")

cor_hyper <- NULL
cor_perpval <- NULL
sample_num <- NULL
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
  y <- rbind(data1[1:2,samples],data2[2,samples],data3[2,samples],data4[2,samples])
  tmp_val <- NULL
  tmp_p <- NULL
  for (idx in 2:5){
    test <- cor.test(as.numeric(y[1,]),as.numeric(y[idx,]))
    tmp_val <- c(tmp_val,test$estimate)
    tmp_p <- c(tmp_p,test$p.value)
    
  }

  sample_num <- c(sample_num,ncol(y))
  cor_hyper <- rbind(cor_hyper,tmp_val)
  cor_perpval <- rbind(cor_perpval,tmp_p)
}

cor_hyper <- as.data.frame(cor_hyper)
rownames(cor_hyper) <- cancerTypes
colnames(cor_hyper) <- c("Cor.CGIhypo","Cor.L1","Cor.WCGW","Cor.Opensea")

cor_perpval <- as.data.frame(cor_perpval)
rownames(cor_perpval) <- cancerTypes
colnames(cor_perpval) <- c("pval.CGIhypo","pval.L1","pval.WCGW","pval.Opensea")

show <- cbind(cor_hyper,cor_perpval)
show$samplenum <- sample_num


show2 <- show[show$samplenum > 19,]
