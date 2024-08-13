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

used_lines <- 1:length(cancerTypes)

mean_rmna <- function(x){
  return(mean(x,na.rm = T)) 
}

cor_mean <- NULL

for (i in used_lines){
  info <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_CCLE/_DiffCGIprobe_",cancerTypes[i],".txt"),
                     sep = "\t",header = T)
  data <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/CellLines/_cgiProbe_450_allprobe_split_",cancerTypes[i],".txt"),
                     sep = "\t",header = T)
  data_colrange <- 5:ncol(data)
  hyper <- info[(info$Mean_T > info$Mean_N) & (info$AberrantFrac_T_0.2 > 0.5)  & (info$Mean_N < 0.3),]
  hypo <- info[(info$Mean_T < info$Mean_N) & (info$AberrantFrac_T_0.2 > 0.5) & (info$Mean_N > 0.7) ,]
  
  hyper_probe <- data[data$probeid %in% hyper$probeid,]
  hypo_probe <- data[data$probeid %in% hypo$probeid,]
  
  
  hyper_mean <- apply(hyper_probe[,data_colrange],2,mean_rmna)- mean(hyper$Mean_N,na.rm = T)
  hypo_mean <- mean(hypo$Mean_N,na.rm = T) - apply(hypo_probe[,data_colrange],2,mean_rmna)
  

  x1 <- cor((hyper_mean),t(hyper_probe[,data_colrange]))
  hist(x1,breaks = seq(-1,1,0.05))
  x2 <- cor((hypo_mean),t(hypo_probe[,data_colrange]))
  hist(x2,breaks = seq(-1,1,0.05))
  x1[is.na(x1)] <- 0
  x2[is.na(x2)] <- 0
  
  cor_mean <- rbind(cor_mean,c(mean(x1),mean(x2)))
  mean_vals <- rbind(hyper_mean,hypo_mean)
  mean_vals <- cbind(c("hyper_mean","hypo_mean"),mean_vals)
  mean_vals <- as.data.frame(mean_vals)
  colnames(mean_vals) <- c("mean_type",colnames(data)[data_colrange])
  write.table(mean_vals,paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_CCLE/_DiffCGIprobeMean_",cancerTypes[i],".txt"),
              sep = "\t",row.names = F,col.names = T, quote = F)
  }
cor_mean <- cbind(cancerTypes,cor_mean)
cor_mean <- as.data.frame(cor_mean)
colnames(cor_mean) <- c("cancer_types","hyper_cor_mean","hypo_cor_mean")

write.table(cor_mean,paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_CCLE/__DiffCGIprobeCorSummary.txt"),
            sep = "\t",row.names = F,col.names = T, quote = F)

