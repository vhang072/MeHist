celltypes <- c("aero_digestive_tract_head and neck",
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

#pdf(file = "Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/_PDF_pic/normal_mean_hist_cellline.pdf",height = 5,width = 7)

path <- c("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_CCLE/_DiffCGIprobe_",
          "Y:/4.basic_data/TCGA_PancanAtlas/methylation_L1/Diff_CCLE/_DiffL1probe_",
          "Y:/4.basic_data/TCGA_PancanAtlas/methylation_soloWCGW/Diff_CCLE/_DiffsoloWCGWprobe_",
          "Y:/4.basic_data/TCGA_PancanAtlas/methylation_opensea/Diff_CCLE/_DiffOpenSeaprobe_")



diff_probes <- NULL
for (x in 1:length(celltypes)){
  tmp <- NULL
  for (y in 1:length(path)){
    data <- read.table(paste0(path[y],celltypes[x],".txt"),sep = "\t",
                     header = T)
    data$label <- "Control"
    data$label[data$Mean_T > data$Mean_N & data$AberrantFrac_T_0.15 > 0.5] <- "Global_Hyper"
    data$label[data$Mean_T < data$Mean_N & data$AberrantFrac_T_0.15 > 0.5] <- "Global_Hypo"
    
    tmp <- c(tmp,sum(data$label == "Global_Hyper"),
                                       sum(data$label == "Global_Hypo"),
                                       nrow(data))
  }
  diff_probes <- rbind(diff_probes,tmp)
}

diff_probes <- cbind(celltypes,diff_probes)
diff_probes <- as.data.frame(diff_probes)
colnames(diff_probes) <- c("cancerTypes","hyper_cgi","hypo_cgi","totalprobe_cgi",
                           "hyper_L1","hypo_L1","totalprobe_L1",
                           "hyper_soloWCGW","hypo_soloWCGW","totalprobe_soloWCGW",
                           "hyper_opensea","hypo_opensea","totalprobe_opensea")


write.table(diff_probes,"Y:/4.basic_data/TCGA_PancanAtlas/methylation_integrated/CCLE/_DiffProbeNumSummary.txt",
            sep = "\t",row.names = F,col.names = T, quote = F)
