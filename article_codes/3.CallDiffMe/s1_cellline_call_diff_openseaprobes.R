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
prefix_file <- "_OpenseaProbe_"
normaltypes <- c(paste0("SolidNormal/",prefix_file,"TCGA-HNSC"),
                 paste0("SolidNormal/",prefix_file,"TCGA-ESCA"),
                 paste0("GEONormal/",prefix_file,"Blood_GSE68456_Bcell"),
                 paste0("GEONormal/",prefix_file,"LAML_GSE58477_bonemarrow10"),
                 paste0("GEONormal/",prefix_file,"Blood_GSE68456_Bcell"),
                 paste0("SolidNormal/",prefix_file,"TCGA-SARC"),
                 paste0("SolidNormal/",prefix_file,"TCGA-SARC"),
                 paste0("SolidNormal/",prefix_file,"TCGA-BRCA"),
                 paste0("SolidNormal/",prefix_file,"TCGA-COAD"),
                 paste0("SolidNormal/",prefix_file,"TCGA-STAD"),#STAD_GSE211580_mucosae9
                 paste0("SolidNormal/",prefix_file,"TCGA-LIHC"),
                 paste0("SolidNormal/",prefix_file,"TCGA-KIRP"),
                 paste0("SolidNormal/",prefix_file,"TCGA-LUAD"),
                 paste0("SolidNormal/",prefix_file,"TCGA-LUAD"),
                 paste0("SolidNormal/",prefix_file,"TCGA-LUAD"),
                 paste0("SolidNormal/",prefix_file,"TCGA-LUAD"),
                 paste0("SolidNormal/",prefix_file,"TCGA-GBM"),#GBM_GSE79122_controlbrain10
                 paste0("SolidNormal/",prefix_file,"TCGA-GBM"),#GBM_GSE79122_controlbrain10
                 paste0("SolidNormal/",prefix_file,"TCGA-PAAD"),
                 paste0("SolidNormal/",prefix_file,"TCGA-SKCM"),#SKCM_GSE120878_nevus73
                 paste0("SolidNormal/",prefix_file,"TCGA-THCA"),
                 paste0("GEONormal/",prefix_file,"OV_GSE146552_ovariesepithelial8"),
                 paste0("SolidNormal/",prefix_file,"TCGA-BLCA"),
                 paste0("SolidNormal/",prefix_file,"TCGA-CESC"),
                 paste0("SolidNormal/",prefix_file,"TCGA-UCEC"),
                 paste0("SolidNormal/",prefix_file,"TCGA-PRAD"))

path1 <- "Y:/4.basic_data/TCGA_PancanAtlas/methylation_opensea/CellLines/_OpenseaProbe_450_allprobe_split_"
path2 <- "Y:/4.basic_data/TCGA_PancanAtlas/methylation_opensea/"
# call differential methylation probes
#diff threshold 
thr_me <- c(0.1, 0.15, 0.2, 0.25, 0.3)#
used_lines <- 1:length(celltypes)#c(10,17,18,20)
for (x in used_lines){
  print(paste0("Processing File: ", celltypes[x]))
  tumor_part <- read.table(paste0(path1,celltypes[x],".txt"),sep = "\t",header = T)
  normal_part <- read.table(paste0(path2,normaltypes[x],".txt"),sep = "\t",header = T)
  rownames(tumor_part) <- tumor_part$probeid
  rownames(normal_part) <- normal_part$probeid
  overlapprobes <- intersect(rownames(tumor_part),rownames(normal_part))
  tumor_select <- tumor_part[overlapprobes,]
  normal_select <- normal_part[overlapprobes,]
  
  diffprobeinfo <- NULL
  for (idx in 1:nrow(tumor_select)){
    tumorme <- tumor_select[idx,2:ncol(tumor_select)]
    normalme <- normal_select[idx,2:ncol(normal_select)]
    tumorselect <- tumorme[is.na(tumorme)==F]
    normalselect <- normalme[is.na(normalme)==F]
    num_tumor <- length(tumorselect)
    num_normal <- length(normalselect)
    if (num_tumor >=2 & num_normal>=2){
      mean_tumor <- mean(tumorselect)
      mean_normal <- mean(normalselect)
      test1 <- wilcox.test(tumorselect, normalselect,exact = F)
      p1 <- test1$p.value
      #test2 <- ks.test(tumorselect, normalselect,alternative = "less",exact = F)
      #p2 <- test2$p.value
      #test3 <- ks.test(tumorselect, normalselect,alternative = "greater",exact = F)
      #p3 <- test3$p.value
      abe_frac_t <- NULL
      abe_frac_n <- NULL
      if (mean_tumor > mean_normal){
        for (thr in thr_me){
          abe_frac_t <- c(abe_frac_t, sum(tumorselect > (mean_normal+thr))/num_tumor)#
          abe_frac_n <- c(abe_frac_n, sum(normalselect > (mean_normal+thr))/num_normal)#
        }
      } else {
        for (thr in thr_me){
          abe_frac_t <- c(abe_frac_t, sum(tumorselect < (mean_normal-thr))/num_tumor)#
          abe_frac_n <- c(abe_frac_n, sum(normalselect < (mean_normal-thr))/num_normal)#
        }
      }
    } else {
      mean_tumor <- NA
      mean_normal <- NA
      p1 <- NA
      #p2 <- NA
      #p3 <- NA
      abe_frac_t <- rep(NA,length(thr_me))#
      abe_frac_n <- rep(NA,length(thr_me))#
    }
    diffprobeinfo <- rbind(diffprobeinfo, c(num_tumor, num_normal, mean_tumor, mean_normal,p1,abe_frac_t,abe_frac_n))#
  }
  diffprobeinfo <- as.data.frame(diffprobeinfo)
  thr_names <- c(paste0("AberrantFrac_T_",thr_me),paste0("AberrantFrac_N_",thr_me))#
  colnames(diffprobeinfo) <- c("Num_T","Num_N","Mean_T","Mean_N","P.wilcox",thr_names)#
  
  probeid <- tumor_select[,1]
  merged <- cbind(probeid,diffprobeinfo)
  write.table(merged,
              paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_opensea/Diff_CCLE/_DiffOpenSeaprobe_",celltypes[x],".txt"),
              sep = "\t", col.names = T, row.names = F, quote = F)
}


