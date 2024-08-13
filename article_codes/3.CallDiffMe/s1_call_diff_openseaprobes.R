path1 <- "Y:/4.basic_data/TCGA_PancanAtlas/methylation_opensea/LAML/"
path2 <- "Y:/4.basic_data/TCGA_PancanAtlas/methylation_opensea/Primary/"
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
tumorFiles <- c(paste0(path1,"_OpenseaProbe_TCGA-",cancerTypes[1],".txt"))
for (x in 2:length(cancerTypes)){
  tumorFiles <- c(tumorFiles,paste0(path2,"_OpenseaProbe_TCGA-",cancerTypes[x],".txt"))
}
lostTypes = c("MESO","TGCT","UCS","UVM")#without matched normal 
uncertainTypes = c("DLBC","GBM","KICH","SKCM","STAD") # used outside normal data
path3 = "Y:/4.basic_data/TCGA_PancanAtlas/methylation_opensea/SolidNormal/"
path4 = "Y:/4.basic_data/TCGA_PancanAtlas/methylation_opensea/GEONormal/"
solidTypes = c("LAML_GSE58477_bonemarrow10",
               "TCGA-PCPG",#ACC_GSE77871_normaladrenal6
               "TCGA-BLCA",
               "TCGA-BRCA",
               "TCGA-CESC",
               "TCGA-CHOL",
               "TCGA-COAD",
               "Blood_GSE68456_Bcell",
               "TCGA-ESCA",
               "TCGA-GBM",#GBM_GSE79122_controlbrain10
               "TCGA-HNSC",
               "TCGA-KIRP",
               "TCGA-KIRC",
               "TCGA-KIRP",
               "TCGA-GBM",#GBM_GSE79122_controlbrain10
               "TCGA-LIHC",
               "TCGA-LUAD",
               "TCGA-LUSC",
               "OV_GSE146552_ovariesepithelial8",
               "TCGA-PAAD",
               "TCGA-PCPG",
               "TCGA-PRAD",
               "TCGA-READ",
               "TCGA-SARC",
               "TCGA-SKCM",#SKCM_GSE120878_nevus73
               "TCGA-STAD",#STAD_GSE211580_mucosae9
               "TCGA-THCA",
               "TCGA-THYM",
               "TCGA-UCEC")
normalFiles <- NULL
for (x in 1:length(solidTypes)){
  if (nchar(solidTypes[x])>10){
    normalFiles <- c(normalFiles,paste0(path4,"_OpenseaProbe_",solidTypes[x],".txt"))
  } else{
    normalFiles <- c(normalFiles,paste0(path3,"_OpenseaProbe_",solidTypes[x],".txt"))
  }
}

# call differential methylation probes
#diff threshold 
thr_me <- c(0.1, 0.15, 0.2, 0.25, 0.3)#
used_lines <- 1:length(cancerTypes)#c(2,10,15,25,26)
for (x in used_lines){
  print(paste0("Processing File: ", tumorFiles[x]))
  tumor_part <- read.table(tumorFiles[x],sep = "\t",header = T)
  normal_part <- read.table(normalFiles[x],sep = "\t",header = T)
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
              paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_opensea/Diff_TCGA/_DiffOpenSeaprobe_",cancerTypes[x],".txt"),
              sep = "\t", col.names = T, row.names = F, quote = F
  )
}