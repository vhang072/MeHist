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

cor_hypo <- NULL
pval_hypo <- NULL
purity <- read.table("Y:/4.basic_data/TCGA_PancanAtlas/tumor_purity/MethylResolver/methylresolver.txt",
                     sep = "\t",header = T)
purity_rmna <- purity[is.na(purity$TumorPurity) == F,]
purity_rownames <- NULL
idx_used <- NULL
library(stringr)
for (idx in 1:nrow(purity_rmna)){
  tmp.name <- substr(purity_rmna$Sample[idx],1,16)
  if (tmp.name %in% purity_rownames == F){
  purity_rownames <- c(purity_rownames,tmp.name)
  idx_used <- c(idx_used,idx)
  }
}
purity_rmna <- purity_rmna[idx_used,]
rownames(purity_rmna) <- purity_rownames



for (x in 1:length(cancerTypes)) {
  data1 <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_TCGA/_DiffCGIprobeMean_",cancerTypes[x],".txt"),
                      sep = "\t",header = T)
  sp_ids <- NULL
  sp_used <- NULL
  for (idx in 2:ncol(data1)){
    sp_now <- substr(colnames(data1)[idx],1,16)
    if (sp_now %in% sp_ids == F){
      sp_ids <- c(sp_ids,sp_now)
      sp_used <- c(sp_used,idx)
    }
  }
  data1_curated <- data1[,sp_used]
  colnames(data1_curated) <- sp_ids
  overlap_sps <- intersect(rownames(purity_rmna),colnames(data1_curated))
  if (length(overlap_sps)>=10){
    purity_val <- purity_rmna[overlap_sps,"TumorPurity"]
    data_val <- data1_curated[,overlap_sps]
    out <- rbind(data_val,purity_val)
    out <- as.data.frame(out)
    write.table(out,paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_TCGA/Methylresolver_MeanCGIprobe/",
                           cancerTypes[x],".CGIMean_Methylresolver.txt"),sep = "\t",row.names = F, col.names = T)
  }
  
}
