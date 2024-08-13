cancerTypes <- c(
  "LAML",
  "ACC",
  "BLCA",
  "BRCA",
  "CESC",
  "CHOL",
  "COAD",
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
fold <- "ABSOLUTE"
#hypoxia signature
genelist1 <- c("ENSG00000149925|ALDOA","ENSG00000240972|MIF","ENSG00000176014|TUBB6","ENSG00000122884|P4HA1",
               "ENSG00000117394|SLC2A1","ENSG00000171314|PGAM1","ENSG00000074800|ENO1","ENSG00000134333|LDHA",
               "ENSG00000100526|CDKN3","ENSG00000111669|TPI1","ENSG00000104419|NDRG1","ENSG00000112715|VEGFA",
               "ENSG00000097021|ACOT7","ENSG00000148926|ADM")
#tumor proliferation signature
genelist2 <- c("ENSG00000148773|MKI67","ENSG00000080986|NDC80","ENSG00000143228|NUF2","ENSG00000164611|PTTG1",
               "ENSG00000171848|RRM2","ENSG00000089685|BIRC5","ENSG00000134057|CCNB1","ENSG00000138180|CEP55",
               "ENSG00000175063|UBE2C","ENSG00000117399|CDC20","ENSG00000176890|TYMS")

#cell cycle signature
library(stringr)
genelist2a <- read.table("I:/Projects/P7.Methy-PanCancer/genome_reference/CellCycle_reactome_R-HSA-164170/Cell_cycle_proteins_[R-HSA-1640170].tsv",
                         sep = "\t",header = T)
cycle <- matrix(unlist(str_split(genelist2a$MoleculeName," ")),ncol = 2,byrow = T)


#methylation enzymes
genelist3 <- c("ENSG00000130816|DNMT1","ENSG00000276043|UHRF1","ENSG00000147854|UHRF2",
               "ENSG00000119772|DNMT3A","ENSG00000088305|DNMT3B",
               "ENSG00000142182|DNMT3L",
               "ENSG00000138336|TET1","ENSG00000168769|TET2","ENSG00000187605|TET3"
               )
#"ENSG00000065060|UHRF1BP1","ENSG00000111647|UHRF1BP1L"
mean_rmna <- function(x){
  return(mean(x,na.rm = T)) 
}
var_rmna <- function(x){
  return(var(x,na.rm = T)) 
}

chose_most_enriched_region <- function(x1){
  x1[is.na(x1)] <- 0
  range_store_x1 <- NULL
  for (x in seq(0,0.99,0.01)){
    for (y in seq(x+0.01,1,0.01)){
      if (sum((x1>=x) & (x1 <= y))>=0.5*length(x1)){
        range_store_x1 <- rbind(range_store_x1,c(x,y))
        break
      }
    }
  }
  range <- range_store_x1[,2] - range_store_x1[,1]
  ind <- order(range,decreasing = T)
  s <- range_store_x1[ind[length(ind)],1]
  e <- range_store_x1[ind[length(ind)],2]
  return(list(s,e))
}

cor_mat <- NULL # probe mean val
pval_mat <- NULL # probe mean val
cor_mat2 <- NULL #probe num
pval_mat2 <- NULL #probe num

used_cancers <- NULL
for (x in 1:length(cancerTypes)){
  data <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_TCGA/",fold,"_MeanCGIprobe/",
                            cancerTypes[x],".CGIMean_",fold,".txt"),sep = "\t",header = T)
  exp <- read.table(paste0("Y:/4.basic_data/TCGA/RNAseq/FPKM/Primary Tumor/TCGA-",cancerTypes[x],".txt"),
                    sep = "\t",header = T)
  rownames(exp) <- exp$gene_id
  exp <- exp[,2:ncol(exp)]
  
  exp_colnames <- NULL
  used_cols <- NULL
  library(stringr)
  for (idx in 1:ncol(exp)){
    sp_name <- substr(colnames(exp)[idx],1,15)
    if (sp_name %in% exp_colnames == F){
      exp_colnames <- c(exp_colnames, sp_name)
      used_cols <- c(used_cols,idx)
    }
  }
  exp_used <- exp[,used_cols]
  colnames(exp_used) <- exp_colnames
  overlap_sp <- intersect(colnames(data),colnames(exp_used))
  
  if (length(overlap_sp)>=20){
    used_cancers <- c(used_cancers,cancerTypes[x] )
    data_now <- data[,overlap_sp]
    exp_now <- exp_used[,overlap_sp]
    
    
    #exp of singatures,hypoxia
    g1 <- exp_now[intersect(genelist1,rownames(exp_now)),]
    mean_mat <- matrix(rep(apply(g1,1,mean_rmna),ncol(g1)),nrow = nrow(g1),byrow = F)
    var_mat <- matrix(rep(apply(g1,1,var_rmna),ncol(g1)),nrow = nrow(g1),byrow = F)
    g1_zscore <- (g1 - mean_mat)/sqrt(var_mat)
    g1_mean <- apply(g1_zscore,2,mean_rmna)
    
    #exp of singatures,proliferation
    g2 <- exp_now[intersect(genelist2,rownames(exp_now)),]
    mean_mat <- matrix(rep(apply(g2,1,mean_rmna),ncol(g2)),nrow = nrow(g2),byrow = F)
    var_mat <- matrix(rep(apply(g2,1,var_rmna),ncol(g2)),nrow = nrow(g2),byrow = F)
    g2_zscore <- (g2 - mean_mat)/sqrt(var_mat)
    g2_mean <- apply(g2_zscore,2,mean_rmna)
    
    #exp of cellcycle
    exp_split <- matrix(unlist(str_split(rownames(exp_now),"\\|")),ncol = 2,byrow = T)
    g2a <- exp_now[exp_split[,2] %in% cycle[,2],]
    mean_mat <- matrix(rep(apply(g2a,1,mean_rmna),ncol(g2a)),nrow = nrow(g2a),byrow = F)
    var_mat <- matrix(rep(apply(g2a,1,var_rmna),ncol(g2a)),nrow = nrow(g2a),byrow = F)
    g2a_zscore <- (g2a - mean_mat)/sqrt(var_mat)
    g2a_mean <- apply(g2a_zscore,2,mean_rmna)    
    
    #exp of singatures
    g3 <- exp_now[genelist3,]
    merge <- rbind(data_now,g1_mean,g2_mean,g2a_mean,g3)
    
    #deal with tumor purity
    interval <- chose_most_enriched_region(merge[5,])
    label <- rep("Excluded",ncol(merge))
    label[(merge[5,]>= interval[[1]][1]) & (merge[5,]<= interval[[2]][1])] <- "Included"
    
    out <- t(merge)
    out <- cbind(out,label)
    
    colnames(out)[1:8] <- c("Hyper","Hypo","Hyper.Num","Hypo.Num","Purity","Hypoxia","Prolif.","CellCycle")
    
    write.table(out, paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_TCGA/",fold,"_MeanCGIprobe/_SigPur.",
                            cancerTypes[x],".CGIMean_",fold,".txt"),sep = "\t",row.names = T, col.names = T,quote = F)
    
    
    #cor test value
    cor_line <- NULL
    pval_line <- NULL
    for (idx in c(2,6:17)){
      corr1 <- cor.test(as.numeric(out[,1]),as.numeric(out[,idx]))
      cor_line <- c(cor_line,corr1$estimate)
      pval_line <- c(pval_line,corr1$p.value)
    }
    out_f <- out[out[,18] == "Included",]
    for (idx in c(2,6:17)){
      corr1 <- cor.test(as.numeric(out_f[,1]),as.numeric(out_f[,idx]))
      cor_line <- c(cor_line,corr1$estimate)
      pval_line <- c(pval_line,corr1$p.value)
    }
    cor_mat <- rbind(cor_mat,cor_line)
    pval_mat <- rbind(pval_mat,pval_line)
    
    #cor test num
    cor_line2 <- NULL
    pval_line2 <- NULL
    for (idx in c(4,6:17)){
      corr1 <- cor.test(as.numeric(out[,3]),as.numeric(out[,idx]))
      cor_line2 <- c(cor_line2,corr1$estimate)
      pval_line2 <- c(pval_line2,corr1$p.value)
    }
    out_f <- out[out[,18] == "Included",]
    for (idx in c(4,6:17)){
      corr1 <- cor.test(as.numeric(out_f[,3]),as.numeric(out_f[,idx]))
      cor_line2 <- c(cor_line2,corr1$estimate)
      pval_line2 <- c(pval_line2,corr1$p.value)
    }
    cor_mat2 <- rbind(cor_mat2,cor_line2)
    pval_mat2 <- rbind(pval_mat2,pval_line2)
    
  }
  
}


cols.t <- c("Hypo","Hypoxia","Proliferation","CellCycle","DNMT1","UHRF1","UHRF2","DNMT3A","DNMT3B","DNMT3L","TET1","TET2","TET3")
cols <- c(cols.t,paste0(cols.t,".f"))          
#mean value
colnames(cor_mat) <- cols
rownames(cor_mat) <- used_cancers
cor_mat <- as.data.frame(cor_mat)

print(t.test(cor_mat$Hypo,cor_mat$Hypo.f,paired = T))

colnames(pval_mat) <- cols
rownames(pval_mat) <- used_cancers
pval_mat <- as.data.frame(pval_mat)


tmp.order <- order(as.numeric(cor_mat[,13]),decreasing = T)
order.cancer <- rownames(cor_mat)[tmp.order]

show <- NULL
for (x in 1:nrow(cor_mat)){
  for (y in 1:ncol(cor_mat)){
    show <- rbind(show,c(rownames(cor_mat)[x],colnames(cor_mat)[y],cor_mat[x,y],pval_mat[x,y]))
  }
}
show <- as.data.frame(show)
colnames(show) <- c("Cancer","Cor.Item","Cor.Pearson","Pvalue")
show$Cor.Pearson <- as.numeric(show$Cor.Pearson)
show$Pvalue <- as.numeric(show$Pvalue)
show$Cancer <- factor(show$Cancer,levels = order.cancer)
show$Cor.Item <- factor(show$Cor.Item,levels = cols)

library(ggplot2)
library(scales)
show$sig <- sign(show$Cor.Pearson)*(-log10(show$Pvalue))
show$sig[show$sig>20] = 20
show$sig[show$sig< -20] = -20

library(paletteer)
pdf("I:/Projects/P7.Methy-PanCancer/pics/tcga_purity_pertumor/meanval_corr_items_absolute_addonlyUHRF1_2.pdf",height = 7,width = 7)
ggplot(show,aes(y=Cancer,x=Cor.Item,fill = sig,size = abs(Cor.Pearson))) + geom_point(shape = 21)+
  scale_color_gradientn(colors = paletteer_d("dichromat::DarkRedtoBlue_12")  ,limits = c(-20,20),aesthetics = "fill") +
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


#num
colnames(cor_mat2) <- cols
rownames(cor_mat2) <- used_cancers
cor_mat2 <- as.data.frame(cor_mat2)
print(t.test(cor_mat2$Hypo,cor_mat2$Hypo.f,paired = T))

colnames(pval_mat2) <- cols
rownames(pval_mat2) <- used_cancers
pval_mat2 <- as.data.frame(pval_mat2)

write.table(cor_mat2, paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_TCGA/",fold,"_MeanCGIprobe/___SigPur.",
                        "corrs",".CGIMean_",fold,".txt"),sep = "\t",row.names = T, col.names = T,quote = F)
write.table(pval_mat2, paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_TCGA/",fold,"_MeanCGIprobe/___SigPur.",
                        "pvals",".CGIMean_",fold,".txt"),sep = "\t",row.names = T, col.names = T,quote = F)

tmp.order <- order(as.numeric(cor_mat2[,13]),decreasing = T)
order.cancer <- rownames(cor_mat2)[tmp.order]

show2 <- NULL
for (x in 1:nrow(cor_mat2)){
  for (y in 1:ncol(cor_mat2)){
    show2 <- rbind(show2,c(rownames(cor_mat2)[x],colnames(cor_mat2)[y],cor_mat2[x,y],pval_mat2[x,y]))
  }
}
show2 <- as.data.frame(show2)
colnames(show2) <- c("Cancer","Cor.Item","Cor.Pearson","Pvalue")
show2$Cor.Pearson <- as.numeric(show2$Cor.Pearson)
show2$Pvalue <- as.numeric(show2$Pvalue)
show2$Cancer <- factor(show2$Cancer,levels = order.cancer)
show2$Cor.Item <- factor(show2$Cor.Item,levels = cols)

library(ggplot2)
library(scales)
show2$sig <- sign(show2$Cor.Pearson)*(-log10(show2$Pvalue))
show2$sig[show2$sig>20] = 20
show2$sig[show2$sig< -20] = -20

library(paletteer)
pdf("I:/Projects/P7.Methy-PanCancer/pics/tcga_purity_pertumor/num_corr_items_absolute_addonlyUHRF1_2.pdf",height = 7,width = 7)
ggplot(show2,aes(y=Cancer,x=Cor.Item,fill = sig,size = abs(Cor.Pearson))) + geom_point(shape = 21)+
  scale_color_gradientn(colors = paletteer_d("dichromat::DarkRedtoBlue_12")  ,limits = c(-20,20),aesthetics = "fill") +
  theme_bw()+theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# ggplot(show,aes(y=Cancer,x=Cor.Item,fill = sig,size = abs(Cor.Pearson))) + geom_point(shape = 21)+
#   scale_color_gradient2(low = "blue", mid = "white", high = "red",midpoint = 0,limits = c(-20,20),aesthetics = "fill") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

