cancerTypes <- c("ESCA","GBM","LAML","LGG","LUSC","STAD","UCEC","PCPG")
library(stringr)

out_muts <- NULL
for (idx in 1:length(cancerTypes)){
  subclu <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Survival_hyper/HyperDichotomy_",cancerTypes[idx],".txt"),
                       sep = "\t",header = T)
  muta0 <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/mutations/",cancerTypes[idx],"_mut_matrix.txt"),sep = "\t",header = T)
  rownames(muta0) <- muta0$genename
  muta <- muta0[,2:ncol(muta0)]
  spname <- substr(colnames(muta),1,12)
  colnames(muta) <- spname
  
  overlaps <- intersect(rownames(subclu),colnames(muta))
  cluf <- subclu[overlaps,]
  mutf <- muta[,overlaps]
  mutf[mutf>1] <- 1
  
  for (x in 1:nrow(mutf)){
    if (sum(mutf[x,])>=3){
      t <- table(factor(as.numeric(mutf[x,])),cluf[,"hyper_label"])
      test1 <- fisher.test(t)
      out_muts <- rbind(out_muts,c(cancerTypes[idx],rownames(mutf)[x],test1$p.value,t[1,1],t[1,2],t[2,1],t[2,2]))
    }
  }
  
}
out_muts <- as.data.frame(out_muts)
colnames(out_muts) <- c("Cancer","Gene","Pvalue","T1.1","T1.2","T2.1","T2.2")

write.table(out_muts,"Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/mutation_hyper/mutation_hyper.txt",sep = "\t",
            row.names = F,quote = F)


##draw
ct<- "STAD"
gz <- "PIK3CA"
subclu <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Survival_hyper/HyperDichotomy_",ct,".txt"),
                     sep = "\t",header = T)
muta0 <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/mutations/",ct,"_mut_matrix.txt"),sep = "\t",header = T)
rownames(muta0) <- muta0$genename
muta <- muta0[,2:ncol(muta0)]
spname <- substr(colnames(muta),1,12)
colnames(muta) <- spname

overlaps <- intersect(rownames(subclu),colnames(muta))
cluf <- subclu[overlaps,]
mutf <- muta[gz,overlaps]
mutf[mutf>1] <- 1
cluf$mut <- factor(as.numeric(mutf))

library(ggplot2)
pdf(paste0("I:/Projects/P7.Methy-PanCancer/pics/mutation_virus/",ct,"-",gz,".pdf"),height = 4,width = 5)
g <- ggplot(cluf,aes(x = hypo_num,y = hyper_num,fill = mut)) + geom_point(shape = 21,color = "black",size = 2) +
  theme_bw() + scale_fill_manual(values=c("#999999", "#800000")) + ggtitle(paste0(ct,"-",gz))
print(g)
dev.off()
