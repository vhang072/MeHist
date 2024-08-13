cancerTypes <- c( "ACC",
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

median_rm <- function(x){
  y <- median(x,na.rm = T)
  return(y)
}

exp_median <- NULL
for (idx in 1:length(cancerTypes)){
  exp <- read.table(paste0("Y:/4.basic_data/TCGA/RNAseq/FPKM/Primary Tumor/TCGA-",cancerTypes[idx],".txt"),
                    sep = "\t",header = T)
  
  exp_median <- cbind(exp_median,apply(exp[,2:ncol(exp)],1,median_rm))
}


me_median <- NULL
for (idx in 1:length(cancerTypes)){
  data0 <- read.table(paste0("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_TCGA/_DiffCGIprobeMean_",cancerTypes[idx],".txt"),
                      sep = "\t",header = T)
  data <- data0[,2:ncol(data0)]
  median_num_hyper <- median_rm(as.numeric(unname(data[3,])))
  median_num_hypo <- median_rm(as.numeric(unname(data[4,])))
  
  me_median <- rbind(me_median,c(median_num_hyper,median_num_hypo))
}


diff_res <- NULL
for (idx in 1:nrow(exp_median)){
  test1 <- cor.test(exp_median[idx,],me_median[,1],method = "spearman")
  test2 <- cor.test(exp_median[idx,],me_median[,2],method = "spearman")
  diff_res <- rbind(diff_res, c(test1$estimate,test1$p.value,test2$estimate,test2$p.value))
}
library(stringr)
gene_mat <- matrix(unlist(str_split(exp[,1],"\\|")),ncol = 2,byrow = T)
out <- cbind(gene_mat,diff_res)

out <- as.data.frame(out)
colnames(out) <- c("geneID","symbol","HyperCor.spearman","HyperCor.P","HypoCor.spearman","HypoCor.P")
for (idx in 3:6){
  out[,idx] <- as.numeric(out[,idx])
}
out <- na.omit(out)
write.table(out,"Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Pancancer_ExpMe_TCGA/pancancer_hyperhypoNumAsso.spearman.txt",
            sep = "\t",row.names = F, col.names = T,quote = T)

cor.test(me_median[,1],me_median[,2])
cor.test(me_median[c(1:12,14:16,18:25),1],me_median[c(1:12,14:16,18:25),2])

genename <- "ENSG00000182481|KPNA2"
genename <- "ENSG00000087586|AURKA"
show <- data.frame(gene = exp_median[exp[,1] == genename,],hyper = me_median[,1],hypo = me_median[,2],cancer = cancerTypes)
library(ggplot2)
pdf("I:/Projects/P7.Methy-PanCancer/pics/Metascape/AURKA.pdf",height = 4,width = 4.5)
g <- ggplot(show, aes(x = gene,y = hyper)) + geom_point(shape = 21, fill = "#6a5acd",size = 2) + 
  geom_text(aes(x = gene,y = hyper + 500,label =cancerTypes),label.size = NA) +
  geom_smooth(method='lm') + xlim(c(0,22)) + 
  theme_bw()

print(g)

g <- ggplot(show, aes(x = gene,y = hypo)) + geom_point(shape = 21, fill = "#6a5acd",size = 2) + 
  geom_text(aes(x = gene,y = hypo + 100,label =cancerTypes),label.size = NA) +
  geom_smooth(method='lm') + xlim(c(0,22)) +
theme_bw()

print(g)
dev.off()
