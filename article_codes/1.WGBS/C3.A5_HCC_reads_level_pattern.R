path <- "I:/Projects/P7.Methy-PanCancer/data_analysis/A5.HCC_WGBS/Mehist/"
sampleID1 <- read.table("I:/Projects/P7.Methy-PanCancer/data_analysis/A5.HCC_WGBS/A5.HCC_normal.txt",sep = "\t",header = F)
sampleID2 <- read.table("I:/Projects/P7.Methy-PanCancer/data_analysis/A5.HCC_WGBS/A5.HCC_tumor.txt",sep = "\t",header = F)


for (idx in 1:nrow(sampleID1)){
  data <- read.table(paste0(path,sampleID1[idx,1],"_diffcgis_Mehist.txt"),sep = "\t",header = T)
  if (idx == 1){
    normalhist <- data
  } else{
    if (sum(data$BedID == normalhist$BedID)==nrow(data)){
      normalhist[,2:11] <- normalhist[,2:11] + data[,2:11]
    }
  }
}

normaldensity <- normalhist
for (idx in 1:nrow(normaldensity)){
  normaldensity[idx,2:11] <- normaldensity[idx,2:11]/sum(normaldensity[idx,2:11])
}

hyper_rg <- 1:1127
hypo_rg <- 1128:2980

hyper.normal.hist <- apply(normaldensity[hyper_rg,2:11],2,mean)
hyper.normal.hist <- hyper.normal.hist/sum(hyper.normal.hist)
hypo.normal.hist <- apply(normaldensity[hypo_rg,2:11],2,mean)
hypo.normal.hist <- hypo.normal.hist/sum(hypo.normal.hist)
normal.hist <- data.frame(hyper = hyper.normal.hist,hypo = hypo.normal.hist,x = seq(1,10),row.names = 1:10)


mean_rmna <- function(x){
  y <- mean(x,na.rm = T)
  return(y)
}

#calculate hyper and hypo distributions in one tumor samples
for (idx in c("SRR15860722","SRR15860734","SRR15860750")){
  data <- read.table(paste0(path,idx,"_diffcgis_Mehist.txt"),sep = "\t",header = T)
  density.t <- data
  for (x in 1:nrow(data)){
    density.t[x,2:11] <- density.t[x,2:11]/sum(density.t[x,2:11])
  }
  hyper.tumor.hist <- apply(density.t[hyper_rg,2:11],2,mean_rmna)
  hyper.tumor.hist <- hyper.tumor.hist/sum(hyper.tumor.hist)
  hypo.tumor.hist <- apply(density.t[hypo_rg,2:11],2,mean_rmna)
  hypo.tumor.hist <- hypo.tumor.hist/sum(hypo.tumor.hist)
  
  tumor.hist <- data.frame(hyper = hyper.tumor.hist,hypo = hypo.tumor.hist,x = seq(1,10),row.names = 1:10)
  library(ggplot2)
  pdf(paste0("I:/Projects/P7.Methy-PanCancer/pics/HCC_WGBS/HCC_",idx,".Mehist.pdf"),height = 3,width = 6)
  text_df <- data.frame(x = c(5,5),y = c(0.3,-0.3),label = c(paste0("Tumor-",sampleID2[idx,1]),"Normal"))
  g1 <- ggplot(tumor.hist,aes(x =x, y = hyper)) + geom_col(width = 1,fill = "#FF8000",color = "black") + theme_bw() + ggtitle("Tumor Hyper")
  g2 <- ggplot(normal.hist,aes(x =x, y = hyper)) + geom_col(width = 1,fill = "#C0C0C0",color = "black") + theme_bw() + ggtitle("Normal Hyper")
    
  g3 <- ggplot(tumor.hist,aes(x =x, y = hypo)) + geom_col(width = 1,fill = "#6A5ACD",color = "black") + theme_bw() + ggtitle("Tumor Hypo")
  g4 <- ggplot(normal.hist,aes(x =x, y = hypo)) + geom_col(width = 1,fill = "#C0C0C0",color = "black") + theme_bw() + ggtitle("Normal Hypo")
  print(g1)
  print(g2)
  print(g3)
  print(g4)
  dev.off()
}

#筛选最多reads数目的CGIs
counts <- NULL
for (idx in 1:nrow(sampleID1)){
  data <- read.table(paste0(path,sampleID1[idx,1],"_diffcgis_Mehist.txt"),sep = "\t",header = T)
  counts <- cbind(counts,apply(data[,2:11],1,sum))
}
tumor.list <- c("SRR15860722","SRR15860734","SRR15860750")
for (idx in tumor.list){
  data <- read.table(paste0(path,idx,"_diffcgis_Mehist.txt"),sep = "\t",header = T)
  counts <- cbind(counts,apply(data[,2:11],1,sum))
}
counts <- as.data.frame(counts)
rownames(counts) <- data$BedID
colnames(counts) <- c(sampleID1$V1,tumor.list)
f <- apply(counts[,1:(ncol(counts)-3)],2,sum)
order1 <- order(f,decreasing = T)
counts_f <- counts[,c(8,ncol(counts)-2,ncol(counts)-1,ncol(counts))]
cgiinfo <- read.table("I:/Projects/P7.Methy-PanCancer/data_analysis/A5.HCC_WGBS/A5.HCC.WGBS.diffcgi.info.txt",sep = "\t",header = T)
rownames(cgiinfo) <- paste0(cgiinfo$chrom,"_",cgiinfo$start,"_",cgiinfo$end)
final.info <- cbind(counts_f,cgiinfo[rownames(counts_f),"NumCpG"],cgiinfo[rownames(counts_f),"label"])
colnames(final.info)[5:6] <- c("Numcpg","Label")
final.info$Sp1 <- final.info[,1]/final.info[,5]
final.info$Sp2 <- final.info[,2]/final.info[,5]
final.info$Sp3 <- final.info[,3]/final.info[,5]
final.info$Sp4 <- final.info[,4]/final.info[,5]
final.info$total <- apply(final.info[,7:10],1,sum)
write.table(final.info,'I:/Projects/P7.Methy-PanCancer/data_analysis/A5.HCC_WGBS/A5.HCC_Examples.ReadsCounts.txt',
            sep = "\t",row.names = T,col.names = T,quote = F)
## compare tumor clusters
cluster1 <- c("SRR15860722","SRR15860709","SRR15860756","SRR15860749","SRR15860735")
cluster2 <- c("SRR15860745","SRR15860734","SRR15860747","SRR15860748","SRR15860742","SRR15860739")
cluster3 <- c("SRR15860790","SRR15860753","SRR15860750","SRR15860741","SRR15860731")

spinfo <- data.frame(spid = c(cluster1,cluster2,cluster3),type = c(rep("Hyper",length(cluster1)),
                                                                   rep("Intermediate",length(cluster2)),
                                                                   rep("Hypo",length(cluster3)))
                     )
vals <- NULL

for (idx in 1:nrow(spinfo)){
  data <- read.table(paste0(path,spinfo[idx,1],"_diffcgis_Mehist.txt"),sep = "\t",header = T)
  density.t <- data
  for (x in 1:nrow(data)){
    density.t[x,2:11] <- density.t[x,2:11]/sum(density.t[x,2:11])
  }
  hyper.tumor.hist <- apply(density.t[hyper_rg,2:11],2,mean_rmna)
  hyper.tumor.hist <- hyper.tumor.hist/sum(hyper.tumor.hist)
  hypo.tumor.hist <- apply(density.t[hypo_rg,2:11],2,mean_rmna)
  hypo.tumor.hist <- hypo.tumor.hist/sum(hypo.tumor.hist)
  
  tumor.hist <- data.frame(hyper = hyper.tumor.hist,hypo = hypo.tumor.hist,x = seq(1,10),row.names = 1:10)
  
  vals <- rbind(vals,c(sum(tumor.hist[6:9,1])/tumor.hist[10,1],sum(tumor.hist[2:5,2])/tumor.hist[1,2]))
}
out <- cbind(spinfo,vals)
out <- as.data.frame(out)
out[,3] <- as.numeric(out[,3])
out[,4] <- as.numeric(out[,4])
colnames(out) <- c("ID","Type","hyper.Ratio","hypo.Ratio")

library(ggplot2)
pdf("I:/Projects/P7.Methy-PanCancer/pics/HCC_WGBS/HCC_RIN_boxplot.pdf",height = 6,width = 3)
out1 <- out[(out$Type == "Hyper") |(out$Type == "Intermediate"),]
ga <- ggplot(out1,aes(x = Type,y=hyper.Ratio)) +geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) + theme_bw()
out2 <- out[(out$Type == "Hypo") |(out$Type == "Intermediate"),]
gb <- ggplot(out2,aes(x = Type,y=hypo.Ratio)) +geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2) + theme_bw()
print(ga)
print(gb)
dev.off()
wilcox.test(out1$hyper.Ratio[out1$Type == "Hyper"],out1$hyper.Ratio[out1$Type == "Intermediate"])
wilcox.test(out2$hypo.Ratio[out2$Type == "Hypo"],out2$hypo.Ratio[out2$Type == "Intermediate"])
