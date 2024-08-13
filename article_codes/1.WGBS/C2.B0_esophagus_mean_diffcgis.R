diffcgis <- read.table("I:/Projects/P7.Methy-PanCancer/data_analysis/B0.esophagus_WGBS/B0.esophagus.WGBS.diffcgi.info.txt",
                       sep = "\t",header = T)
cgiid <- paste0(diffcgis$chrom,"_",diffcgis$start,"_",diffcgis$end)
rownames(diffcgis) <- cgiid
mefile <- read.table("I:/Projects/P7.Methy-PanCancer/data_analysis/B0.esophagus_WGBS/B0.esophagus_cgis_merged.bed",
                     sep = "\t",header = T)
cgiid2 <- paste0(mefile$chrom,"_",mefile$start,"_",mefile$end)
rownames(mefile) <- cgiid2

hypercgis <- diffcgis[diffcgis$label == "Hyper",]
mehyper <- mefile[rownames(hypercgis),]
hypocgis <- diffcgis[diffcgis$label == "Hypo",]
mehypo <- mefile[rownames(hypocgis),]

meanhyper <- apply(mehyper[,47:ncol(mehyper)],2,mean) - mean(hypercgis$Mean_N)
meanhypo <- mean(hypocgis$Mean_N) - apply(mehypo[,47:ncol(mehypo)],2,mean) 
sampleme <- data.frame(meanhyper = meanhyper,meanhypo = meanhypo,name =  colnames(mefile)[47:ncol(mefile)],
                       row.names = colnames(mefile)[47:ncol(mefile)])
write.table(sampleme,"I:/Projects/P7.Methy-PanCancer/data_analysis/B0.esophagus_WGBS/B0.esophagus_MeanMe_persample.txt",
            sep = "\t",row.names = T,quote = F)

cor1 <- cor(t(mehyper[,47:ncol(mehyper)]),meanhyper)
hist(cor1)
cor2 <- cor(t(mehypo[,47:ncol(mehypo)]),-1*meanhypo)
hist(cor2)

spinfo <- read.table("I:/Projects/P7.Methy-PanCancer/data_analysis/B0.esophagus_WGBS/B0.sampleid.txt",
                     sep = "\t",header = F)
library(stringr)
infosp <- matrix(unlist(str_split(spinfo[,2],"CC")),ncol = 2, byrow = T)
spinfo$id <- paste0(spinfo[,3],"_",infosp[,2])
rownames(spinfo) <- spinfo$V1
sampleme$id <- spinfo[sampleme$name,"id"]


library(ggplot2)
pdf("I:/Projects/P7.Methy-PanCancer/pics/Esophagus_WGBS/Esophagus_hyperVShypo.pdf",height = 6,width = 6)
ggplot(sampleme,aes(x = meanhypo,y = meanhyper,label = id)) + geom_point() + geom_text(aes(y = meanhyper + 0.02))+
  xlim(-0.05,0.4) + theme_bw()
dev.off()
cor.test(sampleme$meanhyper,sampleme$meanhypo)
