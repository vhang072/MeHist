diffcgis <- read.table("I:/Projects/P7.Methy-PanCancer/data_analysis/A5.HCC_WGBS/A5.HCC.WGBS.diffcgi.info.txt",
                       sep = "\t",header = T)
cgiid <- paste0(diffcgis$chrom,"_",diffcgis$start,"_",diffcgis$end)
rownames(diffcgis) <- cgiid
mefile <- read.table("I:/Projects/P7.Methy-PanCancer/data_analysis/A5.HCC_WGBS/A5.HCC_cgis_merged.bed",
                     sep = "\t",header = T)
cgiid2 <- paste0(mefile$chrom,"_",mefile$start,"_",mefile$end)
rownames(mefile) <- cgiid2

hypercgis <- diffcgis[diffcgis$label == "Hyper",]
mehyper <- mefile[rownames(hypercgis),]
hypocgis <- diffcgis[diffcgis$label == "Hypo",]
mehypo <- mefile[rownames(hypocgis),]

meanhyper <- apply(mehyper[,38:ncol(mehyper)],2,mean) - mean(hypercgis$Mean_N)
meanhypo <- mean(hypocgis$Mean_N) - apply(mehypo[,38:ncol(mehypo)],2,mean) 
sampleme <- data.frame(meanhyper = meanhyper,meanhypo = meanhypo, name = colnames(mefile)[38:ncol(mefile)],
                       row.names = colnames(mefile)[38:ncol(mefile)])


write.table(sampleme,"I:/Projects/P7.Methy-PanCancer/data_analysis/B0.esophagus_WGBS/A5.HCC_MeanMe_persample.txt",
            sep = "\t",row.names = T,quote = F)

cor1 <- cor(t(mehyper[,38:ncol(mehyper)]),meanhyper)
hist(cor1)
cor2 <- cor(t(mehypo[,38:ncol(mehypo)]), -1*meanhypo)
hist(cor2)


library(ggplot2)
library(stringr)
srr <- read.table("I:/Projects/P7.Methy-PanCancer/data_analysis/A5.HCC_WGBS/A5.HCC.SRR_name.txt",
                  sep = "\t",header = T)
sid <- matrix(unlist(str_split(srr$Sample.Name,"_")),ncol = 3,byrow = T)
sidcat <- paste0(sid[,2],"_",sid[,3])
srr$id <- sidcat
rownames(srr) <- srr$Run
sampleme$id <- srr[sampleme$name,"id"]

pdf("I:/Projects/P7.Methy-PanCancer/pics/HCC_WGBS/HCC_hyperVShypo.pdf",height = 6,width = 6)
ggplot(sampleme,aes(x = meanhypo,y = meanhyper,label = id)) + geom_point() + geom_text(aes(y = meanhyper + 0.01))+
  xlim(-0.05,0.7) + theme_bw()
dev.off()

