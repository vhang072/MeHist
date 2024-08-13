data <- read.table("I:/Projects/P7.Methy-PanCancer/data_analysis/B0.esophagus_WGBS/B0.esophagus_cgis_merged.bed",
                   sep = "\t",header = T)
x <- is.na.data.frame(data[,47:88])
nasrow <- apply(x,1,sum)
#hist(nasrow,breaks = 100)
data_f <- data[nasrow==0,]

out <- NULL
thr_me <- c(0.1, 0.15, 0.2, 0.25, 0.3)#
for (idx in 1:nrow(data_f)){
  me.normal <- as.numeric(data_f[idx,5:46])
  me.normal <- me.normal[is.na(me.normal)== F]
  me.tumor <- as.numeric(data_f[idx,47:88])
  mean.normal <- mean(me.normal)
  mean.tumor <- mean(me.tumor)
  
  test1 <- wilcox.test(me.normal, me.tumor,exact = F)
  p1 <- test1$p.value
  abe_frac_t <- NULL#
  abe_frac_n <- NULL#
  if (mean.tumor > mean.normal){
    for (thr in thr_me){
      abe_frac_t <- c(abe_frac_t, sum(me.tumor > (mean.normal+thr))/length(me.tumor))#
      abe_frac_n <- c(abe_frac_n, sum(me.normal > (mean.normal+thr))/length(me.normal))#
    }
  } else {
    for (thr in thr_me){
      abe_frac_t <- c(abe_frac_t, sum(me.tumor < (mean.normal-thr))/length(me.tumor))#
      abe_frac_n <- c(abe_frac_n, sum(me.normal < (mean.normal-thr))/length(me.normal))#
    }
  }
  out <- rbind(out,c(data_f[idx,1:4],length(me.tumor),length(me.normal),mean.tumor,mean.normal,p1,abe_frac_t,abe_frac_n))
  
}
out <- as.data.frame(out)
thr_names <- c(paste0("AberrantFrac_T_",thr_me),paste0("AberrantFrac_N_",thr_me))#
colnames(out) <- c(colnames(data)[1:4],"Num_T","Num_N","Mean_T","Mean_N","P.wilcox",thr_names)#

for (x in 2:ncol(out)){
  out[,x] <- as.numeric(out[,x])
}
out[,1] <- as.character(out[,1])

hyper_f1 <- (out$Mean_T > out$Mean_N) & (out$AberrantFrac_T_0.2 > 0.25) & (out$Mean_N < 0.3) & (out$AberrantFrac_N_0.2 < 0.05)
hypo_f2 <- (out$Mean_T < out$Mean_N) & (out$AberrantFrac_T_0.2 > 0.25) & (out$Mean_N > 0.7) & (out$AberrantFrac_N_0.2 < 0.05)
sum(hyper_f1)
sum(hypo_f2)

out$label <- "Stable"
out$label[hyper_f1] <- "Hyper"
out$label[hypo_f2] <- "Hypo"

write.table(out,"I:/Projects/P7.Methy-PanCancer/data_analysis/B0.esophagus_WGBS/B0.esophagus.WGBS.diffcgi.info.txt",
            sep = "\t",row.names = F, col.names = T,quote = F)