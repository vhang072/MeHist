data <- read.table("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_TCGA/__DiffCGIprobeCorSummary.txt",
                   sep = "\t",header = T)
show <- data[c(2:13,1,14:27),2:3]
rownames(show) <- data$cancer_types[c(2:13,1,14:27)]
colnames(show) <- c("r.Hyper","r.Hypo")
show$r.Hypo <- -1*show$r.Hypo

library(pheatmap)
library(paletteer)
pdf("I:/Projects/P7.Methy-PanCancer/pics/CorAcrossDiffProbes/TCGA_corindiffprobes.pdf",height = 1.2,width = 9)
p <- pheatmap(t(show), cluster_rows =F ,cluster_cols =F,color = paletteer_d("RColorBrewer::YlOrRd") ,breaks = seq(0,0.8,0.1))
print(p)
dev.off()


data <- read.table("Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Diff_CCLE/__DiffCGIprobeCorSummary.txt",
                   sep = "\t",header = T)
rownames(data) <- data$cancer_types
used <- c("aero_digestive_tract_head and neck",
                         "aero_digestive_tract_oesophagus",
                         "blood_B_cell_lymphoma",
                         "blood_acute_myeloid_leukaemia",
                         "bone_ewings_sarcoma",
                         "breast_breast",
                         "digestive_system_large_intestine",
                         "digestive_system_stomach",
                         "kidney_kidney",
                         "lung_lung_NSCLC_adenocarcinoma",
                         "lung_lung_small_cell_carcinoma",
                         "lung_mesothelioma",
                         "nervous_system_glioma",
                         "nervous_system_neuroblastoma",
                         "pancreas_pancreas",
                         "skin_melanoma",
                         "urogenital_system_ovary",
                         "urogenital_system_bladder")
show <- data[used,2:3]
rownames(show) <- used
colnames(show) <- c("r.Hyper","r.Hypo")
show$r.Hypo <- -1*show$r.Hypo

library(pheatmap)
library(paletteer)
pdf("I:/Projects/P7.Methy-PanCancer/pics/CorAcrossDiffProbes/CCLE_corindiffprobes.pdf",height = 2.7,width = 9)
p <- pheatmap(t(show), cluster_rows =F ,cluster_cols =F,color = paletteer_d("RColorBrewer::YlOrRd") ,breaks = seq(0,0.8,0.1),
              angle_col = 45)
print(p)
dev.off()


