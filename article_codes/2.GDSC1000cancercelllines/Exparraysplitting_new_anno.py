# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 17:39:32 2024

@author: ABC
"""
def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

import pandas as pd 
import os 

anno = pd.read_csv(r"Y:\4.basic_data\TCGA_PancanAtlas\cgi_methy_vs_exp\7.cancer_cell_line\data_The landscape of pharmacogenomic interactions in human cancer\cell_line_annotation_frompaperresult.txt",\
                   header = 0,sep = "\t")
    
exp = pd.read_csv(r"Y:\4.basic_data\TCGA_PancanAtlas\cgi_methy_vs_exp\7.cancer_cell_line\data_The landscape of pharmacogenomic interactions in human cancer\Cell_line_RMA_proc_basalExp.txt",\
                  header = 0,sep = "\t")

path = "Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/CellLines/"
pathout = "Y:/4.basic_data/TCGA_PancanAtlas/RNA_CancerCellLines_rmGeneTitles_anotheranno/"

files = os.listdir(path)

for file in files:
    data = pd.read_csv(path+file,sep = "\t",header = 0)
    cl1 = []
    cl2 = []
    prefix = ["GENE_SYMBOLS"]#,"GENE_title"
    for idx in range(1,data.shape[1]):
        now = data.columns[idx].split("_AVG.Beta")[0]
        if sum(anno.iloc[:,0] == now) == 1:
            cl1.append(now)
            target = "DATA." + str(anno.loc[anno.iloc[:,0] == now,'COSMIC_identifier'].iloc[0])
            cl2.append(target)
    
    ck1 = []
    ck2 = []
    for x in range(len(cl1)):
        if cl2[x] in exp.columns.tolist():
            ck2.append(cl2[x])
            ck1.append(cl1[x])
        else:
            pass
    out = exp[prefix + ck2]
    out.columns = prefix + ck1
    out.to_csv(pathout+"Exp_"+file[19:],sep = "\t",index = False,na_rep = "NA")
        
