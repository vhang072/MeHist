# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 09:46:40 2024

@author: ABC
"""
import os
import pandas as pd 


datalist = ["Y:/4.basic_data/TCGA_PancanAtlas/methylation/450_allprobes/Primary Blood Derived Cancer - Peripheral Blood/",\
            "Y:/4.basic_data/TCGA_PancanAtlas/methylation/450_allprobes/Primary Tumor/",\
            "Y:/4.basic_data/TCGA_PancanAtlas/methylation/450_allprobes/Solid Tissue Normal/",\
            "Y:/4.basic_data/TCGA_PancanAtlas/methylation/450_allprobes_cell_lines/",\
            "Y:/4.basic_data/TCGA_PancanAtlas/methylation/450_GEO_Normal/_PureFile/"]

outlist = ["Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/LAML/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/Primary/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/SolidNormal/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/CellLines/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/GEONormal/"]
    
probeused = pd.read_csv("Y:/4.basic_data/TCGA_PancanAtlas/methylation/info_output/probe_used.txt",\
                        sep = "\t",header = 0)
probelist = probeused.loc[probeused.loc[:,"Ifincluded"]==True,"cgid"].tolist()

for x in range(len(datalist)):
    allfiles = os.listdir(datalist[x])
    for file in allfiles:
        data = pd.read_csv(datalist[x]+file,sep = "\t",header = 0,na_values = ["nan","NA","NaN","N/A"])
        data.index = data.iloc[:,0].tolist()
        outdata = data.loc[probelist,:]
        outdata.to_csv(outlist[x]+file,sep = "\t",index = False)
