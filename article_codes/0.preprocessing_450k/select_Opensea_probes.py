# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 15:03:57 2024

@author: ABC
"""

import pandas as pd
import os

datalist = ["Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/LAML/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/Primary/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/SolidNormal/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/CellLines/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/GEONormal/"]

outlist = ["Y:/4.basic_data/TCGA_PancanAtlas/methylation_opensea/LAML/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_opensea/Primary/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_opensea/SolidNormal/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_opensea/CellLines/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_opensea/GEONormal/"]  
    
probeinfo = pd.read_csv(r"Y:\4.basic_data\TCGA_PancanAtlas\methylation\jhu-usc.edu_BRCA.HumanMethylation450.8.lvl-3.TCGA-E2-A1B4-01A-11D-A12R-05.gdc_hg38.txt",
                        sep = "\t",header = 0)
probeid = []
for x in range(probeinfo.shape[0]):
    probeid.append(probeinfo.iloc[x,2] + "_" + str(probeinfo.iloc[x,3]) + "_" + probeinfo.iloc[x,0])
    
probeinfo.index = probeid    

probesea = probeinfo.loc[probeinfo["Feature_Type"] == ".",:]


for pathnum in range(len(datalist)):
    path = datalist[pathnum]
    files = os.listdir(path)
    for file in files:
        data = pd.read_csv(path + file,sep = "\t",header = 0)
        data.index = data.iloc[:,0]
        print("Processing File:" + file)
        probeused = set(data.iloc[:,0].tolist()).intersection(set(probesea.index.tolist()))
        sortprobe = [x for x in probeused]
        sortprobe.sort()
        data_select = data.loc[sortprobe,:]
        data_select2 = data_select.rename(columns = {data_select.columns[0]:"probeid"})
        data_select2.to_csv(outlist[pathnum] + "_OpenseaProbe_" + file,sep = "\t",na_rep = "NA",index = False)   
        
