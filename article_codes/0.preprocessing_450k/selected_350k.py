# -*- coding: utf-8 -*-
"""
Created on Mon Apr  1 10:25:50 2024

@author: ABC
"""

#select useable probes
import os
import pandas as pd 
import matplotlib
import numpy as np

datalist = ["Y:/4.basic_data/TCGA_PancanAtlas/methylation/450_allprobes/Primary Blood Derived Cancer - Peripheral Blood/",\
            "Y:/4.basic_data/TCGA_PancanAtlas/methylation/450_allprobes/Primary Tumor/",\
            "Y:/4.basic_data/TCGA_PancanAtlas/methylation/450_allprobes/Solid Tissue Normal/",\
            "Y:/4.basic_data/TCGA_PancanAtlas/methylation/450_allprobes_cell_lines/",\
            "Y:/4.basic_data/TCGA_PancanAtlas/methylation/450_GEO_Normal/_PureFile/"]
first_flag = 1
total_samples = []
for x in range(len(datalist)):
    allfiles = os.listdir(datalist[x])
    for file in allfiles:
        data = pd.read_csv(datalist[x]+file,sep = "\t",header = 0,na_values = ["nan","NA","NaN","N/A"])
        datainfo = data.isna().sum(axis = 1)
        outdf = pd.DataFrame(dict(cgid=data.iloc[:,0].tolist(),\
                          file = datainfo.tolist()))
        outdf.columns = ['cgid',file]
        if first_flag == 1:
            totaldf = outdf
            first_flag = 0
        else:
            totaldf = pd.merge(totaldf,outdf,how = "outer",on = "cgid")
        total_samples.append(data.shape[1]-1)
nanumsprobe = totaldf.iloc[:,1:].sum(axis = 1)
nanumsdataset = totaldf.iloc[:,1:].isna().sum(axis = 1)
#matplotlib.pyplot.hist(nanumsprobe,bins = 100,range = [0,100])
#sum(np.array(nanumsprobe<=5) & np.array(nanumsdataset <= 0))
probeinfo = pd.DataFrame(dict(cgid = totaldf.iloc[:,0], \
                              NAnumsample = nanumsprobe,\
                              NAnumdataset = nanumsdataset))
chrom = []
for x in range(probeinfo.shape[0]):
    chrom.append(probeinfo.iloc[x,0].split("_")[0])
probeinfo.loc[:,"chrom"] = chrom
probeinfo.loc[:,"Ifincluded"] = (probeinfo.loc[:,"chrom"] != "chrX") & \
    (probeinfo.loc[:,"chrom"] != "chrY") & \
        (probeinfo.loc[:,"NAnumsample"] <=5) & \
           (probeinfo.loc[:,"NAnumdataset"] <=0) 
totaldf.to_csv("Y:/4.basic_data/TCGA_PancanAtlas/methylation/info_output/probe_NA_times.txt",\
               sep = "\t",index = False)

probeinfo.to_csv("Y:/4.basic_data/TCGA_PancanAtlas/methylation/info_output/probe_used.txt",\
               sep = "\t",index = False)
    
#checkedprobe = totaldf.loc[nanumsprobe>=1,:]

sum(probeinfo.loc[:,"Ifincluded"] == True)