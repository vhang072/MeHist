# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 10:46:01 2024

@author: ABC
"""

import pandas as pd
import os

def cgprobeIntoL1(data_chrom,l1_chrom):
    idx1 = 0
    idx2 = 0
    num1 = data_chrom.shape[0]
    num2 = l1_chrom.shape[0]
    probewithinl1 = []
    probel1name = []
    while idx1 < num1 and idx2 <num2:
        l1_name = l1_chrom.iloc[idx2,0] + "_" + str(l1_chrom.iloc[idx2,1]) + "_" + \
                            str(l1_chrom.iloc[idx2,2])
        if data_chrom.iloc[idx1,1] < l1_chrom.iloc[idx2,1] + 1:
            idx1 += 1
        elif data_chrom.iloc[idx1,1] >= l1_chrom.iloc[idx2,1] + 1 and \
            data_chrom.iloc[idx1,1] <= l1_chrom.iloc[idx2,2]:
            probewithinl1.append(idx1)
            probel1name.append(l1_name)
            idx1 += 1
        elif data_chrom.iloc[idx1,1] > l1_chrom.iloc[idx2,2]:
            idx2 += 1
    selectprobe = data_chrom.iloc[probewithinl1,:]
    selectprobe.index = range(selectprobe.shape[0])
    l1probes  = pd.concat([pd.DataFrame({"cgi":probel1name}),selectprobe],axis = 1)
    l1probes = l1probes.rename(columns = {l1probes.columns[3]:"probeid"})
    return l1probes
                
            


datalist = ["Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/LAML/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/Primary/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/SolidNormal/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/CellLines/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/GEONormal/"]

outlist = ["Y:/4.basic_data/TCGA_PancanAtlas/methylation_L1/LAML/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_L1/Primary/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_L1/SolidNormal/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_L1/CellLines/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_L1/GEONormal/"]    
    
l1 = pd.read_csv("I:/Projects/P7.Methy-PanCancer/genome_reference/Repeats/hg38_L1_Mainchrom.txt",\
                   header = None, sep = "\t")
chroms = ['chr' + str(x) for x in range(1,23)]

for pathnum in range(len(datalist)):
    path = datalist[pathnum]
    files = os.listdir(path)
    for file in files:
        data = pd.read_csv(path + file,sep = "\t",header = 0)
        print("Processing File:" + file)
        #data = pd.read_csv(datalist[0] + "TCGA-LAML.txt",sep = "\t",header = 0)
        data_chrom = []
        data_pos = []
        for x in range(data.shape[0]):
            data_chrom.append(data.iloc[x,0].split("_")[0])
            data_pos.append(int(data.iloc[x,0].split("_")[1]))
        probe_info = pd.DataFrame({"chrom":data_chrom,"pos":data_pos},index = data.iloc[:,0].tolist())
        data.index = data.iloc[:,0].tolist()
        data_combined = pd.merge(probe_info,data,how = "outer",left_index=True,right_index=True)
        first_flag = 1
        for chrom in chroms:
            print(chrom)
            data_chrom0 = data_combined.loc[data_combined.iloc[:,0]==chrom,:]
            l1_chrom0 = l1.loc[l1.iloc[:,0] == chrom,:]
            data_chrom = data_chrom0.sort_values(by="pos",inplace = False)
            l1_chrom = l1_chrom0.sort_values(by=1,inplace = False)
            l1probes = cgprobeIntoL1(data_chrom,l1_chrom)
            if first_flag == 1:
                l1probe_total = l1probes
                first_flag = 0
            else:
                l1probe_total = pd.concat([l1probe_total,l1probes],axis = 0)
        #file1:l1 probe data
        l1probe_total.to_csv(outlist[pathnum] + "_L1Probe_" + file,sep = "\t",na_rep = "NA",index = False)    






    