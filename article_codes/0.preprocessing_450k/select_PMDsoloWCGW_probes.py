# -*- coding: utf-8 -*-
"""
Created on Wed Apr 10 09:07:20 2024

@author: ABC
"""

import pandas as pd
import os

def cgprobeIntosolowcgw(data_chrom,solowcgw_chrom):
    idx1 = 0
    idx2 = 0
    num1 = data_chrom.shape[0]
    num2 = solowcgw_chrom.shape[0]
    probeselected = []
    while idx1 < num1 and idx2 <num2:
        if data_chrom.iloc[idx1,1] < solowcgw_chrom.iloc[idx2,1] + 1:
            idx1 += 1
        elif data_chrom.iloc[idx1,1] == solowcgw_chrom.iloc[idx2,1] + 1:
            probeselected.append(idx1)
            idx1 += 1
        elif data_chrom.iloc[idx1,1] > solowcgw_chrom.iloc[idx2,1] + 1:
            idx2 += 1
    selectedprobe = data_chrom.iloc[probeselected,:]
    selectedprobe.index = range(selectedprobe.shape[0])
    solowcgwprobes = selectedprobe.rename(columns = {selectedprobe.columns[2]:"probeid"})
    return solowcgwprobes
                
            


datalist = ["Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/LAML/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/Primary/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/SolidNormal/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/CellLines/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/GEONormal/"]

outlist = ["Y:/4.basic_data/TCGA_PancanAtlas/methylation_soloWCGW/LAML/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_soloWCGW/Primary/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_soloWCGW/SolidNormal/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_soloWCGW/CellLines/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_soloWCGW/GEONormal/"]    
    
solowcgw = pd.read_csv("I:/Projects/P7.Methy-PanCancer/genome_reference/PMD_soloWCGW/solo_WCGW_inCommonPMDs_hg38.bed",\
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
            solowcgw_chrom0 = solowcgw.loc[solowcgw.iloc[:,0] == chrom,:]
            data_chrom = data_chrom0.sort_values(by="pos",inplace = False)
            solowcgw_chrom = solowcgw_chrom0.sort_values(by=1,inplace = False)
            
            solowcgwprobes = cgprobeIntosolowcgw(data_chrom,solowcgw_chrom)
            if first_flag == 1:
                solowcgwprobe_total = solowcgwprobes
                first_flag = 0
            else:
                solowcgwprobe_total = pd.concat([solowcgwprobe_total,solowcgwprobes],axis = 0)
        #file1:soloWCGW probe data
        solowcgwprobe_total.to_csv(outlist[pathnum] + "_PMDsoloWCGWProbe_" + file,sep = "\t",na_rep = "NA",index = False)   

    