# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 08:25:42 2024

@author: ABC
"""

import pandas as pd
import os

def cgprobeIntoCGI(data_chrom,cgi_chrom):
    idx1 = 0
    idx2 = 0
    num1 = data_chrom.shape[0]
    num2 = cgi_chrom.shape[0]
    probewithincgi = []
    probecginame = []
    cgiprobenum = []
    cgime = pd.DataFrame(columns = data_chrom.columns)
    cgisres = pd.DataFrame(columns = ['cgi'] + data_chrom.columns[3:].tolist())
    while idx1 < num1 and idx2 <num2:
        cgi_name = cgi_chrom.iloc[idx2,0] + "_" + str(cgi_chrom.iloc[idx2,1]) + "_" + \
                            str(cgi_chrom.iloc[idx2,2])
        if data_chrom.iloc[idx1,1] < cgi_chrom.iloc[idx2,1] + 1:
            idx1 += 1
        elif data_chrom.iloc[idx1,1] >= cgi_chrom.iloc[idx2,1] + 1 and \
            data_chrom.iloc[idx1,1] <= cgi_chrom.iloc[idx2,2]:
            probewithincgi.append(idx1)
            probecginame.append(cgi_name)
            cgime = pd.concat([cgime,data_chrom.iloc[idx1,:].to_frame().T],axis =0 )
            idx1 += 1
        elif data_chrom.iloc[idx1,1] > cgi_chrom.iloc[idx2,2]:
            if cgime.shape[0] >= 2:
                cgimerge = cgime.iloc[:,3:].mean(axis = 0,skipna = False)
                cgivalue = cgimerge.to_frame().T
                cgivalue.index = [0]
                cgidf = pd.concat([pd.DataFrame({"cgi":[cgi_name]}),cgivalue],axis = 1)
                cgisres = pd.concat([cgisres,cgidf],axis = 0)
            elif cgime.shape[0] == 1:
                cgivalue = cgime.iloc[:,3:]
                cgivalue.index = [0]
                cgidf = pd.concat([pd.DataFrame({"cgi":[cgi_name]}),cgivalue],axis = 1)
                cgisres = pd.concat([cgisres,cgidf],axis = 0)
            cgiprobenum.append(cgime.shape[0])
            cgime = pd.DataFrame(columns = data_chrom.columns)
            idx2 += 1
    while idx2 < num2:
        if cgime.shape[0] >= 2:
            cgimerge = cgime.iloc[:,3:].mean(axis = 0,skipna = False,numeric_only = True)
            cgivalue = cgimerge.to_frame().T
            cgivalue.index = [0]
            cgidf = pd.concat([pd.DataFrame({"cgi":[cgi_name]}),cgivalue],axis = 1)
            cgisres = pd.concat([cgisres,cgidf],axis = 0)
        elif cgime.shape[0] == 1:
            cgivalue = cgime.iloc[:,3:]
            cgivalue.index = [0]
            cgidf = pd.concat([pd.DataFrame({"cgi":[cgi_name]}),cgivalue],axis = 1)
            cgisres = pd.concat([cgisres,cgidf],axis = 0)
        cgiprobenum.append(cgime.shape[0])
        cgime = pd.DataFrame(columns = data_chrom.columns)
        idx2 += 1
    while idx1 < num1:
        idx1 += 1
    selectprobe = data_chrom.iloc[probewithincgi,:]
    selectprobe.index = range(selectprobe.shape[0])
    cgiprobes  = pd.concat([pd.DataFrame({"cgi":probecginame}),selectprobe],axis = 1)
    cgiprobes = cgiprobes.rename(columns = {cgiprobes.columns[3]:"probeid"})
    cgi_chrom.index = range(len(cgi_chrom))
    cgi_chrom_improved = pd.concat([cgi_chrom,pd.DataFrame({"probenumin350k":cgiprobenum})],axis = 1)
    return cgiprobes,cgi_chrom_improved,cgisres
                
            


datalist = ["Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/LAML/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/Primary/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/SolidNormal/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/CellLines/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_selected350k/GEONormal/"]

outlist = ["Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/LAML/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/Primary/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/SolidNormal/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/CellLines/",\
           "Y:/4.basic_data/TCGA_PancanAtlas/methylation_cgiprobe/GEONormal/"]    
    
cgis = pd.read_csv("I:/Projects/P7.Methy-PanCancer/genome_reference/CpG_islands/hg38_regulation_CpG_islands_ucsctable_onlyMainChroms_sort.bed",\
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
            cgi_chrom0 = cgis.loc[cgis.iloc[:,0] == chrom,:]
            data_chrom = data_chrom0.sort_values(by="pos",inplace = False)
            cgi_chrom = cgi_chrom0.sort_values(by=1,inplace = False)
            cgiprobes,cgi_chrom_improved,cgisres = cgprobeIntoCGI(data_chrom,cgi_chrom)
            if first_flag == 1:
                cgiprobe_total = cgiprobes
                cgiinfo_total = cgi_chrom_improved
                cgime_total = cgisres
                first_flag = 0
            else:
                cgiprobe_total = pd.concat([cgiprobe_total,cgiprobes],axis = 0)
                cgiinfo_total = pd.concat([cgiinfo_total,cgi_chrom_improved],axis = 0)
                cgime_total = pd.concat([cgime_total,cgisres],axis = 0)
        #file1:cgi probe data
        cgiprobe_total.to_csv(outlist[pathnum] + "_cgiProbe_" + file,sep = "\t",na_rep = "NA",index = False)   
        #file2: cgi file info
        cgiinfo_total2 = cgiinfo_total.astype({1:'int32',2:'int32','probenumin350k':'int32'})
        cgiinfo_total2.to_csv(outlist[pathnum] + "_cgiInfo_" + file,sep = "\t",na_rep = "NA",index = False)  
        #file3: cgi Me
        dt = cgime_total.dtypes
        dt[1:] = "float"
        cgime_total2 = cgime_total.astype(dt)
        cgime_total2.to_csv(outlist[pathnum] + "_cgiMergedMe_" + file,sep = "\t",na_rep = "NA",index = False,\
                           float_format = '%.3f')  






    