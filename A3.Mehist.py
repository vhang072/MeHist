# -*- coding: utf-8 -*-
"""
"""

import h5py
import getopt,sys
import numpy as np
import pandas as pd

options = "i:o:b:"
long_options = ["input=","output=","bed="]
#path_in = "I:/tmp_download/SRR9070182_cgi_MeMatrix.hdf5"
#bedfile = "I:/Projects/P7.Methy-PanCancer/data_analysis/A5.HCC_WGBS/hyper_hypo_cgis_A5HCC.bed"
arguments, values = getopt.getopt(sys.argv[1:], options, long_options)
for currentArgument, currentValue in arguments:
    if currentArgument in ("-i", "--input"):
        path_in = currentValue
    elif currentArgument in ("-o", "--output"):
        path_out = currentValue
    elif currentArgument in ("-b", "--bed"):
        bedfile = currentValue    
        
h = h5py.File(path_in,"r")
carebeds = pd.read_csv(bedfile,sep = "\t",header = None)
bins = np.arange(0,1.1,0.1,dtype = float)
bins[0] = -0.01
cols = ["BedID"] + [str(x) + "tenth" for x in range(10)]
metricdata = pd.DataFrame([],columns= cols)
for idx in range(carebeds.shape[0]):
    uid = carebeds.iloc[idx,0] + "_" + str(carebeds.iloc[idx,1]) + "_" + str(carebeds.iloc[idx,2])
    memat = h[uid][...]
    tmpbed = np.zeros((1,10),dtype = int)
    for x in range(2,memat.shape[0]):
        mes = sum(memat[x,] == 1)
        unmes = sum(memat[x,] == -1)
        if mes + unmes >= 3:
            fracme = mes*1.0/(mes+unmes)
            tmpbed[0,(fracme>bins[0:10]) & (fracme<=bins[1:11])] += 1
    metricdata.loc[len(metricdata)] = [uid]+tmpbed[0,].tolist()
metricdata.to_csv(path_out + "_Mehist.txt",sep = "\t",index = False,header = True)       
        