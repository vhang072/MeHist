# -*- coding: utf-8 -*-
"""
"""

import pysam
import getopt,sys
import numpy as np
import pandas as pd
import multiprocessing
import h5py

options = "i:o:c:b:t:"
long_options = ["input=","output=","cpgfile=","bed=","cores="]
arguments, values = getopt.getopt(sys.argv[1:], options, long_options)
for currentArgument, currentValue in arguments:
    if currentArgument in ("-i", "--input"):
        path_in = currentValue
    elif currentArgument in ("-o", "--output"):
        path_out = currentValue
    elif currentArgument in ("-c", "--cpgfile"):
        path_cpgfile = currentValue
    elif currentArgument in ("-b", "--bed"):
        path_bedfile = currentValue 
    elif currentArgument in ("-t", "--cores"):
        cores = int(currentValue) 

global cpgs
global pat
global bed 
cpgs = pysam.TabixFile(path_cpgfile)
pat = pysam.TabixFile(path_in)
bed = pd.read_csv(path_bedfile,sep="\t",header = None)

def readintegrate(read, bedcpg_info):
    read_strand = read.split("\t")[6]
    read_start = int(read.split("\t")[1])
    read_info = read.split("\t")[3]
    onereadme = np.zeros((1,bedcpg_info.shape[0]),dtype = np.int8)
    if read_strand == "+":
        for idx in range(0,len(read_info)):
            if idx+read_start in bedcpg_info["cpg_rank"].tolist():
                if read_info[idx] == "C":
                    onereadme[0,bedcpg_info.iloc[:,1] == idx+read_start] = 1
                elif read_info[idx] == "T":
                    onereadme[0,bedcpg_info.iloc[:,1] == idx+read_start] = -1
    elif read_strand == "-":
        for idx in range(0,len(read_info)):
            if idx+read_start in bedcpg_info["cpg_rank"].tolist():
                if read_info[idx] == "G":
                    onereadme[0,bedcpg_info.iloc[:,1] == idx+read_start] = 1
                elif read_info[idx] == "A":
                    onereadme[0,bedcpg_info.iloc[:,1] == idx+read_start] = -1
    return onereadme

def eachbedline(idx):
    chrom_query = bed.iloc[idx,0]
    start_query = bed.iloc[idx,1]
    end_query = bed.iloc[idx,2]
    uid = chrom_query + "_" + str(start_query) + "_" + str(end_query)
    cpg_ranks = []
    cpg_poss = []
    for site in cpgs.fetch(chrom_query,start_query,end_query,multiple_iterators=True):
        cpg_ranks.append(int(site.split("\t")[2]))
        cpg_poss.append(int(site.split("\t")[1]))
    if len(cpg_ranks) > 0:
        chroms = [chrom_query]*len(cpg_ranks)
        bedcpg_info = pd.DataFrame(data = {"chrom":chroms,"cpg_rank":cpg_ranks,"cpg_pos":cpg_poss})
        readsme = np.array([cpg_ranks,cpg_poss])
        for read in pat.fetch(chrom_query,min(cpg_ranks)-1,max(cpg_ranks),multiple_iterators=True):
            onereadme = readintegrate(read,bedcpg_info)
            readsme = np.vstack((readsme,onereadme))
    else:
        readsme = np.zeros(0,dtype = np.int8)
    return readsme,uid
   
        
if __name__ == '__main__':
    pool  = multiprocessing.Pool(processes=cores)
    
    f = pool.map(eachbedline,range(len(bed)),chunksize=60)
    pool.close()
    
    h = h5py.File(path_out + "_MeMatrix.hdf5","w")
    for val in f:
        uid = val[1]
        memat = val[0]
        datset = h.create_dataset(uid,data = memat)
    h.close()
    


