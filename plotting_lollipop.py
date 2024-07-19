# -*- coding: utf-8 -*-
"""
"""

import h5py
import getopt,sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


options = "i:o:b:n:"
long_options = ["input=","output=","bed=","name="]
#path_in = "I:/tmp_download/SRR9070182_cgi_MeMatrix.hdf5"
#titlename = "SRR9070182"
#bedfile = "I:/Projects/P7.Methy-PanCancer/data_analysis/A5.HCC_WGBS/hyper_hypo_cgis_A5HCC.bed"
arguments, values = getopt.getopt(sys.argv[1:], options, long_options)
for currentArgument, currentValue in arguments:
    if currentArgument in ("-i", "--input"):
        path_in = currentValue
    elif currentArgument in ("-o", "--output"):
        path_out = currentValue
    elif currentArgument in ("-b", "--bed"):
        bedfile = currentValue    
    elif currentArgument in ("-n", "--name"):
        titlename = currentValue

h = h5py.File(path_in,"r")
carebeds = pd.read_csv(bedfile,sep = "\t",header = None)
for idx in range(carebeds.shape[0]):
    uid = carebeds.iloc[idx,0] + "_" + str(carebeds.iloc[idx,1]) + "_" + str(carebeds.iloc[idx,2])
    memat = h[uid][...]
    if memat.shape[0] >= 3:
        cpgs_num = memat.shape[1]
        reads_num = memat.shape[0] - 2
        plot_mat = np.zeros((reads_num,cpgs_num),dtype = int)
        reads_mat = memat[2:,:]
        reads_info = np.zeros((reads_num,3),dtype = int)
        for readi in range(reads_num):
            pos1 = np.min(np.nonzero(reads_mat[readi,]))
            pos2 = np.max(np.nonzero(reads_mat[readi,]))
            tmp_line = np.zeros((1,cpgs_num),dtype = int)
            tmp_line[0,pos1:(pos2 + 1)] = 1
            for plotlinei in range(reads_num):
                if np.max(plot_mat[plotlinei,:] + tmp_line)>=2:
                    pass
                else:
                    reads_info[readi,:] = [pos1,pos2,plotlinei]
                    plot_mat[plotlinei,:] += tmp_line[0,:]
                    break
        plot_line_num = np.max(reads_info[:,2])+1
        #fig
        fig = plt.figure(figsize = [cpgs_num/3 + 1,plot_line_num/3 + 3])
        
        ax1 = fig.add_subplot(plot_line_num + 9,1,(1,6))
        dis = memat[1,1:] - memat[1,:-1]
        x = np.arange(1.5,len(dis)+1.5,1)
        ax1.plot(x,dis,'k',linewidth = 1.0,zorder = 0)
        ax1.set(xlim=(0, cpgs_num +1), xticks=np.arange(1, cpgs_num))
        ax1.set_title(titlename +"_"+ uid)
        
        ax2 = fig.add_subplot(plot_line_num + 9,1,(7,plot_line_num + 9))
        for readi in range(reads_num):
            pos1 = reads_info[readi,0]
            pos2 = reads_info[readi,1]
            placeline = reads_info[readi,2]
            x = np.array([pos1+1, pos2+1])
            y = np.array([placeline+1, placeline+1])
            yt = plot_line_num -y + 1
            
            x2 = np.arange(pos1+1, pos2 + 2)
            y2 = np.repeat(placeline+1,pos2 - pos1 + 1)
            y2t = plot_line_num - y2 + 1
            
            c2 = reads_mat[readi,pos1:(pos2+1)]
            cz2 = np.zeros((pos2-pos1 + 1,3))
            cz2[c2==-1,:] = np.repeat([0.68,0.92,1],sum(c2==-1),axis = 0).reshape((sum(c2==-1),3),order = 'F')
            cz2[c2==0,:] = np.repeat([0.8,0.8,0.8],sum(c2==0),axis = 0).reshape((sum(c2==0),3),order = 'F')
            cz2[c2==1,:] = np.repeat([0.6,0.2,0],sum(c2==1),axis = 0).reshape((sum(c2==1),3),order = 'F')
            
            ax2.scatter(x2,y2t,s = 150,c = cz2,edgecolors = 'k',alpha = 1)
            ax2.plot(x,yt,'k',linewidth=1.0,zorder= 0)
        ax2.set(xlim=(0, cpgs_num+1), xticks=np.arange(1, cpgs_num+1),\
                ylim=(0, plot_line_num +1), yticks=np.arange(1, plot_line_num + 1))
        
        ax2.set_xlabel("# CpG sites")
        ax2.set_ylabel("# Reads")
        fig.savefig(path_out + titlename + "_" + uid + "_lollipop.pdf",dpi=300)
    else:
        print(titlename + ":" + uid + ": has no enough reads (no read covers 3 target CpGs).")
        