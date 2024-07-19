# MeHist
epiallele analysis for WGBS (whole genome bisulfite sequencing) data

# Installation
* Modules of python are required:`pysam`, `numpy`, `multiprocessing`,`matplotlib`,`pandas`,`h5py`
* tools used: bismark, bgzip, tabix, samtools

# Citation
The paper is not published. It will be updated as soon as it is published.

# Contact
zhangxl.2015@tsinghua.org.cn (Xianglin Zhang, Shandong University)

# Usage
## preparation 
* bam files output from [Bismark] (https://github.com/FelixKrueger/Bismark), should be sorted by name if paired-end reads in order to merge paired-end reads into a single record (fragment)
* reference genome (.fasta) used to build the index of Bismark
* bed file which contains the genomic regions that you care, such as CpG islands

## 

