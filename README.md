# MeHist
Epiallele analysis for WGBS (whole genome bisulfite sequencing) data


# Citation and Contact
Xianglin Zhang<sup>†</sup>, Wei Zhang<sup>†</sup>, Jinyi Zhang, Xiuhong Lyu, Haoran Pan, Tianwei Jia, Ting Wang*,Xiaowo Wang*, Haiyang Guo* Simulation-Guided Pan-Cancer Analysis Identifies a Novel Regulator of CpG Island Hypermethylation Heterogeneity, Briefings in Bioinformatics (Accepted)

zhangxl.2015@tsinghua.org.cn (Xianglin Zhang, Shandong University)

# Usage
## Preparation 
* bam files output from [Bismark](https://github.com/FelixKrueger/Bismark), should be sorted by name if paired-end reads in order to merge paired-end reads into a single record (fragment)
* reference genome (.fasta) used to build the index of Bismark
* bed file which contains the genomic regions that you care, such as CpG islands
## Installation
* Modules of python are required:`pysam`, `numpy`, `multiprocessing`,`matplotlib`,`pandas`,`h5py`.
* tools used: `bismark`, `bgzip`, `tabix`, `samtools`.

## Scrtips
### 1. preparing CpG positions for the genome
Usage: `python building_cpg_index.py -i hg38.fa -o ./`
* `i`,  The path to reference sequences (.fa);
* `o`,  The path that you want to deposit the positions of CpG sites, files of each chromosome and whole-genome;  
Followed processing:    
* `bgzip cpgpos_wholegenome.pos;`
* `tabix -b 2 -e 2 cpgpos_wholegenome.pos.gz`

### 2. converting bam files to fragment records
Usage: `python A1.bam2record.py -i *.bam -o tmp_post.txt -c cpgpos_wholegenome.pos.gz -t 0 ./`
* `i`,  The path to bam file from Bismark;
* `o`,  The path and file name that you want to deposit;
* `c`,  CpG index that built in step 1;
* `t`,  Trim reads (bp);  
Followed processing:     
* `sort -k1,1 -k2,2n tmp_post.txt > tmp_post_sort.txt`
* `bgzip tmp_post_sort.txt`
* `tabix -p bed tmp_post_sort.txt.gz`

### 3. fragment DNAme records to hdf5 containing read-CpG matrix of genomic regions in bed file
Usage: `python A2.record2hdf5.py -i tmp_post.txt.gz -o tmp_bed_ -c cpgpos_wholegenome.pos.gz -b region.bed -t 24 ./`
* `i`,  The path to bam file from Bismark;
* `o`,  The path and prefixed hdf5 file name that you want to deposit;
* `c`,  CpG index that built in step 1;
* `b`,  Bed file containing regions that you care, such CpG islands;
* `t`,  threads for multiprocessing;  

### 4. calculating the distribution of methylation levels of reads for each bed region
Usage: `python A3.Mehist.py -i tmp_bed_MeMatrix.hdf5 -o tmp -b region.bed ./`
* `i`,  The path to hdf5 file ;
* `o`,  The path that you want to deposit;
* `b`,  Bed file containing regions that you care, such CpG islands;  
  
* Output: tmp_Mehist.txt: columns : [BedID,0tenth,1tenth,2tenth,3tenth,4tenth,5tenth,6tenth,7tenth,8tenth,9tenth]  

### 5. drawing lollipop
Usage: `python plotting_lollipop.py -i tmp_bed_MeMatrix.hdf5 -o tmp -b region.bed -n sample1./`
* `i`,  The path to hdf5 file ;
* `o`,  The path that you want to deposit;
* `b`,  Bed file that you want to draw lollipop;
* `n`,  Prefix sample name of lollipop plot;  
  
* Output: plot in PDF file.  
Each line could cotain multiple fragments which are not overlapped. Top panel shows distance(bp) between two adjacent CpG sites.  
<img src="https://github.com/vhang072/MeHist/blob/main/pic/Lollipop_example.png" width="1750" height="400">

