# -*- coding: utf-8 -*-
"""
"""

import pysam,math
import getopt,sys
import numpy as np


options = "i:o:c:t:"
long_options = ["input=","output_txt=","cpgfile=","trim5="]
arguments, values = getopt.getopt(sys.argv[1:], options, long_options)
for currentArgument, currentValue in arguments:
    if currentArgument in ("-i", "--input"):
        path_in = currentValue
    elif currentArgument in ("-o", "--output_txt"):
        path_out = currentValue
    elif currentArgument in ("-c", "--cpgfile"):
        path_cpgfile = currentValue
    elif currentArgument in ("-t", "--trim5"):
        trim_bp = int(currentValue)        
        
data = pysam.AlignmentFile(path_in,"rb")
cgpos = pysam.TabixFile(path_cpgfile)

out = open(path_out,"w")
out_info = open(path_out+".info","w")


def merge_paired(pos1,pos2,seq1,seq2):
    return dict(map(lambda m,n: (m,n), pos1 + pos2, [x for x in seq1] + [x for x in seq2]))


def merge_single(pos1,seq1):
    return dict(map(lambda m,n: (m,n), pos1 , [x for x in seq1] ))

def trim5(pos1,seq1,trim_bp,whichend):
    if trim_bp < len(pos1)  and trim_bp < len(seq1) and trim_bp !=0:
        if whichend == "5'":
            pos1 = pos1[trim_bp:]
            seq1 = seq1[trim_bp:]
        elif whichend == "3'":
            pos1 = pos1[:-1*trim_bp]
            seq1 = seq1[:-1*trim_bp]
    return pos1,seq1

def fetch_cpg(record,cg_wanted):
    cgstate = ""
    poss = []
    cgranks = []
    for row in cg_wanted:
        pos = int(row.split("\t")[1])
        cgrank = int(row.split("\t")[2])
        if pos in record.keys():
            cgstate += record[pos]
        else:
            cgstate += "."
        poss.append(pos)
        cgranks.append(cgrank)
    if len(cgstate) == 0:
        pos_start = None
        cgrank_start = None
        pos_end = None
    else: 
        pos_start = min(poss)
        pos_end = max(poss)
        cgrank_start = min(cgranks)
    return cgrank_start,pos_start,pos_end,cgstate

def rm_del(seq1,cigar):
    used = []
    seq2 = ""
    for row in cigar:
        if row[0] == 0:
            used += [True]*row[1]
        elif row[0] == 1:
            used += [False]*row[1]
    if len(seq1) == len(used):
        seq2 = "".join(np.array([x for x in seq1])[np.array(used) == True].tolist())
    else:
        seq2 = seq1
    return seq2


idx1_total = 0
idx2_noCpG = 0
idx3_stored = 0
pairmerge_flag = 0
for row in data:
    idx1_total +=1 
    if math.ceil(idx1_total/100000)*100000 == idx1_total:
        print("Processing " + str(idx1_total/1000000) + "M reads!\n")
    if row.is_read1 == True and row.is_reverse == False:
        seq1 = row.query_alignment_sequence
        pos1 = [x+1 for x in row.get_reference_positions()]
        start1 = min(pos1)
        end1 = max(pos1)
        if "I" in row.cigarstring:
            seq1 = rm_del(seq1,row.cigartuples)
        pos1,seq1 = trim5(pos1,seq1,trim_bp, "5'")
        read1 = row.query_name
        chr1 = row.reference_name############################
        forward1 = True
        if len(seq1) == len(pos1):
            mate1 = True
        else:
            mate1 = False
    elif row.is_read1 == True and row.is_reverse == True:
        seq1 = row.query_alignment_sequence
        pos1 = row.get_reference_positions()
        start1 = min(pos1)+1
        end1 = max(pos1)+1
        if "I" in row.cigarstring:
            seq1 = rm_del(seq1,row.cigartuples)
        pos1,seq1 = trim5(pos1,seq1,trim_bp, "3'")
        read1 = row.query_name
        chr1 = row.reference_name########################
        forward1 = False
        if len(seq1) == len(pos1):
            mate1 = True
        else:
            mate1 = False
    elif row.is_read2 == True and row.is_reverse == True:
        seq2 = row.query_alignment_sequence
        pos2 = [x+1 for x in row.get_reference_positions()]
        start2 = min(pos2)
        end2 = max(pos2)
        if "I" in row.cigarstring:
            seq2 = rm_del(seq2,row.cigartuples)
        pos2, seq2 = trim5(pos2,seq2,trim_bp,"3'")
        read2 = row.query_name
        chr2 = row.reference_name#############################
        forward2 = False
        if len(seq2) == len(pos2):
            mate2 = True
        else:
            mate2 = False
    elif row.is_read2 == True and row.is_reverse == False:    
        seq2 = row.query_alignment_sequence
        pos2 = [x for x in row.get_reference_positions()]
        start2 = min(pos2)+1
        end2 = max(pos2)+1
        if "I" in row.cigarstring:
            seq2 = rm_del(seq2,row.cigartuples)
        pos2, seq2 = trim5(pos2,seq2,trim_bp,"5'")  
        read2 = row.query_name
        chr2 = row.reference_name##########################
        forward2 = True
        if len(seq2) == len(pos2):
            mate2 = True
        else:
            mate2 = False                
    if pairmerge_flag == 1:
        if read1 == read2 and chr1 == chr2 and forward1 != forward2 and mate1 == True and mate2 == True:
            record = merge_paired(pos1,pos2,seq1,seq2)
            allpos = list(record.keys())
            cg_wanted = cgpos.fetch(chr1,min(allpos)-1,max(allpos))
            cgrank_start,pos_start,pos_end,cgstate = fetch_cpg(record,cg_wanted)
            if len(cgstate) > 0:
                idx3_stored += 1
                cgrank_end = cgrank_start + len(cgstate) -1
                if forward1 == True:
                    strand = "+"
                else:
                    strand = "-"
                out.write(chr1 + "\t" +str(cgrank_start)+"\t"+ str(cgrank_end) + "\t" +\
                          cgstate + "\t" + str(pos_start)+ "\t"+ str(pos_end)+ "\t" + strand + "\n")
            else:
                idx2_noCpG += 1
            pairmerge_flag = 0
        else:
            pairmerge_flag = 1
    else:
        pairmerge_flag = 1
out_info.write("Total read records\t"+str(idx1_total)+"\n")
out_info.write("Read pairs without CpGs\t"+str(idx2_noCpG)+"\n")
out_info.write("Read pairs with CpGs\t"+str(idx3_stored)+"\n")
out_info.write("discordant read records\t"+str(idx1_total - idx3_stored*2 - idx2_noCpG*2)+"\n")
out_info.close()
out.close()
data.close()
cgpos.close()
