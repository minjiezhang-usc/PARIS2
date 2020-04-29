"""
windowmeta.py
produce the metaprofile of windows in a genome, considering the strand.

python ~/bin/windowmeta.py 400 AMT_Stress_trim_nodup_bc07-bc12_\
hg38mask14addAligned_tohssnRNA_U1_1to12_15M_splice200nt_splice3_non0.bed \
AMT_Stress_trim_nodup_bc07-bc12_\
hg38mask14addAligned_tohssnRNA_U1_1to12_15M_splice200nt_splice3_sum.bed



"""


import sys
import numpy as np

if len(sys.argv) < 3:
    print("Usage: python windowmeta.py 400 inputbed bedgraph")
    sys.exit()

windowlen = int(sys.argv[1])
inputbed = open(sys.argv[2], 'r')
bedgraph = open(sys.argv[3], 'w')
output = []
tempout = open("tempout", 'w')

#collect windows that have at least 1 read. 
windowdict = {}
for line in inputbed:
    chrom, start, end, name, anchor, strand, relpos, cov = tuple(line.split())
    window = '_'.join([chrom, start, end])
    if window not in windowdict: windowdict[window] = [strand, anchor]
    windowdict[window].append(int(cov))
inputbed.close()


#remove 0-windows.
sumcov = 0
windowdict1 = {}
for window in windowdict:
    #print(window, sum(windowdict[window][2:]))
    if sum(windowdict[window][2:])>0 and len(windowdict[window])-2==windowlen:
        windowdict1[window]=windowdict[window]
        sumcov += sum(windowdict[window][2:])
print("All and non-zero windows:", len(windowdict), len(windowdict1))
print("Total coverage:", sumcov)



#output non-zero windows
inputbed = open(sys.argv[2], 'r') 
for line in inputbed:
    chrom, start, end = tuple(line.split()[:3])
    window = '_'.join([chrom, start, end])
    if window in windowdict1: tempout.write(line)
tempout.close()


#summarize each position for the windowlist
windowlist = []
covdict = {}
for window in windowdict1:
    strand, anchor = tuple(windowdict1[window][:2])
    values = windowdict1[window][2:]
    if strand == '-': values = values[::-1]
    windowlist.append(values)
    covdict[window] = sum(values)

w = [str(sum(i)) for i in np.transpose(windowlist)]
windowsum = ["chr1\t"+str(i)+"\t"+str(i+1)+"\t"+w[i] for i in range(len(w))]
covlist = sorted([(v,k) for (k,v) in covdict.items()])
for i in covlist[-20:]: print(i[1], i[0]) #export top coverage RNAs. 

inputbed.close()
bedgraph.write('\n'.join(windowsum))
bedgraph.close()













    
