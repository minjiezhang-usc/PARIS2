"""
minjiezhang123@gmail.com    2020-01-16
Extract the chimeric alignments that one arm mapped to rRNA and other arm mapped to protein coding gene.
Gene annotation bed file is needed.
"""


#1. input and output setup
################################################################################
################################################################################
#this section sets up the input and output
import sys, numpy, os, re, itertools, random
from datetime import datetime

if len(sys.argv) < 3:
    print("Usage:       python sam2mRNArRNAchimera.py inputsam gene.bed hs45S/mm45S outputsam")
    print("inputsam:    two alignmnets is needed for chimeric reads.")
    print("anno.bed:    Gene annotation info in bed format")
    print("             can be downloaded from https://drive.google.com/open?id=1wHSC-mf1jNNClXrVqMugqVmDVT4Crxzz")
    print("hs45S/mm45S: target region interacted with mRNA")
    print("             hs45S: to extract human rRNA-mRNA chimera")
    print("             mm45S: to extract mouse rRNA-mRNA chimera")
    print("outputsam:   Output sam file")
    sys.exit()

inputfile = sys.argv[1]
inputsam = open(inputfile, 'r'); del inputfile
annobed = open(sys.argv[2],'r')
chimerictarget = sys.argv[3]
outputsam = open(sys.argv[4],'w')
inputcount = 0
################################################################################


#2. function definitions
################################################################################

def timenow(): return str(datetime.now())[:-7]

###############getOverlap with more than 50% of Query
##a, list of reference bed file:
##chr1 36082 35276 FAM138A 1000 +
##b, query bed file:
##chr1 36082 35276 -
def getOverlap(a, b):
    if a[0] == b[0]:
        overlap = int(min(int(a[2]),int(b[2]))) - int(max(int(a[1]),int(b[1])))
        readslen = (int(b[2]) - int(b[1]))/2
        if overlap >1 and overlap >= readslen: seganno = "Ture"
        else: seganno = "False"
        return seganno
 

def mergeCIGAR(CIGAR): 
    #merge all operations that consume the reference, i.e. MI=X
    #example: 1S2M3N4M5I6M7S -> 1S2M3N10M7S 
    ops = re.findall('\d+[MNISH=X]', CIGAR) #all that consume query
    newops = [ops[0]]
    for op in ops[1:]: #concatenate all internal ops that consume query [MIS=X
        if op[-1] not in "I=X":
            if newops[-1][-1]=="M" and op[-1]=="M":
                newops[-1] = str(int(newops[-1][:-1])+int(op[:-1]))+"M"
            else: newops.append(op)
    newCIGAR = ("".join(str(i) for i in newops))
    return newCIGAR

  
### remove duplicated element form list      
def deleteDuplicatedElementFromList(listA):
    return sorted(set(listA), key = listA.index)


#3. Processing sam file
################################################################################################
########################## Read ref gene bed file ###########################
print(str(datetime.now())[:-7], "Reading ref gene bed file ...")
refs=[]
## hs45S   4253    4277    H18     1000 +
for line in annobed:
    align = line.split()
    genename = align[4].split(',')[0]
    if align[4].split(',')[2] == "protein_coding":
        refs.append((align[0],align[1],align[2],align[3]))
annobed.close()
#print(refs)


########################## Process sam file to annotate each segment ###########################
print(str(datetime.now())[:-7], "Reading sam file ...")
for line in inputsam:
    align = line.split()
    header=''; firstline=''; segs=[]
    if line[0]=="@": outputsam.write(line); continue
    inputcount += 1
    if not inputcount%1000000:
        print(str(datetime.now())[:-7], "Processed", inputcount, "reads ...")
    
    if align[-1].split(',')[0] == "SA:Z:"+chimerictarget:
        CHR, POS, CIGAR = align[2], int(align[3]), align[5]
        STRAND = '-' if '{0:012b}'.format(int(align[1]))[-5] == '1' else '+'
        segs = [] #store all gaps from this CIGAR string, each as a 3-tuple.
        seglens = [int(i[:-1]) for i in re.findall('[0-9]+M', mergeCIGAR(CIGAR))] #seg lengths
        Ns =[i.rstrip('0123456789') for i in mergeCIGAR(CIGAR).split('M')]
        rx = [] #reference consumed: MD=X
        for N in Ns:
            rx.append(sum([int(i[:-1]) for i in re.findall('[0-9]+[ND=X]', N)]))
        for i in range(len(seglens)): #combine ref and segment lengths to make the junctions
            l, r = POS+sum(rx[:i+1])+sum(seglens[:i]), POS+sum(rx[:i+1])+sum(seglens[:i+1])-1
            segs.append((CHR, str(l), str(r), STRAND))
                
    if align[2] == chimerictarget:
        chemisegs = align[-1].split(',')
        CHR, POS_chemi, STRAND_chemi, CIGAR_chemi = chemisegs[0].split(':')[-1], int(chemisegs[1]), chemisegs[2], chemisegs[3]
        seglens_chemi = [int(i[:-1]) for i in re.findall('[0-9]+M', mergeCIGAR(CIGAR_chemi))] #seg lengths
        Ns_chemi =[i.rstrip('0123456789') for i in mergeCIGAR(CIGAR_chemi).split('M')]
        rx_chemi = [] #reference consumed: MD=X
        for N in Ns_chemi:
            rx_chemi.append(sum([int(i[:-1]) for i in re.findall('[0-9]+[ND=X]', N)]))
        for i in range(len(seglens_chemi)): #combine ref and segment lengths to make the junctions
            l, r = POS_chemi+sum(rx_chemi[:i+1])+sum(seglens_chemi[:i]), POS_chemi+sum(rx_chemi[:i+1])+sum(seglens_chemi[:i+1])-1
            segs.append((CHR, str(l), str(r), STRAND_chemi))
    
    for seg in segs:
        for ref in refs:
            if getOverlap(ref, seg) == "Ture": outputsam.write(line); break

outputsam.close()
inputsam.close()
    
   
