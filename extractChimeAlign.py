"""
minjiezhang123@gmail.com    2020-01-16

Extract the gapped or chimeric alignments that one arm mapped to rRNA and other arm mapped to protein coding gene.
Gene bed file bed file for protein coding gene is needed.


"""


#1. input and output setup
################################################################################
################################################################################
#this section sets up the input and output
import sys, numpy, os, re, itertools, random
from datetime import datetime
from itertools import combinations
from multiprocessing import Process, Lock, Manager
from collections import Counter

if len(sys.argv) < 3:
    print("Usage:       python extractChimeAligno.py inputsam gene.bed outputsam")
    print("inputsam:    two alignmnets is needed for chimeric reads.")
    print("gene.bed:    Gene info in bed format")
    sys.exit()

inputfile = sys.argv[1]
inputsam = open(inputfile, 'r'); del inputfile
annobed = open(sys.argv[2],'r')
outputsam = open(sys.argv[3],'w')
inputcount = 0
################################################################################


#2. function definitions
################################################################################

def timenow(): return str(datetime.now())[:-7]

###############getOverlap with more than 50% of Query
##a, reference bed file:
##chr1 - 36082 35276 FAM138A intron1 lncRNA
##b, query bed file:
##chr1 - 36082 35276
def getOverlap(a, b): 
    overlap = int(min(int(a[3]),int(b[3]))) - int(max(int(a[2]),int(b[2])))
    readslen = (int(b[3]) - int(b[2]))/2
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
## hs45S   4253    4277    H18     +
for line in annobed:
    align = line.split()
    ref = [align[0],align[4],align[1],align[2]]
annobed.close()
#print(ref)


########################## Process sam file to annotate each segment ###########################
print(str(datetime.now())[:-7], "Reading sam file ...")

for line in inputsam:
    align = line.split()
    header=''; firstline=''
    if line[0]=="@": outputsam.write(line); continue
    inputcount += 1
    if not inputcount%1000000:
        print(str(datetime.now())[:-7], "Processed", inputcount, "reads ...")
    #    outputbed.write(outstring); outstring='' #write output to free up memory
    
    if align[-1].split(',')[0] == "SA:Z:mm45S": #SA:Z:chr3,69434139,+,8M22H,255,0;
        chemisegs = align[-1].split(',')
        CHR, POS_chemi, STRAND_chemi, CIGAR_chemi = chemisegs[0].split(':')[-1], int(chemisegs[1]), chemisegs[2], chemisegs[3]
        seglens_chemi = [int(i[:-1]) for i in re.findall('[0-9]+M', mergeCIGAR(CIGAR_chemi))] #seg lengths
        l, r = POS_chemi, POS_chemi+sum(seglens_chemi)-1
        seg = [CHR, STRAND_chemi, str(l), str(r)]
        if getOverlap(ref, seg) == "Ture" and align[2].startswith('chr'):
            outputsam.write(line)

outputsam.close()
inputsam.close()
    
   
