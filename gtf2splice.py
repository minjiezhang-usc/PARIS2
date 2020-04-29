"""
gtf2splice.py.

Only take the genomic regions, without putting together different exons,
not even the micro-exons.

For each splice site, determine the 5' or 3'. 

cd Downloads
head -n 10000 gencode.v33.annotation.gtf > gencode.v33.annotation10000.gtf
python ~/Documents/scripts/duplex/gtf2splice.py gencode.v33.annotation10000.gtf\
 gencode.v33.annotation_10000_splice200nt 200

python ~/Documents/scripts/duplex/gtf2splice.py gencode.v33.annotation.gtf \
gencode.v33.annotation_splice200nt 200
"""

import sys

if len(sys.argv)<4:
    print("Usage: python gtf2splice.py gtf bed window")
    sys.exit()

gtf = open(sys.argv[1], 'r')
bedprefix = sys.argv[2]
bedsplice5 = open(bedprefix+"_splice5.bed", 'w')
bedsplice3 = open(bedprefix+"_splice3.bed", 'w')
window = int(sys.argv[3])


PAR = ["AKAP17A","ASMT","ASMTL","CD99","CRLF2","CSF2RA","DHRSX","GTPBP6",
       "IL3RA","P2RY8","PLCXD1","PPP2R3B","SHOX","SLC25A6","XG","ZBED1",
       "IL9R","SPRY3","VAMP7","CXYorf1"] #pseudoautosomal (PAR) genes to ignore

#1. take gtf and extract one primary transcript for each protein-coding gene.
#for each exon, except the first and last, take the information for each splice
exonsdict = {}
#gene:{transcript:[strand,exon1,exon2,...]}
for line in gtf:
    record = line.split()
    if line[0] != "#" and "appris_principal_1" in line and \
       record[13][1:-2] == record[17][1:-2] == "protein_coding":
        gene = record[15][1:-2]
        if gene in PAR: continue #ignore pseudoautosomal genes
        transcript = record[19][1:-2]
        coords = [record[0], int(record[3]), int(record[4])]
        if gene not in exonsdict: exonsdict[gene] = {}
        if transcript not in exonsdict[gene]: #initialize exon,start,stop,strand
            exonsdict[gene][transcript] = [record[6]] 
        if record[2] == "exon":
            exonsdict[gene][transcript].append(coords)   
print("Number of protein_coding genes:", len(exonsdict))
#print(list(exonsdict))
#for i in list(exonsdict.items()): print(i)
            
#['OR4F5' +, 'OR4F29' -, 'OR4F16' -, 'NOC2L' -, 'KLHL17' +, 'HES4' -]
#example for 'OR4F5', which has one transcript 'OR4F5-201': 
#[[['chr1',69055,70108]], ['chr1',69091,69093,1], ['chr1',70006,70008,1]]
#refGene model for OR4F5 is incorrect. The Gencode model is chr1:69055-70108

metasplice3 = []
metasplice5 = []
for gene in exonsdict:
    transcript1 = list(exonsdict[gene].keys())[0]
    strand = exonsdict[gene][transcript1][0]
    exoncount = len(exonsdict[gene][transcript1])-1
    if exoncount == 1: continue #ignore single exon genes
    exons = sorted(exonsdict[gene][transcript1][1:], key=lambda x: x[1])
    genesplice3, genesplice5 = [], []
    for exon in exons:
        splice3, splice5 = [], []
        if strand == '+':
            splice3 = [exon[0],exon[1]-200,exon[1]+200,gene,"splice3",'+']
            splice5 = [exon[0],exon[2]-200,exon[2]+200,gene,"splice5",'+']
        elif strand == '-':
            splice3 = [exon[0],exon[2]-200,exon[2]+200,gene,"splice3",'-']
            splice5 = [exon[0],exon[1]-200,exon[1]+200,gene,"splice5",'-']
        genesplice3.append(splice3)
        genesplice5.append(splice5)

    #take the internal splice sites (transcript starts and ends are removed)
    if strand=='+': metasplice3+=genesplice3[1:];metasplice5+=genesplice5[:-1]
    elif strand=='-': metasplice3+=genesplice3[:-1];metasplice5+=genesplice5[1:]

    #print(strand, len(exons), exons)
    #print(len(genesplice3), genesplice3)
    #print(len(genesplice5), genesplice5)
    #if gene=="PLCXD1": print(exonsdict[gene])

#print('\n', metasplice3, '\n')
#print('\n', metasplice5, '\n')

str1 = '\n'.join(['\t'.join([str(j) for j in i]) for i in metasplice3])
str2 = '\n'.join(['\t'.join([str(j) for j in i]) for i in metasplice5])

gtf.close()
bedsplice3.write(str1)
bedsplice3.close()
bedsplice5.write(str2)
bedsplice5.close()




















        
