"""
minjiezhang123@gmail.com    2020-02-07
Python version 3.6
Generate bed file from genecode gtf file.
Transcript selection
    for protein coding genes: 	longest appris principal transcript
    for non coding genes:	    longest transcript
    
"""


#1. input and output setup
################################################################################
################################################################################
#this section sets up the input and output
import sys, numpy, os, re, itertools, random
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

if len(sys.argv) < 2:
    print("Usage:           python mRNAmegaCoverage.py inputbam UTRCDS_anno.bed window output")
    print("inputbam:        bam file with target genes")
    print("UTRCDS_anno.bed: bed file contains UTR and CDS(exon) infos. Can be downloaded from https://drive.google.com/open?id=1wHSC-mf1jNNClXrVqMugqVmDVT4Crxzz.")
    print("window:          window bins for target region")
    print("output:          output file")
    sys.exit()


inputbam = sys.argv[1]
annobed = open(sys.argv[2],'r')
winbin = int(sys.argv[3])
output = open(sys.argv[4],'w')

fiveutr = {};   startdict = {}; enddict = {};   threeutr = {}
genes = []


#2. Subfunctions
##############################################################################
def getTargetRegion(list,startsite,winbin):
    ## dict:{genename:[chr,start,end,strand,len], [chr,start,end,strand,len]..}
    regions = [];   sum = 0
    if startsite == "right":
        for i in range(len(list)-1,-1,-1):
            chr = list[i][0];   strand = list[i][3]
            if int(list[i][4]) + sum >= int(winbin):
                #print(list[i][2],winbin,str(sum))
                start = str(int(list[i][2])-int(winbin)+sum);    end = str(list[i][2]);   
                regions.append([chr,start,end,strand])
                break
            else:
                start = str(int(list[i][1])-1);    end = str(list[i][2]);
                regions.append([chr,start,end,strand])
                sum += int(list[i][4])
        return(regions)
    
    if startsite == "left":
        for i in range(0,len(list),1):
            chr = list[i][0];   strand = list[i][3]
            if int(list[i][4]) + sum >= int(winbin):
                #print(list[i][2],winbin,str(sum))
                end = str(int(list[i][1])+int(winbin)-sum-1);    start = str(int(list[i][1])-1);   
                regions.append([chr,start,end,strand])
                break
            else:
                start = str(int(list[i][1])-1);    end = str(list[i][2]);
                regions.append([chr,start,end,strand])
                sum += int(list[i][4])
        return(regions)

#list=[["chr1","8","28","+","21"],["chr1","38","57","+","20"]]
#print(getTargetRegion(list,"left",60))
    


def callCoverage(inputbam,TargetRegion,startsite,winbin):
    ## TargetRegion: [['chr1', '48', '57', '+'], ['chr1', '19', '28', '+']]
    ## output bed file
    coverage = []
    tempbed = open("temp.bed",'w')
    for line in TargetRegion:
        outstring = line[0]+'\t'+line[1]+'\t'+line[2]+'\t'+line[3]+'\n'
        tempbed.write(outstring)
    tempbed.close()
    os.system("sort -k1,1 -k2,2n -o %s %s" % ("temp.bed", "temp.bed"))
    os.system("bedtools coverage -split -d -a %s -b %s > %s" % ("temp.bed", inputbam, "temp_cover.bg"))
    ## get the coverage list
    count = len(open("temp_cover.bg",'rU').readlines())
    if startsite == "right":
        if count < winbin:
            for i in range(winbin-count):   coverage.append("0") 
        coveragebg = open("temp_cover.bg",'r')
        for line in coveragebg:
            coverage.append(str(line.split()[5]))
        coveragebg.close()
    
    if startsite == "left":
        coveragebg = open("temp_cover.bg",'r')
        for line in coveragebg:
            coverage.append(str(line.split()[5]))
        coveragebg.close()
        if count < winbin:
            for i in range(winbin-count):   coverage.append("0") 
    return(coverage)
    
#TargetRegion=[['chr1', '81021', '81030', '+'], ['chr1', '81035', '81040', '+']]
#print(callCoverage(inputbam,TargetRegion,"left",20))      
   
   
## change the position of list: back -> head; head -> back   
def listSwitch(list):
    listnew = []
    for i in range(len(list)-1,-1,-1):
        listnew.append(list[i])
    return(listnew)


## merger the coveage of to list
def mergerCov(list1,list2):
    listnew = []
    for i in range(len(list1)):
        sum = int(list1[i]) + int(list2[i])
        listnew.append(str(sum))
    return(listnew)


## using matplotlib to plot line chart
def plotGraph(size,xlab,ylab,output,dist):
    figure = plt.figure(figsize=(8,2))
    axes = plt.Axes(figure, [.3,.3,.6,.6])
    figure.add_axes(axes)
    plt.bar(range(0, size), dist, color='k')
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.yaxis.set_ticks_position('left')
    axes.xaxis.set_ticks_position('bottom')
    plt.title("Reads coverage of "+output, fontsize=10)
    plt.xlim(0, size)
    sizey = max(dist)
    plt.ylim(0,sizey)
    axes.yaxis.set_major_locator(ticker.LinearLocator(numticks=2))
    
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.savefig(output+'.pdf')
    plt.show()


    
#2. Processing anno.bed file
##############################################################################
##########################  Reading the anno.bed database  ###################
print(str(datetime.now())[:-7], "Reading the anno.bed database ...")
"""
chr1    1324802 1325102 +       CPTP,five_utr,protein_coding
chr1    1326836 1326910 +       CPTP,five_utr,protein_coding
chr1    1326911 1327032 +       CPTP,exon1,protein_coding
chr1    1327241 1327760 +       CPTP,exon2,protein_coding
chr1    1327761 1328896 +       CPTP,three_utr,protein_coding
"""
for line in annobed:
    chr, start, end, strand = line.split()[0],line.split()[1],line.split()[2],line.split()[3]
    genename, region, biotype = line.split()[4].split(',')[0], line.split()[4].split(',')[1], line.split()[4].split(',')[2]
    length = int(end) - int(start) + 1
    if biotype == "protein_coding":
        if genename not in genes: genes.append(genename)
        if region == "five_utr":
            if genename not in fiveutr: fiveutr[genename] = []
            coord = [chr,start, end, strand, length]
            fiveutr[genename].append(coord)
        elif region == "three_utr":
            if genename not in threeutr: threeutr[genename] = []
            coord = [chr,start, end, strand, length]
            threeutr[genename].append(coord)
        else:
            if genename not in startdict: startdict[genename] = []
            if genename not in enddict: enddict[genename] = []
            coord = [chr,start, end, strand, length]
            startdict[genename].append(coord)
            enddict[genename].append(coord)
annobed.close()
#print(fiveutr)


##########################  Get the target regions  ##########################
print(str(datetime.now())[:-7], "get the target region ...")
FiveUtrCov = []; ThreeUtrCov = []; StartCov = []; EndCov = []
regions = [];   coverage = []
for i in range(winbin): FiveUtrCov.append("0"); ThreeUtrCov.append("0"); StartCov.append("0"); EndCov.append("0"); 

for gene in genes:
    ## 5'UTR (fiveutr:{genename:[chr,start,end,strand,len]..})
    if gene in fiveutr:
        strand = fiveutr[gene][0][3]
        if strand == "+":
            regions = getTargetRegion(fiveutr[gene],"right",winbin)
            coverage = callCoverage(inputbam,regions,"right",winbin)
    
        if strand == "-":
            regions = getTargetRegion(fiveutr[gene],"left",winbin)
            coverage = listSwitch(callCoverage(inputbam,regions,"left",winbin))
        FiveUtrCov = mergerCov(FiveUtrCov,coverage)    
        #print(coverage)
        #print(FiveUtrCov)
    
    
    ## Start (startdict:{genename:[chr,start,end,strand,len]..})
    if gene in startdict:
        strand = startdict[gene][0][3]
        if strand == "+":
            regions = getTargetRegion(startdict[gene],"left",winbin)
            coverage = callCoverage(inputbam,regions,"left",winbin)
    
        if strand == "-":
            regions = getTargetRegion(startdict[gene],"right",winbin)
            coverage = listSwitch(callCoverage(inputbam,regions,"right",winbin))
        StartCov = mergerCov(StartCov,coverage)    
        #print(coverage)
        #print(StartCov)
    
    
    ## End (enddict:{genename:[chr,start,end,strand,len]..})
    if gene in enddict:
        strand = enddict[gene][0][3]
        if strand == "+":
            regions = getTargetRegion(enddict[gene],"right",winbin)
            coverage = callCoverage(inputbam,regions,"right",winbin)
    
        if strand == "-":
            regions = getTargetRegion(enddict[gene],"left",winbin)
            coverage = listSwitch(callCoverage(inputbam,regions,"left",winbin))
        EndCov = mergerCov(EndCov,coverage)    
        #print(coverage)
        #print(EndCov)
    
    
    ## 3'UTR (enddict:{genename:[chr,start,end,strand,len]..})
    if gene in threeutr:
        strand = threeutr[gene][0][3]
        if strand == "+":
            regions = getTargetRegion(threeutr[gene],"left",winbin)
            coverage = callCoverage(inputbam,regions,"left",winbin)
    
        if strand == "-":
            regions = getTargetRegion(threeutr[gene],"right",winbin)
            coverage = listSwitch(callCoverage(inputbam,regions,"right",winbin))
        ThreeUtrCov = mergerCov(ThreeUtrCov,coverage)    
        #print(coverage)
        #print(ThreeUtrCov)

#print(FiveUtrCov)
#print(StartCov)
#print(EndCov)
#print(ThreeUtrCov)

     
##########################  Plot line chart   ##########################
print(str(datetime.now())[:-7], "plot line chart ...")

## 5'UTR + Start
size = int(len(FiveUtrCov)+len(StartCov))
xlab = "Start Codon"
ylab = "Coverage"
out = "5UTR_START"
dist = FiveUtrCov + StartCov
dist = list(map(int,dist))
#print(dist)
plotGraph(size,xlab,ylab,out,dist)
dist = []

## End + 3'UTR
size = int(len(EndCov)+len(ThreeUtrCov))
xlab = "End Codon"
ylab = "Coverage"
out = "END_3UTR"
dist = EndCov + ThreeUtrCov
dist = list(map(int,dist))
#print(dist)
plotGraph(size,xlab,ylab,out,dist)
dist = []

output.write('5UTR:\n' + '-'.join(FiveUtrCov) + '\n')
output.write('Start:\n'+ '-'.join(StartCov) + '\n')
output.write('End:\n'+ '-'.join(EndCov) + '\n')
output.write('Three:\n'+ '-'.join(ThreeUtrCov) + '\n')
output.close()

os.system("rm -f %s %s" % ("temp.bed", "temp_cover.bg"))
