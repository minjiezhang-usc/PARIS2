"""
sam2chimera.py
extract chimera so that the supplement is anchored to a specific location
output to both sam and bed formats, but the bed format contains the interaction. 

example command:
python ~/Documents/scripts/duplex/sam2chimera.py \
AMT_Stress_trim_nodup_bc07_hg38genrefdfamaddAligned_ACTB_hs45S.sam \
AMT_Stress_trim_nodup_bc07_hg38genrefdfamaddAligned_ACTB_hs45Sh18new.sam \
hs45S:4220-4300 15

python ~/bin/sam2chimera.py \
AMT_Stress_trim_nodup_bc07_hg38genrefdfamaddAligned_sorted.sam \
AMT_Stress_trim_nodup_bc07_hg38genrefdfamaddAligned_other_tohs45S.sam \
hs45S 15

python ~/bin/sam2chimera.py \
AMT_Stress_trim_nodup_bc07_hg38genrefdfamaddAligned_other_tohs45S.sam \
AMT_Stress_trim_nodup_bc07_hg38genrefdfamaddAligned_other_tohs45Sh18.sam \
hs45S:4220-4300 15 &

Artifacts in the datasets:
SHISAL2B

                    5'UTR
('GAPDH', 900)      0 (strong internal peak, high crosslinking efficiency?)
('RACK1', 900)     >10
('RPS10-NUDT3', 903) 3
('ACTB', 940)
('RPS16', 969)     >10
('RPL27', 976)      1
('RPS4X', 1048)     2
('H1-2', 1084)      0, HIST1H1C
('SLC1A5', 1093)   >20
('EEF1A1', 1215)    2
('RPL23A', 1302)    0
('RPL15', 1513)     1
('RPLP1', 1533)     8 ******
('RPL37', 1540)     0
('TUBA1B', 1545)   >15
('H3C2', 1548)      0, chr6:26031589-26032099
('RPS19', 1637)     3 reads in 5'UTR
('H1-4', 1718)      8 (HIST1H1E)
('ACTG1', 2682)     2, a few in 3'UTR, chr17:81509971-81512799
('SHISAL2B', 3032)

Analysis of PARIS data after mRNA enrichment, mapped to hg38genrefdfamadd
mv 293mRNAAligned.sortedByCoord.out.bam 293mRNAAligned_sorted.bam
python ~/bin/sam2chimera.py 293mRNAAligned_sorted.sam \
293mRNAAligned_other_tohs45S.sam hs45S 15 &

python ~/bin/sam2chimera.py 293mRNAAligned_other_tohs45S.sam \
293mRNAAligned_other_tohs45Sh18.sam hs45S:4220-4300 15 &

The data are much lower than what we had from the HEK1 total RNAs. Next, check
the mouse mRNA data, which should be much better quality, based on mapping info.

Also analyzed the virus data as follows:
/home/rcf-40/zhipengl/storage/kongpanl/VirusProj/R707_HeLa_VR/
R707_HeLa_VRAligned.sortedByCoord.out.bam

python ~/bin/sam2chimera.py R707_HeLa_VRAligned_sorted.sam \
R707_HeLa_VRAligned_toVR1197.sam VR1197 15


"""


import sys, re

if len(sys.argv) < 5:
    print("python sam2chimera.py insam outsam chr:start-end length")
    print("extracts chimera with supplements (SA:Z) in a particular RNA")
    print("now only process simple chimera, one match (M) on each side")
    print("modified bed output: chrom1 start1 end1 chrom2 start2 end2")
    sys.exit()

insam = open(sys.argv[1], 'r')
outsam = open(sys.argv[2], 'w')
outbed = open(sys.argv[2][:-3]+"bed", 'w')

RNAinfo = sys.argv[3]
length = int(sys.argv[4])

RNA, start, end = '', "unset", "unset"
if ":" in RNAinfo and "-" in RNAinfo:
    RNA, position = tuple(RNAinfo.split(":"))
    print(RNA,position)
    start, end = tuple([int(i) for i in position.split("-")])
else: RNA = RNAinfo


for line in insam:
    if line[0] == "@": outsam.write(line); continue
    record = line.split()
    RNAME1,POS1,CIGAR1,lasttag = record[2],int(record[3]),record[5],record[-1]
    if RNAME1 == RNA or lasttag[:2] != "SA": continue #remove nonchimeric
    RNAME2,POS2,STRAND2,CIGAR2 = tuple(lasttag.split(':')[-1].split(',')[:4])
    POS2 = int(POS2)
    #tag: SA:Z:ref,29,-,6H5M,17,0;

    if RNAME2 != RNA: continue

    M1, M2 = re.findall('\d+M', CIGAR1), re.findall('\d+M', CIGAR2)
    M1len, M2len = 0, 0
    if len(M1)==len(M2)==1: M1len, M2len = int(M1[0][:-1]), int(M2[0][:-1])
    if M1len>=length and M2len>=length:
        if start == "unset" or POS2>=start-M2len and POS2<=end:
            bed=[RNAME1,str(POS1),str(POS1+M1len),
                 RNAME2,str(POS2),str(POS2+M2len)]
            outsam.write(line)
            outbed.write('\t'.join(bed)+'\n')

insam.close()
outsam.close()
outbed.close()








