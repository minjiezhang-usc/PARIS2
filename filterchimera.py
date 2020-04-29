"""
filterchimera.py

filter chimera so that each segment is at least a certain length. 
"""

import sys, re

if len(sys.argv) < 4:
    print("Usage: python filterchimera.py insam outsam length")
    sys.exit()

insam = open(sys.argv[1], 'r')
outsam = open(sys.argv[2], 'w')
length = int(sys.argv[3])

for line in insam:
    if line[0] == "@": outsam.write(line); continue
    record = line.split()
    CIGAR1, CIGAR2 = record[5], record[20].split(',')[3]
    M1 = [int(i[:-1])for i in re.findall('\d+M', CIGAR1)]
    M2 = [int(i[:-1])for i in re.findall('\d+M', CIGAR2)]
    print(M1, M2)
    if max(M1)>=length and max(M2)>=length: outsam.write(line)

insam.close()
outsam.close()
        
