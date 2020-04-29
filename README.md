# PARIS2
1, Peparing the masked genome indicies
In order to accurately and easily analyze PARIS data, pseudogenes and multicopy genes from gencode, refGene and Dfam were masked from hg38/mm10 genome. And then single copy of them was added back as a separated “chromosome”. For example, multicopy of snRNAs were masked from the basic hg38/mm10 assembly genome, and 9 snRNAs (U1, U2, U4, U5, U6, U11, U12, U4atac and U6atac) were concatenated into one reference, separated by 100nt “N”s, was added back. The curated hg38/mm10 genome contained 25 reference sequences, or “chromosomes”, masked the multicopy genes and added back single copies. This reference is best suited for the PARIS analysis. 
The EV_D68 viral genome (GenBank, KM851225.1) was downloaded from NCBI and manually corrected based on our viral sequencing data. After mutation identifying using GATK software, three variant sites on EV_D68 genome were corrected (2023:G->A; 2647:G->A; 3242:A->G). The curated EV_D68 genome was added to hg38 refence as an independent chromosome. 

2, 


2, Global profiling of spliceosomal snRNP binding sites

2.1, Extracting snRNA-mRNA interaction chimeric alignments:
samtools view -h AMT_Stress_trim_nodup_bc07-bc12_hg38mask14addAligned_sorted.bam hssnRNA | awk  '$1~/^@/ || $21!~/hssnRNA/ && $21~"SA:Z" && $21~",0;"' > AMT_Stress_trim_nodup_bc07-bc12_hg38mask14addAligned_hssnRNA_toother.sam

snRNA-target interaction alignments with at least 15 nt matches for the snRNA targets were filtered using filterchimera.py. 200 nt windows around splice sites was extracted for gencode gtf file using gtf2splice.py. Chimera connecting specific snRNA regions were further extracted using sam2chimera.py script. Coverage along the 200nt windows was calculated using bedtools coverage. Meta-analysis for all windows around start 5’ and 3’ splice sites were performed with windowmeta.py. The Output.bedgraph can be loaded to IGV for visualization.
