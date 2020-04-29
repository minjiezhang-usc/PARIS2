# PARIS2 analysis strategy

1, The masked genome indicies

In order to accurately and easily analyze PARIS data, pseudogenes and multicopy genes from gencode, refGene and Dfam were masked from hg38/mm10 genome. And then single copy of them was added back as a separated “chromosome”. For example, multicopy of snRNAs were masked from the basic hg38/mm10 assembly genome, and 9 snRNAs (U1, U2, U4, U5, U6, U11, U12, U4atac and U6atac) were concatenated into one reference, separated by 100nt “N”s, was added back. The curated hg38/mm10 genome contained 25 reference sequences, or “chromosomes”, masked the multicopy genes and added back single copies. This reference is best suited for the PARIS analysis. 
The curated genome of hg38/mm10 can be downloaded from https://drive.google.com/open?id=1wHSC-mf1jNNClXrVqMugqVmDVT4Crxzz


2, Mapping

Reads were mapped to manually curated hg38 or mm10 genome using STAR program(Dobin, Davis et al. 2013). 

Global profiling of ribosome small subunit (SSU) analysis:
     
     STAR --runThreadN 8 --runMode alignReads --genomeDir OuputPath --readFilesIn SampleFastq  --outFileNamePrefix Outprefix --genomeLoad NoSharedMemory outReadsUnmapped Fastx  --outFilterMultimapNmax 10 --outFilterScoreMinOverLread 0 --outSAMattributes All --outSAMtype BAM Unsorted SortedByCoordinate --alignIntronMin 1 --scoreGap 0 --scoreGapNoncan 0 --scoreGapGCAG 0 --scoreGapATAC 0 --scoreGenomicLengthLog2scale -1 --chimOutType WithinBAM HardClip --chimSegmentMin 5 --chimJunctionOverhangMin 5 --chimScoreJunctionNonGTAG 0 --chimScoreDropMax 80 --chimNonchimScoreDropMin 20

Global profiling of spliceosomal snRNP binding sites:

    STAR --runThreadN 8 --runMode alignReads --genomeDir OuputPath --readFilesIn SampleFastq --outFileNamePrefix Outprefix --genomeLoad NoSharedMemory --outReadsUnmapped Fastx  --outFilterMultimapNmax 100 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --outSAMattributes All --outSAMtype BAM Unsorted SortedByCoordinate --alignIntronMin 1 --scoreGap 0 --scoreGapNoncan 0 --scoreGapGCAG 0 --scoreGapATAC 0 --scoreGenomicLengthLog2scale -1 --chimOutType WithinBAM HardClip Junctions --chimSegmentMin 5 --chimJunctionOverhangMin 5 --chimScoreJunctionNonGTAG 0 --chimScoreDropMax 80 --chimNonchimScoreDropMin 80 --chimScoreSeparation 0 --chimSegmentReadGapMax 30 --limitBAMsortRAM 20000000000



3, Global profiling of ribosome small subunit (SSU) analysis

3.1 Extracting mRNA-rRNA chimeric alignments using awk and sam2mRNArRNAchimera.py:

    samtools view -h SampleAligned.sortedByCoord.out.bam | awk  '$1~/^@/ || $0~/hs45S/' > Sample_rRNA.sam
    python sam2mRNArRNAchimera.py Sample_rRNA.sam hg38_UTRCDS_anno.bed hs45S Sample_mRNArRNA.sam

3.2 Analyzing the binding sites of mRNAs on the hs45S:
    
    Sample_mRNArRNA.sam can be converted to bam file and loaded to IGV to check the mRNA-rRNA binding sites on rRNA.

3.3 Analyzing he binding sites of h18 and h26 on the meta mRNA:
    
    python mRNAmegaCoverage.py Sample_mRNArRNA.bam hg38_UTRCDS_anno.bed 200 mmH26.dist  


4, Global profiling of spliceosomal snRNP binding sites

4.1 Extracting snRNA-mRNA interaction chimeric alignments:
     
    samtools view -h SampleAligned_sorted.bam hssnRNA | awk  '$1~/^@/ || $21!~/hssnRNA/ && $21~"SA:Z" && $21~",0;"' > Sample_tohssnRNA.sam

4.2 snRNA-target interaction alignments with at least 15 nt matches for the snRNA targets were filtered using filterchimera.py

    python filterchimera.py Sample_tohssnRNA.sam Sample_tohssnRNA_15M.sam 15

4.3 200 nt windows around splice sites was extracted for gencode gtf file using gtf2splice.py. 

    python gtf2splice.py gencode.v33.annotation.gtf gencode.v33.annotation_splice200nt 200

4.4 Chimera connecting specific snRNA regions were further extracted using sam2chimera.py script. 

    python sam2chimera.py Sample_tohssnRNA_15M.sam Sample_tohssnRNA_15M_U1_1to12_15M.sam hssnRNA:1-12 15

4.4 Coverage along the 200nt windows was calculated using bedtools coverage. 

    bedtools coverage -d -a gencode.v33.annotation_splice200nt_splice3.bed -b Sample_tohssnRNA_15M_U1_1to12_15M.sam > Sample_tohssnRNA_15M_U1_1to12_15M_splice200nt_splice3.bed
    bedtools coverage -d -a gencode.v33.annotation_splice200nt_splice5.bed -b Sample_tohssnRNA_15M_U1_1to12_15M.sam > Sample_tohssnRNA_15M_U1_1to12_15M_splice200nt_splice5.bed

4.5 Meta-analysis for all windows around start 5’ and 3’ splice sites were performed with windowmeta.py. 

    python windowmeta.py 400 Sample_tohssnRNA_15M_U1_1to12_15M_splice200nt_splice3.bed Sample_tohssnRNA_15M_U1_1to12_15M_splice200nt_splice5_sum.bedgraph
    python windowmeta.py 400 Sample_tohssnRNA_15M_U1_1to12_15M_splice200nt_splice5.bed Sample_tohssnRNA_15M_U1_1to12_15M_splice200nt_splice5_sum.bedgraph

4.6 The Output.bedgraph can be loaded to IGV for visualization.
