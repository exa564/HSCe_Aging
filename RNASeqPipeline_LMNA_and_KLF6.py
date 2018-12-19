#!/usr/bin/env python

"""This pipeline is designed to align, run fastQC, and generate counts for human, first-stranded, paired-end, 75 bp RNA-seq libraries.
    This is the pipeline used to prcoess all CD34+ shLMNA and CD34+ sgKLF6 RNAseq data.
    One input configuration file 'RNASampleID.txt' is required. It must be a 3 column tab seperated file with a header, containing the following:
        Column 1 contains Sample_Name
        Column 2 contains /path/R1_fastq.gz
        Column 3 contain /path/R2_fastq/gz
    The following versions were used for the analysis:
        Bowtie v2.2.6
        Cutadapt v1.12
        Samtools v1.3.1
        STAR v2.5.2b
        QoRTs v1.0.7
    Data was aligned to gencode v19
    Prior to alignment, generate indices for STAR:
        STAR --runThreadN 30 --runMode genomeGenerate --genomeDir /PathToIndices/ --genomeFastaFiles hg19ERCC.fa  --sjdbGTFfile  gencode.v19.annotation.complete.gtf  --sjdbOverhang 51
"""


import sys
import pandas as pd
from multiprocessing.pool import Pool
import os



df = pd.read_csv('RNASampleID.txt', sep='\t')


for r1, r2 in df[['R1', 'R2']].values.tolist():
    if os.path.isfile(r1) and os.path.isfile(r2):
        pass
    else:
        print 'File not found', r1, r2
        exit()

create_counts_main = 'mkdir counts'

os.system(create_counts_main)

def run_rnaseq(x):
    cmd1 =  'cutadapt -a AGATCGGAAGAGC -A AGATCGGAAGAGC -m 72 --length 73 -o {name}_R1_trim.fastq.gz  -p {name}_R2_trim.fastq.gz  {R1} {R2}'

    cmd1 = cmd1.format(name=x[0], R1=x[1], R2=x[2])
    print '=' * 30
    print cmd1
    os.system(cmd1)
    cmd2 = 'fastqc  --noextract -f fastq -t 5 {name_R1_trim} ;fastqc  --noextract -f fastq -t 5 {name_R2_trim}'

    trim_name_r1 = cmd1.split('-o ')[1].split(' ')[0]
    trim_name_r2 = cmd1.split('-p ')[1].split(' ')[0]

    cmd2 = cmd2.format(name_R1_trim=trim_name_r1, name_R2_trim=trim_name_r2)
    print '=' * 30
    print cmd2
    os.system(cmd2)
    cmd3 = 'STAR --runThreadN 25  --genomeDir /PathToIndices/ --readFilesIn {name_R1_trim} {name_R2_trim} --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFileNamePrefix {name}_ --alignEndsType EndToEnd --readFilesCommand gunzip -c  --outSAMtype BAM SortedByCoordinate --outWigType bedGraph --outWigStrand Stranded --outWigNorm RPM '


    cmd3 = cmd3.format(name=x[0] ,name_R1_trim=trim_name_r1, name_R2_trim=trim_name_r2)
    print '=' * 30
    print cmd3
    os.system(cmd3)
    
    cmd4 = 'samtools sort -n {name}_Aligned.sortedByCoord.out.bam -@ 20 -o {name}.sorted.bam'

    cmd4 = cmd4.format(name=x[0])
    print '=' * 30
    print cmd4
    os.system(cmd4)


    # Prior to calling counts, create gtf file 'gencode.v19_wo_rRNA.gtf' that does not contain rRNA, mitoRNA, or tRNA
    cmd5 = 'java -jar QoRTs.jar QC  --nameSorted --maxReadLength 75 --generatePlots --stranded --stranded_fr_firststrand {cmd4_out}  gencode.v19_wo_rRNA.gtf {counts_dir}/'

    sorted_bam = cmd4.split('-o ')[1]

    dir_count = 'counts/{name}'.format(name= x[0])
    cmd_create_dir = 'cd counts; mkdir {dir_var}; cd ../'.format(dir_var =  x[0])
    print '=' * 30
    print cmd_create_dir
    os.system(cmd_create_dir)

    cmd5 = cmd5.format(cmd4_out=sorted_bam,
                counts_dir = 'counts/{}'.format(x[0]) )
    print '=' * 30
    print cmd5
    os.system(cmd5)
    
     
p =Pool(int(sys.argv[1]))

p.map(run_rnaseq, df.values.tolist())  