#!/usr/bin/env python

'''This script uses a tab delimted input .txt file that contains three columns. Column1= sampleName, Column2 = pathToInputFile/Input_file.bam, and Column3 = pathToIPFile/IP_file.bam.
This will run macs2 SPMR and then generate the bigwig track of the fold-enrichment of the IP/Input
For H3K27me3, use --broad option for macs2 callpeak 
To run script, do: python macs2SPMR.py #ofProcessors'''

import sys
import pandas as pd
from multiprocessing.pool import Pool
import os


df = pd.read_csv('Sample_Names.txt', sep='\t')


def run_macs2(x):
    cmd = 'macs2 callpeak -t {ip} -c {input_f} -f BAM -B -n {name} -g hs --nomodel -q 0.0001 --bw 250 --SPMR'
    #cmd = 'macs2 callpeak -t {ip} -c {input_f} -f BAM -B -n {name} -g hs --nomodel --bw 250 --SPMR' --broad #For H3K27me3
    cmd2 = 'macs2 bdgcmp -t {name}_treat_pileup.bdg -c {name}_control_lambda.bdg -o {name}_FE.bdg -m FE'.format(name=x[0])

    out_cmd2 = cmd2.split('-o ')[1].split(' -m')[0]
    cmd3 = '/bdg2bw {cmd2_out} /GenomeReference/hg19.chrom.sizes '.format(cmd2_out=out_cmd2)
    cmd_to_run = 'echo {};'.format(x[0])+ cmd.format(ip=x[2], input_f=x[1], name=x[0]) +' ; '+cmd2 +' ; '+cmd3
    
    os.system(cmd_to_run)
    return cmd_to_run




p =Pool(int(sys.argv[1]))

p.map(run_macs2, df.values.tolist())
    


