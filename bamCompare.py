#!/usr/bin/env python

''' Using deeptools2 bamCompare, this script will generate a log2(IP/Input) read normalized bigwig file.
    One input configuration file 'SampleID.txt' is required. It must be a 3 column tab seperated file containing the following:
        Column 1 contains Sample_Name
        Column 2 contains /path/Input.bam
        Column 3 contain /path/IP.bam
    To run script, do: python bamCompare.py 'number of threads'
'''

import sys
import pandas as pd
from multiprocessing.pool import Pool
import os


df = pd.read_csv('SampleID.txt', sep='\t')


def run_bamCompare(x):
    cmd = 'bamCompare --bamfile1 {ip} --bamfile2 {input_f} --outFileName {name}_log2.bw --ratio log2  -p 30 --verbose --outFileFormat bigwig --scaleFactorsMethod readCount'
    cmd_to_run = 'echo {};'.format(x[0])+ cmd.format(ip=x[2], input_f=x[1], name=x[0]) 
    
    os.system(cmd_to_run)
    return cmd_to_run

p =Pool(int(sys.argv[1]))

p.map(run_bamCompare, df.values.tolist())
    