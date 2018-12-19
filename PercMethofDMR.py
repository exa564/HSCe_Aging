#!/usr/bin/env python

'''
For a directory containing myCpG/methCall files, this script will find the average percent methylation for DMRs (inputed as a seperate bedfile), by doing the following:
1) Preparing an empty table indexed by the DMRs (which are inputed as a bed-3 file without a header)
2) Bedtools intersect of myCpG/methCall and DMR file to obtain a table of DMC that overlap DMR
3) Find percent methylation for each CpG that overlaps a DMR
3) Find average percentage methylation for each DMR
4) Create a dataframe that contains the average percent meth for each sample. Each row is a DMR, and each column is a sample.
if there is no overlap, a NA value will be printed.
'''

from __future__ import division, print_function
import sys
import pybedtools
import glob
import os
import pandas as pd
import argparse
from collections import Counter
import numpy as np
from datetime import datetime



parser = argparse.ArgumentParser(description='This program finds the average percentage of methylation of myCpG files for DMRs')
parser.add_argument('-i', '--inputDirectory', help='Directory that contains myCpG/methCall files', required=True)
parser.add_argument('-d', '--DMRFile', help='Bedfile that contains chr, start, and end for DMR')

args = parser.parse_args()

### Convert DMR to a bedtool object and form empty dataframe thats index is the DMR coordinates
DMRs=pd.read_csv(args.DMRFile, header=None, index_col=None, sep="\t")
DMRs.columns = ['chrom','start','end']
DMRs.index = DMRs.chrom + '-' + DMRs.start.astype(str) + '-' + DMRs.end.astype(str)
FinalDF = pd.DataFrame(index=DMRs.index)
a = pybedtools.BedTool.from_dataframe(DMRs)

### Find average percent methylation for each DMR

path = args.inputDirectory
listing =os.listdir(path)

for methCall in listing:
    fullpath= os.path.join(path, methCall)
    outputName=methCall.split('.')[0]
    with open(fullpath,'r') as b:
        myDF = pd.read_table(b)
    myDF=myDF[['chr', 'base', 'base', 'strand', 'coverage', 'freqC', 'freqT']]
    d=pybedtools.BedTool.from_dataframe(myDF) ##this is adding decimal places to freqC and freqT
    x=a.intersect(d, wa=True, wb=True)
    ovlpDF=pybedtools.BedTool.to_dataframe(x) # freqC and freqT-columns 8/9 are float64, need to convert to int
   
    #Find percent meth for each CpG
    ovlpDF.iloc[:,8]=ovlpDF.iloc[:,8].astype(int)
    ovlpDF.iloc[:,9]=ovlpDF.iloc[:,9].astype(int)
    numCs=pd.Series(ovlpDF.iloc[:,7]*ovlpDF.iloc[:,8]/100, name='numCs').astype(int)
    numTs=pd.Series(ovlpDF.iloc[:,7]*ovlpDF.iloc[:,9]/100, name='numTs').astype(int)
    percMeth=pd.Series((numCs/(numCs+numTs)) *100, name="PercMeth").astype(int)
    percMethDF=pd.concat([ovlpDF.iloc[:,0:3], numCs, numTs, percMeth], axis=1).reset_index()
    
    #Take mean percent meth for each DMR:
    avgMeth=percMethDF.groupby(['chrom', 'start', 'end'])['PercMeth'].mean().reset_index()
    avgMeth.index=avgMeth.chrom + '-' + avgMeth.start.astype(str) + '-' + avgMeth.end.astype(str) 
    avgMeth.rename(columns={'PercMeth':outputName}, inplace=True)
    FinalDF=FinalDF.join(avgMeth.iloc[:,3], how='left')

   
### Export final matrix that contains percent methylation
FinalDF['ID']=FinalDF.index
FinalDF[['chr', 'start', 'end']]=FinalDF['ID'].str.split('-', expand=True)
cols=FinalDF.columns.tolist()
cols=cols[-3:] + cols[:-3]
FinalDF=FinalDF.loc[:,cols]
FinalDF.drop(['ID'], axis=1, inplace=True)
FinalDF.to_csv('PercMethMatrix.txt', sep='\t', na_rep='NA', index=False)




    


    

    



