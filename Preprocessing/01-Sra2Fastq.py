import pandas as pd
import os

def Sra2Fastq(inF, inF2=''):
    D = {}
    if inF2:
        inFile = open(inF2)
        for line in inFile:
            line = line.strip()
            fields = line.split('\t')
            sample = ''.join(fields[1].split())
            D[fields[0]] = sample
        inFile.close()

    ouFile = open(inF.split('.txt')[0] + '.sh', 'w')
    ouFile.write('#!/usr/bin/bash\n\n')
    df = pd.read_table(inF, header=0, sep=',')
    for n in range(df.shape[0]):
        SRR = df.loc[n, 'Run']
        SRP = df.loc[n, 'SRA Study']
        GSM = df.loc[n, 'Sample Name']
        LL = df.loc[n, 'LibraryLayout']
        #CL = df.loc[n, 'Cell_Line']
        AssayType = df.loc[n, 'Assay Type']
        print(GSM)
        if GSM in D:
            sample = D[GSM]
            ouFile.write('wget `srapath %s` -O %s.sra\n'%(SRR, sample))
            if LL == 'PAIRED':
                ouFile.write('fastq-dump --split-files --gzip %s.sra\n'%sample)
            else:
                ouFile.write('fastq-dump --gzip %s.sra\n'%sample)

            ouFile.write('rm %s.sra\n'%sample)
        else:
            sample = GSM + '_' + SRR 
            ouFile.write('wget `srapath %s` -O %s.sra\n'%(SRR, sample))
            if LL == 'PAIRED':
                ouFile.write('fastq-dump --split-files --gzip %s.sra\n'%sample)
            else:
                ouFile.write('fastq-dump --gzip %s.sra\n'%sample)

            ouFile.write('rm %s.sra\n'%sample)

    ouFile.close()


Sra2Fastq('SraRunTable.txt', 'SampleAnnot.txt')
