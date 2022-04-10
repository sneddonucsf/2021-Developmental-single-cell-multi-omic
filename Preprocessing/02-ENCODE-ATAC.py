import os
import json

def PerBioRep(inFs, WDL='atac-seq-pipeline/atac.wdl', TSV='atac-seq-pipeline/data/hg38/hg38.tsv'):
        ouFile = open('CB.sh', 'w')
        ouFile.write('#!/usr/bin/bash\n')
        ouFile.write('cd CB\n')
        ouFile.write('source activate encode-atac-seq-pipeline\n')
        ouFile.write('caper run %s -i %s\n'%(WDL, '../' +  'CB.json'))
        ouFile.close()

        CONFIG = {}
        CONFIG['atac.pipeline_type'] = 'atac'
        CONFIG['atac.genome_tsv'] = TSV
        for n in range(len(inFs)):
            CONFIG['atac.fastqs_rep%s_R1'%(n+1)] = ['%s_1.fastq.gz'%inFs[n]]
            CONFIG['atac.fastqs_rep%s_R2'%(n+1)] = ['%s_2.fastq.gz'%inFs[n]]

        CONFIG['atac.paired_end'] = True
        CONFIG['atac.auto_detect_adapter'] = True
        CONFIG['atac.enable_xcor'] = False
        CONFIG['atac.title'] = 'snATAC-Seq'

        ouFile = open('CB.json', 'w')
        json.dump(CONFIG, ouFile, indent = 4)
        ouFile.close()

Fs = set([os.getcwd() + '/' + x.split('_1')[0].split('_2')[0] for x in os.listdir('.') if x.endswith('.gz') and x.find('10X') == -1])
inFs = sorted(list(Fs))
PerBioRep(inFs)
