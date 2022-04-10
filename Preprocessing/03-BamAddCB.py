import pysam
import gzip

def BamAddCB(inF):
    D = {}
    fq1 = '../' + inF.split('_1.trim')[0] + '_1.fastq.gz'
    inFile = gzip.open(fq1)
    while True:
        line1 = inFile.readline()
        line2 = inFile.readline()
        line3 = inFile.readline()
        line4 = inFile.readline()
        if line1:
            fds = line1.decode().split()
            qname = fds[0][1:]
            cb = fds[1].split(':')[0]
            D[qname] = cb
        else:
            break
    inFile.close()


    inFile = pysam.AlignmentFile(inF)
    ouFile = pysam.AlignmentFile(inF.split('.bam')[0] + '_CB.bam', 'w', template=inFile)

    cb = 'NA'
    for item in inFile.fetch():
        if item.qname in D:
            cb = D[item.qname]
        item.set_tag('CB', cb)
        ouFile.write(item)



BamAddCB('Islet1_1.trim.srt.bam')
BamAddCB('Islet2_1.trim.srt.bam')
BamAddCB('Islet3_1.trim.srt.bam')
