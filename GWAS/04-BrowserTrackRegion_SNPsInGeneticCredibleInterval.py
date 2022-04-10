def BrowserTrackRegion(inF, inF2='Homo_sapiens.GRCh38.104.GeneRegionType', FLANK=10000):
    D = {}
    inFile = open(inF2)
    for line in inFile:
        line = line.strip()
        fields = line.split('\t')
        gene = fields[1]
        if gene not in D:
            D[gene] = fields
    inFile.close()

    inFile = open(inF)
    ouFile = open(inF.split('.txt')[0] + '_TrackRegion.txt', 'w')
    ouFile.write('\t'.join(['gene', 'rs', 'ch', 'region_start', 'region_end', 'index_SNP', 'CredibleSNPsPeaks', 'CredibleSNPsAll']) + '\n')
    head = inFile.readline()
    for line in inFile:
        line = line.strip()
        fields = line.split('\t')
        index_SNP_ch = fields[2]
        rs = fields[1]

        index_SNP = fields[17]
        fds = fields[-2].split('#')
        SNPsWithPeak = [x.split('|')[0].split(':')[1] for x in fds]
        print(SNPsWithPeak)

        fds2 = fields[-3].split('|')
        SNPsAll = [x.split(':')[1] for x in fds2]
        print(SNPsAll)

        gene = fields[0].split('/')[0]
        if gene in D:
            gene_ch = D[gene][2]
            gene_start = D[gene][3]
            gene_end = D[gene][4]

            if index_SNP_ch == gene_ch:
                L = [int(x) for x in [gene_start, gene_end] + SNPsAll + [index_SNP]]
                left = str(min(L) - FLANK)
                right = str(max(L) + FLANK)
                ouFile.write('\t'.join([gene, rs, index_SNP_ch, left, right, index_SNP, ','.join(SNPsWithPeak), ','.join(SNPsAll)]) + '\n')


    inFile.close()
    ouFile.close()



BrowserTrackRegion('Adult_Beta-vs-Fetal_Beta_FDR0.05Log2FC1_down_T2D_SNPsInGeneticCredibleInterval_ByLoci.txt')
BrowserTrackRegion('Adult_Hormone-vs-Fetal_EPs_FDR0.05Log2FC1_down_T2D_SNPsInGeneticCredibleInterval_ByLoci.txt')

