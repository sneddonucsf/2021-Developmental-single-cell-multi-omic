import pandas as pd

def UpDown(inF):
    ouF1 = inF.split('.txt')[0] + '_up.txt'
    ouF2 = inF.split('.txt')[0] + '_down.txt'

    df = pd.read_table(inF, header=0, sep='\t')
    df_up = df.loc[df['Log2FC'] > 0, ]
    df_down = df.loc[df['Log2FC'] < 0, ]

    df_up.to_csv(ouF1, header=True, index=False, sep='\t')
    df_down.to_csv(ouF2, header=True, index=False, sep='\t')

UpDown('Adult_Beta-vs-Fetal_Beta_FDR0.05Log2FC1.txt')
UpDown('Adult_Hormone-vs-Fetal_EPs_FDR0.05Log2FC1.txt')

def Overlap(inF1, inF2='T2D_AssociatedSignalsCredibleSets_NoCommas_GRCh38_LD_r2-0.8_SNPsInGeneticCredibleInterval.txt'):
    D = {}
    inFile = open(inF1)
    head1 = inFile.readline().strip().split('\t')
    for line in inFile:
        line = line.strip('\n')
        fields = line.split('\t')
        ch = fields[0]
        D.setdefault(ch, [])
        D[ch].append('|'.join(fields))
    inFile.close()

    ouFile = open(inF1.split('.txt')[0] + '_T2D_SNPsInGeneticCredibleInterval_ByLoci.txt', 'w')
        
    inFile = open(inF2)
    head2 = inFile.readline().strip()
    ouFile.write(head2 + '\t' + 'SNPsInGeneticCredibleIntervalWithPeak\tNumSNPsInGeneticCredibleIntervalWithPeak' + '\n')
    Lx = []
    for line in inFile:
        line = line.strip()
        fields = line.split('\t')
        ch = 'chr' + fields[2]
        L = []
        Ps = [x for x in fields[-1].split('|') if x.find('NA') == -1]
        
        if ch in D:
            for rs in Ps:
                pos = int(rs.split(':')[1])
                for item in D[ch]:
                    start2 = int(item.split('|')[2])
                    end2 = int(item.split('|')[3])
                    if start2 < pos < end2:
                        L.append(rs + '|' + item)
        if L:
            Lx.append([line] + ['#'.join(L)] + [str(len(L))])

    Lx = sorted(Lx, key=lambda x:float(x[-1]), reverse=True)

    for item in Lx:
        ouFile.write('\t'.join(item) + '\n')


    inFile.close()
    ouFile.close()



Overlap('Adult_Beta-vs-Fetal_Beta_FDR0.05Log2FC1_down.txt')
Overlap('Adult_Hormone-vs-Fetal_EPs_FDR0.05Log2FC1_down.txt')
