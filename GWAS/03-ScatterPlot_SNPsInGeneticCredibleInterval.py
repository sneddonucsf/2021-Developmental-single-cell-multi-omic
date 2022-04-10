from pyst import *

def ScatterPlot(inF, Title, yTickSize=8):
    inFile = open(inF)
    head = inFile.readline()
    L = []
    for line in inFile:
        line = line.strip()
        fields = line.split('\t')
        L.append(fields[0:2] + [int(fields[-4])] + [int(fields[-1])])
    inFile.close()

    df = pd.DataFrame(L)
    df.columns = ['Gene', 'dbSNP', 'No. variants', 'NumSNPsInGeneticCredibleIntervalWithPeaks']
    df['Percentage of variants overlapping with differential peaks'] = df['NumSNPsInGeneticCredibleIntervalWithPeaks']/df['No. variants']*100.0
    df['Loci'] = df['dbSNP'] + '_' + df['Gene']
    df.sort_values(by=['Percentage of variants overlapping with differential peaks', 'Gene'], ascending=False, inplace=True)
    print(df.head())
    print('Percentage > 0: %s loci'%df.shape[0])

    df_sub = df.loc[df['Percentage of variants overlapping with differential peaks'] > 5, ]
    print('Percentage > 5: %s loci'%df_sub.shape[0])

    df_sub = df.loc[df['Percentage of variants overlapping with differential peaks'] > 10, ]
    print('Percentage > 10: %s loci'%df_sub.shape[0])


    df = df.loc[df['Percentage of variants overlapping with differential peaks'] > 20, ]
    print('Percentage > 20: %s loci'%df.shape[0])

    fig = plt.figure()
    ax = fig.add_subplot()
    sns.scatterplot(y='Loci', x='Percentage of variants overlapping with differential peaks', ax=ax, data=df, size='No. variants')
    ax.set_title(Title)
    ax.set_ylabel('')
    plt.yticks(fontsize=yTickSize)
    plt.legend(loc='lower right', title='No. variants')
    plt.tight_layout()
    plt.savefig(inF.split('.txt')[0] + '_ScatterPlot.pdf')




ScatterPlot('Adult_Beta-vs-Fetal_Beta_FDR0.05Log2FC1_down_T2D_SNPsInGeneticCredibleInterval_ByLoci.txt', Title='T2D-risk variants overlapping with chromatin-accessible peaks\nfrom fetal Beta vs. adult Beta cells', yTickSize=6)
ScatterPlot('Adult_Hormone-vs-Fetal_EPs_FDR0.05Log2FC1_down_T2D_SNPsInGeneticCredibleInterval_ByLoci.txt', Title='T2D-risk variants overlapping with chromatin-accessible peaks\nfrom fetal EPs vs. adult Hormone+ endocrine cells', yTickSize=8)
