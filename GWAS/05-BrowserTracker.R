library(ArchR)
set.seed(1234)

proj = loadArchRProject(path='snATACseq')

DiffPeaks = function(g1='', g2=c(), gby='CellTypeMerged'){
mT <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = gby,
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = g1,
  bgdGroups = g2,
)

if(grepl('Fetal', g1)){
### should be adult vs fetus, do the switching
assays(mT)$Log2FC = - assays(mT)$Log2FC
}
return(mT)
}

mT1 = DiffPeaks('Fetal_Beta', 'Adult_beta')
mT5 = DiffPeaks('Fetal_EPs', c('Adult_alpha', 'Adult_beta', 'Adult_delta', 'Adult_gamma'))

BrowserTrack = function(mT, g1, g2, gby='CellTypeMerged', inF=''){

sig_peaks_df = read.table(gsub('_down_T2D_SNPsInGeneticCredibleInterval_ByLoci_TrackRegion', '', inF), header=T, sep='\t')
sig_peaks = GRanges(seqnames=sig_peaks_df$seqnames, ranges = IRanges(sig_peaks_df$start, sig_peaks_df$end))

mat = read.table(inF, header=T)

for(n in 1:nrow(mat)){

    gene = as.character(mat[n, 'gene'])
    rs = as.character(mat[n, 'rs'])
    ch = paste0('chr', as.character(mat[n, 'ch']))
    region_start = mat[n, 'region_start']
    region_end = mat[n, 'region_end']
    loci_points = strsplit(as.character(mat[n, 'CredibleSNPsPeaks']), ',')[[1]]
    snp_points = strsplit(as.character(mat[n, 'CredibleSNPsAll']), ',')[[1]]

    loci = GRanges(seqnames=Rle(ch, length(loci_points)), ranges=IRanges(as.integer(loci_points) - 100, as.integer(loci_points) + 100))
    snp = GRanges(seqnames=Rle(ch, length(snp_points)), ranges=IRanges(as.integer(snp_points) - 100, as.integer(snp_points) + 100))

    p <- plotBrowserTrack(
    ArchRProj = proj,
    region = GRanges(ch, region_start:region_end),
    groupBy = gby,
    useGroups = c(g1, g2),
    features =  GRangesList(Differential_Peaks=sig_peaks, T2D_SNPsInCredibleInterval=snp, T2D_SNPsInCredibleIntervalWithPeaks=loci),
    title = paste0(gene, ' ', rs),
    baseSize = 40,
    )

    if(length(g2) > 1){
        g3 = 'Adult_Hormone'
        plotPDF(p, name = sprintf("BrowserTrack_%s-vs-%s_%s_%s.pdf",g3, g1, gene, rs), ArchRProj = proj, addDOC = FALSE, width = 20, height = 20)
        }else{
        plotPDF(p, name = sprintf("BrowserTrack_%s-vs-%s_%s_%s.pdf",g2, g1, gene, rs), ArchRProj = proj, addDOC = FALSE, width = 20, height = 20)
        }

}
}

BrowserTrack(mT1, g1='Fetal_Beta', g2='Adult_beta', inF='Adult_Beta-vs-Fetal_Beta_FDR0.05Log2FC1_down_T2D_SNPsInGeneticCredibleInterval_ByLoci_TrackRegion.txt')
BrowserTrack(mT5, g1='Fetal_EPs', g2=c('Adult_alpha', 'Adult_beta', 'Adult_delta', 'Adult_gamma'), inF='Adult_Hormone-vs-Fetal_EPs_FDR0.05Log2FC1_down_T2D_SNPsInGeneticCredibleInterval_ByLoci_TrackRegion.txt')

