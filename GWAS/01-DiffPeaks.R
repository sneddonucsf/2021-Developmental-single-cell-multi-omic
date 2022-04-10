library(ArchR)

set.seed(1234)

proj = loadArchRProject(path='snATACseq')

ccd = proj@cellColData
ccd$CellTypeMerged = gsub('_\\d+', '', ccd$CellType)
proj@cellColData = ccd

p1 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "Sample", embedding = "UMAP", rastr=F, baseSize=40)
p2 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "CellType", embedding = "UMAP", rastr=F, baseSize=40)
p3 <- plotEmbedding(ArchRProj = proj, colorBy = "cellColData", name = "CellTypeMerged", embedding = "UMAP", rastr=F, baseSize=40)
plotPDF(p1, p2, p3, name = "Plot-UMAP-Sample-CellType.pdf", ArchRProj = proj, addDOC = FALSE, width = 20, height = 20)

proj = addGroupCoverages(ArchRProj = proj, groupBy = "CellTypeMerged", force=TRUE)

proj <- addReproduciblePeakSet(ArchRProj = proj, groupBy = "CellTypeMerged", force=TRUE)

proj = addPeakMatrix(proj)
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", name = "Motif", force=TRUE)

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
p1 <- plotMarkers(seMarker = mT, name = g1, cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
p2 <- plotMarkers(seMarker = mT, name = g1, cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "Volcano")
plotPDF(p1, p2, name = sprintf("%s-vs-%s-Markers-MA-Volcano", paste0(g2, collapse='_'), g1), width = 20, height = 20, ArchRProj = proj, addDOC = FALSE)
}else{
p1 <- plotMarkers(seMarker = mT, name = g1, cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "MA")
p2 <- plotMarkers(seMarker = mT, name = g1, cutOff = "FDR <= 0.05 & abs(Log2FC) >= 1", plotAs = "Volcano")
plotPDF(p1, p2, name = sprintf("%s-vs-%s-Markers-MA-Volcano", g1, paste0(g2, collapse='_')), width = 20, height = 20, ArchRProj = proj, addDOC = FALSE)

}

m = assays(mT)
mdf = c(rowData(mT), m$AUC, m$FDR, m$Log2FC, m$Mean, m$MeanBGD, m$MeanDiff, m$Pval)
colnames(mdf) = c('seqnames', 'idx', 'start', 'end', 'AUC', 'FDR', 'Log2FC', 'Mean', 'MeanBGD', 'MeanDiff', 'Pval')
mdf = mdf[order(mdf$Pval), ]

if(grepl('Fetal', g1)){
ouF1 = sprintf('%s-vs-%s.txt', paste0(g2, collapse='_'), g1)
ouF2 = sprintf('%s-vs-%s_FDR0.05Log2FC1.txt', paste0(g2, collapse='_'), g1)
}else{
ouF1 = sprintf('%s-vs-%s.txt', g1, paste0(g2, collapse='_'))
ouF2 = sprintf('%s-vs-%s_FDR0.05Log2FC1.txt', g1, paste0(g2, collapse='_'))
}
write.table(mdf, ouF1, row.names=F, col.names=T, sep='\t', quote=F)
wh = mdf$FDR <= 0.05 & abs(mdf$Log2FC) >= 1
mdf_sub = mdf[wh, ]
write.table(mdf_sub, ouF2, row.names=F, col.names=T, sep='\t', quote=F)
return(mT)
}

mT1 = DiffPeaks('Adult_beta', 'Fetal_Beta')
mT2 = DiffPeaks('Adult_alpha', 'Fetal_Alpha')
mT3 = DiffPeaks('Adult_delta', 'Fetal_Delta')
mT4 = DiffPeaks('Fetal_EPs', c('Adult_alpha', 'Adult_beta', 'Adult_delta', 'Adult_gamma', 'Fetal_Alpha', 'Fetal_Beta', 'Fetal_Delta', 'Fetal_Epsilon'))
mT5 = DiffPeaks('Fetal_EPs', c('Adult_alpha', 'Adult_beta', 'Adult_delta', 'Adult_gamma'))

if(FALSE){

MotifEnrich = function(mT, ouF){
motifsUp <- peakAnnoEnrichment(
    seMarker = mT,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.05 & Log2FC >= 1"
  )

mfa = assays(motifsUp)
mdf = cbind(mfa$mlog10Padj, mfa$mlog10p, mfa$Enrichment, mfa$BackgroundProporition, mfa$nBackground, mfa$BackgroundFrequency, mfa$CompareProportion, mfa$nCompare, mfa$CompareFrequency, mfa$feature)
names(mdf) = c("mlog10Padj", "mlog10p", "Enrichment", "BackgroundProporition", "nBackground", "BackgroundFrequency", "CompareProportion", "nCompare", "CompareFrequency", "feature")
mdf = mdf[order(mdf$mlog10Padj, decreasing = TRUE), ]
write.table(mdf, paste0(ouF, '_motifUp.txt'), row.names=T, col.names=NA, sep='\t', quote=F)

dfUp <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
dfUp <- dfUp[order(dfUp$mlog10Padj, decreasing = TRUE),]
dfUp$rank <- seq_len(nrow(dfUp))

ggUp <- ggplot(dfUp, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = dfUp[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched in up peaks") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

motifsDo <- peakAnnoEnrichment(
    seMarker = mT,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.05 & Log2FC <= -1"
  )

mfa = assays(motifsDo)
mdf = cbind(mfa$mlog10Padj, mfa$mlog10p, mfa$Enrichment, mfa$BackgroundProporition, mfa$nBackground, mfa$BackgroundFrequency, mfa$CompareProportion, mfa$nCompare, mfa$CompareFrequency, mfa$feature)
names(mdf) = c("mlog10Padj", "mlog10p", "Enrichment", "BackgroundProporition", "nBackground", "BackgroundFrequency", "CompareProportion", "nCompare", "CompareFrequency", "feature")
mdf = mdf[order(mdf$mlog10Padj, decreasing = TRUE), ]
write.table(mdf, paste0(ouF, '_motifDown.txt'), row.names=T, col.names=NA, sep='\t', quote=F)


dfDo <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
dfDo <- dfDo[order(dfDo$mlog10Padj, decreasing = TRUE),]
dfDo$rank <- seq_len(nrow(dfDo))

ggDo <- ggplot(dfDo, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = dfDo[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched in down peaks") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))

plotPDF(ggUp, ggDo, name = ouF, width = 20, height = 20, ArchRProj = proj, addDOC = FALSE)

}

MotifEnrich(mT1, 'Adult_beta-vs-Fetal_Beta')
MotifEnrich(mT2, 'Adult_alpha-vs-Fetal_Alpha')
MotifEnrich(mT3, 'Adult_delta-vs-Fetal_Delta')
MotifEnrich(mT4, 'AdultFetusHormone-vs-Fetal_EP')
MotifEnrich(mT5, 'AdultHormone-vs-Fetal_EP')

}

proj <- saveArchRProject(ArchRProj = proj)
