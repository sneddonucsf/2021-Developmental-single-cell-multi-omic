#load packages
library("DESeq2")
library("limma")
library("edgeR")
library("gdata")
library("gplots")
library("Hmisc")
library("ggrepel")
library("org.Hs.eg.db")
library("EnsDb.Hsapiens.v79")
library("clusterProfiler")
library("DOSE")
library("EnhancedVolcano")
library("viridis")

#Load the counts file 
data <- read.csv("/Users/seandelao/Desktop/readcount_genename.csv", header=T, check.names=F,
                 row.names=1)

#Subset D10 data
colnames(data)
myvars <- c("WT07S6D10", "WT10S6D10", "WT11S6D10", "KO07S6D10", "KO10S6D10", "KO11S6D10")
data <- data[myvars]
nrow(data) #57,905 rows
data = as.data.frame(data)

#Remove rows with zero values to shrink data
keep <-rowSums(cpm(data)>0) >= 6
keep <-rowSums(cpm(data)>0) >=6 #(here ,>=6, 6 means in 3 samples, can be adjusted)
data.keep <-data[keep,]
nrow(data.keep) #20,174
data.keep = as.matrix(data.keep)
data.keep

#Create a sample data column
sampleTable <- data.frame(row.names = colnames(data),
                          Genotype = c("WT", "WT", "WT",  "KO", "KO", "KO"))

#Transform to DESeq matrix
DE_matrix <- DESeqDataSetFromMatrix(countData = data, colData =
                                      sampleTable, design = ~Genotype)

#Log transform
DE_matrix_log <- rlog(DE_matrix)

#Correlation and PCA clustering
exprsR <- assay(DE_matrix_log)
data_corr <- cor(exprsR, method="pearson")
dev.off()
heatmap(data_corr)

#PCA
PCA<-plotPCA(DE_matrix_log, intgroup = "Genotype", returnData=T)
PCA
percentVar <- round(100 * attr(PCA, "percentVar"))
qplot(PC1, PC2, data=PCA) + 
  theme(text = element_text(size=15), axis.text.x = element_text(angle=90, vjust=1)) + 
  geom_point(colour = "black", size = 1, stroke = 1) + 
  theme(legend.title=element_blank()) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  xlim(c(-15,30))+ylim(c(-16,15)) + geom_text_repel(aes(PC1, PC2, label = rownames(PCA))) +
  theme_classic(base_size = 16)

#Differential gene expression
matrix1 <- DESeqDataSetFromMatrix(countData = data, colData =
                                    sampleTable, design = ~ Genotype)
DEmatrix1 <- DESeq(DE_matrix)
toptable1 <- results(DEmatrix1)
summary(toptable1) #gives summary of DESeq analysis
head(toptable1) #see what the data looks like
toptable1

#Significant genes 
toptable1_matrix <- as.data.frame(toptable1, geneid = rownames(toptable1))
nrow(toptable1_matrix) #20,174 genes

#Convert from ensembl.gene to gene.symbol
ensembl.genes <- rownames(toptable1_matrix)
geneIDs1 <- ensembldb::select(EnsDb.Hsapiens.v79, keys= ensembl.genes, keytype = "GENEID", columns = c("SYMBOL","GENEID"))
toptable1_matrix$Ensembl_ID <- rownames(toptable1_matrix)
toptable1_matrix = merge(toptable1_matrix, geneIDs1, by.x="Ensembl_ID",by.y="GENEID", row.names(toptable1_matrix))
toptable1_matrix

# Remove genes that did not get an adjusted p value
toptable1_matrix <- dplyr::filter(toptable1_matrix, padj != 'NA')
toptable1_matrix

#Significant genes based on p-value
significant <- dplyr::filter(toptable1_matrix, padj < 0.05)
significant
nrow(significant) #2,008 genes
#write.csv(significant, file = "/Users/seandelao/Desktop/FEV WT_KO Bulk/S6D10/FEV WT_KO S6D10 DESeq2 Significant.csv")

#Number of up-regulated (higher in WT)
upregulated <- dplyr::filter(significant, log2FoldChange > 0.5)
upregulated <- dplyr::arrange(upregulated, -log2FoldChange)
nrow(upregulated)
write.csv(upregulated, file = "/Users/seandelao/Desktop/FEV WT_KO S6D10 DESeq2 Upregulated in WT log2FC > 0.5.csv")
#Upregulated:964

#Number of downregulated (higher in KO)
downregulated <- dplyr::filter(significant, log2FoldChange < -0.5)
downregulated <- dplyr::arrange (downregulated, desc(-log2FoldChange))
nrow(downregulated)
downregulated
write.csv(downregulated, file = "/Users/seandelao/Desktop/FEV WT_KO S6D10 DESeq2 Upregulated in KO log2FC < -0.5.csv")
#Downregulated: 873

#### Plotting results ####
#Volcano plot of differentially expressed genes 
EnhancedVolcano(toptable1_matrix, lab = toptable1_matrix$SYMBOL, 
                selectLab = "",
                x = 'log2FoldChange', 
                y = 'pvalue', 
                xlim = c(-7, 6), 
                title = 'FEV WT vs KO', 
                pCutoff = 10e-2, 
                FCcutoff = 0.5, 
                ylim = c(0, 31), 
                drawConnectors = TRUE,
                colAlpha = 0.8, 
                boxedLabels = FALSE, 
                col = c('black', 'black', 'black', 'red3'))

#### Heatmap of selected genes ####
#Annotate exprsR results with gene names for plotting 
anno <- select(EnsDb.Hsapiens.v79,keys = as.character(rownames(exprsR)), keytype = "GENEID", columns = c("ENSEMBL", "SYMBOL", "GENENAME"))
pick.first <- anno[match(row.names(exprsR), anno$GENEID),]
pick.first
cbind(exprsR, pick.first)
rownames(exprsR) <- pick.first$SYMBOL

#Subset heatmap for genes you want to plot
Beta = c(
  #Beta
  "INS",
  "DLK1",
  "IAPP",
  "ACVR1C",
  "ASB9",
  "ROBO2",
  "ASPH",
  "SLC30A8",
  "LMO2",
  #FEV High
  'RASGRP1',
  'DCDC1',
  'SSTR1',
  'GAD2',
  'NKX2-2',
  'NEUROD1',
  'RASSF6', 
  #Pre-Beta
  'STMN2',
  'HECTD4',
  'RAB3B',
  'SIM1',
  'KCNK16',
  'MNX1',
  'SCGN',
  'KCNB2',
  'CRYBA2',
  'PAX6')
genes <- subset((exprsR), rownames(exprsR) %in% Beta)
my_palette2 <- colorRampPalette(c("blue4", "white", "red4"))(n=299)
pcorr<-function(x) as.dist(1-cor(t(x),method="pearson"))
dev.off()
heatmap.2(as.matrix(genes), 
          dendrogram="both",col=my_palette2, 
          scale="row", key=T, density.info="none", 
          trace="none", cexCol=1, cexRow=1, Colv=T, 
          Rowv=F, labRow=NULL, distfun=pcorr)

#Alpha
Alpha = c(
  #Beta
  "GCG",
  "TTR",
  "IRX2",
  "CRYBA2",
  "CHGA",
  "C5orf38",
  "GC",
  "SCG5",
  "IRX1",
  'F10',
  'PTPRT',
  'SMC2',
  'SLC7A2',
  'KCNK16',
  #Pre-Alpha
  'NKX2-2',
  'NEUROD1',
  'FGF14',
  'IGFBPL1',
  'CACNA2D1',
  'INSM1',
  'RFX6',
  'STMN2',
  'GLUL',
  'TP53INP1')
genes <- subset((exprsR), rownames(exprsR) %in% Alpha)
heatmap.2(as.matrix(genes), 
          dendrogram="none",col=my_palette2, 
          scale="row", key=T, density.info="none", 
          trace="none", cexCol=1, cexRow=1, Colv=F, 
          Rowv=F, labRow=NULL, distfun=pcorr)

#Delta and Epsilon
Delta_Epsilon = c(
  #Delta
  "SST",
  "CRH",
  "PPP1R1A",
  "SEC11C",
  "C3orf14",
  "AMIGO2",
  "SERPINA1",
  "GAD2",
  "NPW",
  'KCTD8',
  #Epsilon
  'SPTSSB',
  'FGF12',
  'FGF14',
  'CNNM1',
  'ELL2',
  'QDPR',
  'IL20RA',
  'RGS17',
  'NKX2-2',
  'APLP1', 
  "F10", 
  "P4HTM", 
  "GCNT2")
genes <- subset((exprsR), rownames(exprsR) %in% Delta_Epsilon)
heatmap.2(as.matrix(genes), 
          dendrogram="none",col=my_palette2, 
          scale="row", key=T, density.info="none", 
          trace="none", cexCol=1, cexRow=1, Colv=F, 
          Rowv=F, labRow=NULL, distfun=pcorr)

#FEV targets
FEV_targets = c("TPH1", "TPH2", "GCK", "Slc2a2", "LMX1B", "DDC", "SLC6A4", "SLC18A2")
genes <- subset((exprsR), rownames(exprsR) %in% FEV_targets)
heatmap.2(as.matrix(genes), 
          dendrogram="none",col=my_palette2, 
          scale="row", key=T, density.info="none", 
          trace="none", cexCol=1, cexRow=1, Colv=F, 
          Rowv=F, labRow=NULL, distfun=pcorr)

#Top DEGs
WT = head(upregulated$SYMBOL, n = 10)
KO = head(downregulated$SYMBOL, n = 10)
DEGs = c(WT, KO)
genes <- subset((exprsR), rownames(exprsR) %in% DEGs)
heatmap.2(as.matrix(genes), 
          dendrogram="none",col=my_palette2, 
          scale="row", key=T, density.info="none", 
          trace="none", cexCol=1, cexRow=1, Colv=F, 
          Rowv=F, labRow=NULL, distfun=pcorr)

#### Gene ontology analysis ####
df = toptable1
original_gene_list <- df$log2FoldChange
names(original_gene_list) <- rownames(df)
gene_list<-na.omit(original_gene_list)
gene_list = sort(gene_list, decreasing = TRUE)
gse <- gseGO(geneList=gene_list,
             ont ="ALL",
             keyType = "ENSEMBL",
             nPerm = 10000,
             minGSSize = 3,
             maxGSSize = 800,
             pvalueCutoff = 0.05,
             verbose = TRUE,
             OrgDb = org.Hs.eg.db, pAdjustMethod = "none")
dotplot(gse, showCategory=20, split=".sign") + facet_grid(.~.sign)
ridgeplot(gse) + labs(x = "enrichment distribution")
gseaplot(gse, by = "all", title = gse$Description[3], geneSetID = 3)

#### KEGG analysis ####
ids<-bitr(names(original_gene_list), fromType = "ENSEMBL", toType =
            "ENTREZID", OrgDb=org.Hs.eg.db)
dedup_ids = ids[!duplicated(ids[c("ENSEMBL")]),]
df2 = df[rownames(df) %in% dedup_ids$ENSEMBL,]
df2$Y = dedup_ids$ENTREZID
kegg_gene_list <- df2$log2FoldChange
names(kegg_gene_list) <- df2$Y
kegg_gene_list<-na.omit(kegg_gene_list)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)
kegg_organism = "hsa"
kk2 <- gseKEGG(geneList = kegg_gene_list,
               organism = kegg_organism,
               nPerm = 10000,
               minGSSize = 3,
               maxGSSize = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType = "ncbi-geneid")
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" ,
        split=".sign") + facet_grid(.~.sign)
ridgeplot(kk2) + labs(x = "enrichment distribution")