#Initialize libraries
library(Seurat)
library(dplyr)
library(SCENIC)
BiocManager::install(c("AUCell", "RcisTarget", "SingleCellExperiment"), ask = FALSE)

#### First part of SCENIC analysis ####
#Get raw data from endocrine dataset 
HumanFetal_combined_IntegrationAnchors_Object_Endocrine_labeled.raw.data <- GetAssayData(object = HumanFetal_combined_IntegrationAnchors_Object_Endocrine_labeled, 
                                                                                               assay = "RNA", slot = "data")

#Extract cellInfo for plotting, etc. 
cellInfo <- colnames(HumanFetal_combined_IntegrationAnchors_Object_Endocrine_labeled.raw.data)
cellInfo <- as.data.frame(cbind(cellInfo))
cellInfo = cbind(cellInfo, HumanFetal_combined_IntegrationAnchors_Object_Endocrine_labeled.raw.data@active.ident, 
                 HumanFetal_combined_IntegrationAnchors_Object_Endocrine_labeled.raw.data@meta.data$nFeature_RNA, 
                 HumanFetal_combined_IntegrationAnchors_Object_Endocrine_labeled.raw.data@meta.data$nCount_RNA)
colnames(cellInfo) <- c("rowname", "Class", "nGene", "nUMI")
cellInfo <- cellInfo[, 2:4]
exprMat <- as.matrix(HumanFetal_combined_IntegrationAnchors_Object_Endocrine_labeled.raw.data)

#Create colVars
DimPlot(HumanFetal_combined_IntegrationAnchors_Object_Endocrine_labeled)
colVars <- list(CellType=c("0"="#F8766D", 
                           "1"="#D39200", 
                           "3.0.0"="#00BA38", 
                           "3.0.1"="#00C19F", 
                           "3.1.0"="#00B9E3", 
                           "3.1.1"="#619CFF", 
                           "4"="#DB72FB", 
                           "5"="#FF61C3"))
colVars$CellType <- colVars$CellType[intersect(names(colVars$CellType), cellInfo$CellType)]

#Create scenicOptions
org = "hgnc" # or hgnc, or dmel
dbDir = directory_to_RcisTarget_database
myDatasetTitle = "SCENIC sc lib"
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 
scenicOptions@inputDatasetInfo$cellInfo <- cellInfo
scenicOptions@inputDatasetInfo$colVars <- colVars

#Filter gene matrix 
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)
exprMat_filtered <- exprMat[genesKept, ]

#run correlation
runCorrelation(exprMat_filtered, scenicOptions)

#Export filtered matrix for Arboreto step
exportsForArboreto(exprMat_filtered, scenicOptions, dir = output_directory)

#### run this next section in command line, NOT in R ####
#import os
#import pandas as pd
#from arboreto.algo import grnboost2, genie3
#from arboreto.utils import load_tf_names
#if __name__ == '__main__':
#  ex_matrix = pd.read_csv('output_directory/1.1_exprMatrix_filtered_t.txt', sep='\t')
#tf_names = load_tf_names('output_directory/1.1_inputTFs.txt')
#network = grnboost2(expression_data=ex_matrix,
#                    tf_names=tf_names)
#network.to_csv('output_directory/network.tsv', sep='\t', header=False, index=False)

#### Part 3 of SCENIC, this is run in R ####
#Load SCENIC Options object 
GRNBoost_output <- read.delim("output_directory/sc_01_network.tsv", header=FALSE)
colnames(GRNBoost_output) <- c("TF","Target","weight")
runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
exprMat <- as.matrix(HumanFetal_combined_IntegrationAnchors_Object_Endocrine_labeled.raw.data)
logMat <- log2(exprMat+1)
runSCENIC_3_scoreCells(scenicOptions, logMat)

#### Analyze SCENIC results ####
#scaled relative regulon activity by cluster - heatmap - all regulons
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType_human_endocrine <- sapply(split(rownames(cellInfo), cellInfo$seuratCluster),
                                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_human_endocrine <- regulonActivity_byCellType_human_endocrine
regulonActivity_byCellType_human_endocrine_scaled <- t(scale(t(regulonActivity_byCellType_human_endocrine), center = T, scale=T))
regulonActivity_byCellType_human_endocrine_scaled <- regulonActivity_byCellType_human_endocrine_scaled

#top regulons by cluster of scaled data
topRegulators <- reshape2::melt(regulonActivity_byCellType_human_endocrine_scaled)
colnames(topRegulators) <- c("Regulon", "seuratClusters", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]

