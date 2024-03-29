
###########################################################################################################
### This script details how to take the 10x snATAC Seq output and analyze with ArchR. The dataset       ###
### generated is shown in figure 3 of the paper.                                                        ###
###########################################################################################################

library(ArchR)
library(cowplot)
library(Seurat)
library(viridis)
library(GenomicRanges)
library(stringr)

addArchRThreads()

#### Create the arrow file from the fragments file generated by CellRanger ####
inputFiles = directory_to_fragments_file
addArchRGenome("hg38")
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 4, #Dont set this too high because you can always increase later
  filterFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

#### Create dataset and filter ####
#Load 12wpc_Epcam Arrow File to create Arrow Project 
outdir = output_directory
Twelve_Epcam <- ArchRProject(ArrowFiles = ArrowFiles, outputDirectory = outdir, copyArrows = FALSE)
Twelve_Epcam
#6,429 cells total 

#Perform cutoffs and do downstream analysis
plotGroups(
  ArchRProj = Twelve_Epcam, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

plotGroups( 
  ArchRProj = Twelve_Epcam, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "nFrags",
  plotAs = "violin",
  alpha = 0.4,
  addBoxPlot = TRUE
)

#Perform the analysis with doublet filtering
doubScores <- addDoubletScores(
  input = Twelve_Epcam,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search with doublet projection.
  LSIMethod = 1)
#UMAP Projection R^2 = 0.99742

Twelve_Epcam <- filterDoublets(Twelve_Epcam)
#Filtering 413 cells from ArchRProject!
#Epcam-12wpc : 413 of 6429 (6.4%)

Twelve_Epcam
#6016 cells total 

#### Dimensional reduction, clustering and UMAP ####
#Dimensional reduction
Twelve_Epcam <- addIterativeLSI(
  ArchRProj = Twelve_Epcam,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.5), 
    sampleCells = 2500, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30, 
  force = TRUE
)

#LSI clustering
Twelve_Epcam <- addClusters(
  input = Twelve_Epcam,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.1, 
  force = TRUE
)

#UMAP 
Twelve_Epcam <- addUMAP(
  ArchRProj = Twelve_Epcam, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine"
)

#### Subset endocrine cells and perform dimensional reduction, clustering and UMAP ####
idxPass <- which(Twelve_Epcam@cellColData$Clusters %in% c("C7", "C8", "C9"))
cellsPass <- Twelve_Epcam$cellNames[idxPass]
Twelve_Epcam_Endocrine = Twelve_Epcam[cellsPass, ]
Twelve_Epcam_Endocrine #1730 cells; ~4k cells filtered

#Dimensional reduction
Twelve_Epcam_Endocrine <- addIterativeLSI(
  ArchRProj = Twelve_Epcam_Endocrine,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.3), 
    sampleCells = 1500, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 2:10, 
  force = TRUE
)

#LSI clustering
Twelve_Epcam_Endocrine <- addClusters(
  input = Twelve_Epcam_Endocrine,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 1, 
  force = TRUE
)

#UMAP 
Twelve_Epcam_Endocrine <- addUMAP(
  ArchRProj = Twelve_Epcam_Endocrine, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine", 
  dimsToUse = 2:10, 
  force = TRUE
)

#### Perform unconstrained integration and pseudo-scRNA-Seq profile ####
library(Seurat)
load(directory_to_endocrine_dataset_from_figure_2)
seRNA <- SetIdent(seRNA, value = 'integrated_snn_res.0.1') #This is the tier 1 clustering from CellFindR. 
#We chose this level of clustering due to the low cell numbers of the heterogenous EP populations if we use the Tier 3 
#clustering from the scRNA Endocrine dataset. Essentially we collapsed all four EP populations into 1. 

Twelve_Epcam_Endocrine <- addGeneIntegrationMatrix(
  ArchRProj = Twelve_Epcam_Endocrine, 
  useMatrix = "GeneScoreMatrix",
  matrixName = "GeneIntegrationMatrix",
  reducedDims = "IterativeLSI",
  seRNA = seRNA,
  addToArrow = FALSE,
  force = TRUE,
  groupRNA = "integrated_snn_res.0.1",
  nameCell = "predictedCell_Un",
  nameGroup = "predictedGroup_Un",
  nameScore = "predictedScore_Un")

