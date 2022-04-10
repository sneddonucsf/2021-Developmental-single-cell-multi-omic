library(ArchR)
library(ggplot2)
set.seed(1234)

addArchRThreads(threads = 8) 
inputFiles = c('Islet1_1.trim.srt_CB_MAPQ30.bam', 'Islet2_1.trim.srt_CB_MAPQ30.bam', 'Islet3_1.trim.srt_CB_MAPQ30.bam')
names(inputFiles) = c('Islet1', 'Islet2', 'Islet3')

addArchRGenome("hg38")

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, #Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
  bcTag = "CB",
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "CB",
  copyArrows = FALSE 
)

proj <- saveArchRProject(ArchRProj = proj)


