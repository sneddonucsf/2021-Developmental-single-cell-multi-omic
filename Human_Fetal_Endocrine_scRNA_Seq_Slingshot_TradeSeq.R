#Initialize libraries
library(slingshot)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)
library(SingleCellExperiment)
library(RColorBrewer)
library(tradeSeq)
library(BiocParallel)


#### Trajectory and pseudotime inference on endocrine dataset with Slingshot ####
slingshot <- slingshot(Embeddings(HumanFetal_combined_IntegrationAnchors_Object_Endocrine_labeled, "umap"), 
                       clusterLabels = HumanFetal_combined_IntegrationAnchors_Object_Endocrine_labeled$CellfindR, 
                 start.clus = "3.1.0")

#Function to extract the pallete for plotting
cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

#Plot the resulting trajectories
cell_colors_clust <- cell_pal(HumanFetal_combined_IntegrationAnchors_Object_Endocrine_labeled$CellfindR, hue_pal())
plot(reducedDim(slingshot), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(slingshot, lwd = 2, type = 'lineages', col = 'black')
#Smooth the curves
plot(reducedDim(slingshot), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(slingshot, lwd = 2, col = 'black')
dev.off()

#Extract the resulting pseudotime and cellweight values for tradeSeq analysis 
pseudotime <- slingPseudotime(slingshot, na = FALSE) #get the calculated pseudotime 
cellWeights <- slingCurveWeights(slingshot) #get the calculated cellWeights (likelihood it belongs to a particular lineage)
exprs_human<- GetAssayData(HumanFetal_combined_IntegrationAnchors_Object_Endocrine_labeled, slot="counts", assay = "RNA") #get the gene counts

#### Differential gene expression analysis across lineages with TradeSeq ####
#Fit smoothers for all genes 
tradeseq_obj <- fitGAM(counts = exprs_human, pseudotime = pseudotime, cellWeights = cellWeights,
              nknots = 7, verbose = TRUE)

#Perform DE tests across the lineages; lineages = TRUE to measure each lineage separately
assoRes <- associationTest(tradeseq_obj, lineages = TRUE)
startRes <- startVsEndTest(tradeseq_obj, lineages = TRUE)
endRes <- diffEndTest(tradeseq_obj, lineages = TRUE) 
patternRes <- patternTest(tradeseq_obj, lineages = TRUE)
compare <- inner_join(patternRes %>% mutate(Gene = rownames(patternRes),
                                            pattern = waldStat) %>%
                        select(Gene, pattern),
                      endRes %>% mutate(Gene = rownames(endRes),
                                        end = waldStat) %>%
                        select(Gene, end),
                      by = c("Gene" = "Gene")) %>%
  mutate(transientScore = (min_rank(desc(end)))^2 +
           (dense_rank(pattern))^2)



