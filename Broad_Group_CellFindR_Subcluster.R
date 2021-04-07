
###########################################################################################################
### This script details the cleaning up of the dataset and and running CellFindR clustering on each     ###
### Broad Group.                                                                                        ###
###########################################################################################################

library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggpubr)

#### Remove clusters that seem to be empty droplets ####
HumanFetal_combined_IntegrationAnchors_Object <- subset(HumanFetal_combined_IntegrationAnchors_Object, 
                                                              idents = c("26", "50"), invert = TRUE)

#### Identify each Broad Group by marker expression and subset ####
#Mesenchymal cells (COL3A1+)
HumanFetal_combined_IntegrationAnchors_Object_Mesenchymal <- subset(
  HumanFetal_combined_IntegrationAnchors_Object, 
  idents = c(1,2,5,6,7,8,9,10,
             13,18,19,25,28,29,
             30,33,34,39,40,43,62))

#Neuronal cells (SOX10+)
HumanFetal_combined_IntegrationAnchors_Object_Neuronal <- subset(
  HumanFetal_combined_IntegrationAnchors_Object, 
  idents = c(21, 58, 54), invert = FALSE)

#Exocrine cells (CPA1/CFTR+)
HumanFetal_combined_IntegrationAnchors_Object_Exocrine <- subset(
  HumanFetal_combined_IntegrationAnchors_Object, 
  idents = c(4, 59, 27, 3, 23, 12, 15, 44, 56, 24, 51, 53), invert = FALSE)

#Endocrine cells (CHGA+)
HumanFetal_combined_IntegrationAnchors_Object_Endocrine <- subset(
  HumanFetal_combined_IntegrationAnchors_Object, 
  idents = c(49, 48, 35, 36, 55, 17, 42, 0, 20, 57), invert = FALSE)

#Immune cells
HumanFetal_combined_IntegrationAnchors_Object_Immune <- subset(
  HumanFetal_combined_IntegrationAnchors_Object, 
  idents = c(31, 37, 45, 64, 11, 61, 47, 32, 63), invert = FALSE)

#Endothelial Populations
HumanFetal_combined_IntegrationAnchors_Object_Endothelial <- subset(
  HumanFetal_combined_IntegrationAnchors_Object, 
  idents = c(46, 14, 16, 60, 52), invert = FALSE)

#Proliferating Populations
HumanFetal_combined_IntegrationAnchors_Object_Proliferating <- subset(
  HumanFetal_combined_IntegrationAnchors_Object, 
  idents = c(38, 41, 22), invert = FALSE)

#### Run CellFindR on each Broad Group object created above. NOTE: must run CellFindR functions script for in order to get the functions ####
#Mesenchymal cells
HumanFetal_combined_IntegrationAnchors_Object_Mesenchymal <- RunPCA(object = HumanFetal_combined_IntegrationAnchors_Object_Mesenchymal, 
                                                                    npcs = 30, 
                                                                    verbose = FALSE)
HumanFetal_combined_IntegrationAnchors_Object_Mesenchymal <- RunUMAP(object = HumanFetal_combined_IntegrationAnchors_Object_Mesenchymal, 
                                                                     reduction = "pca", dims = 1:18)
res <- find_res(HumanFetal_combined_IntegrationAnchors_Object_Mesenchymal, 
                initial_res = 0.1, jump = 0.1)
HumanFetal_combined_IntegrationAnchors_Object_Mesenchymal <- FindNeighbors(object = HumanFetal_combined_IntegrationAnchors_Object_Mesenchymal, 
                                                                           reduction = "pca", dims = 1:18)
HumanFetal_combined_IntegrationAnchors_Object_Mesenchymal <- FindClusters(HumanFetal_combined_IntegrationAnchors_Object_Mesenchymal, 
                                                                          resolution = res)
output_folder <- paste("/your_directory/All Integrated Mesenchymal CellFindR")
dir.create(output_folder)
HumanFetal_combined_IntegrationAnchors_Object_Mesenchymal_labeled <- sub_clustering(HumanFetal_combined_IntegrationAnchors_Object_Mesenchymal, 
                                                                                    output_folder)
rm(HumanFetal_combined_IntegrationAnchors_Object_Mesenchymal)

#Endocrine CellFindR
HumanFetal_combined_IntegrationAnchors_Object_Endocrine <- RunPCA(object = HumanFetal_combined_IntegrationAnchors_Object_Endocrine, 
                                                                  npcs = 30, 
                                                                  verbose = FALSE) 
HumanFetal_combined_IntegrationAnchors_Object_Endocrine <- RunUMAP(object = HumanFetal_combined_IntegrationAnchors_Object_Endocrine, 
                                                                   reduction = "pca", 
                                                                   dims = 1:10)
res <- find_res(HumanFetal_combined_IntegrationAnchors_Object_Endocrine,
                initial_res = 0.1, 
                jump = 0.1)
HumanFetal_combined_IntegrationAnchors_Object_Endocrine <- FindNeighbors(HumanFetal_combined_IntegrationAnchors_Object_Endocrine, 
                                                                         dims = 1:10)
HumanFetal_combined_IntegrationAnchors_Object_Endocrine <- FindClusters(HumanFetal_combined_IntegrationAnchors_Object_Endocrine, 
                                                                        resolution = res)
output_folder <- paste("/your_directory/All Integrated Endocrine CellFindR")
dir.create(output_folder)
HumanFetal_combined_IntegrationAnchors_Object_Endocrine_labeled <- sub_clustering(HumanFetal_combined_IntegrationAnchors_Object_Endocrine, 
                                                                                  output_folder)
rm(HumanFetal_combined_IntegrationAnchors_Object_Endocrine)

#Exocrine CellFindR
HumanFetal_combined_IntegrationAnchors_Object_Exocrine <- RunPCA(object = HumanFetal_combined_IntegrationAnchors_Object_Exocrine, 
                                                                 npcs = 30, 
                                                                 verbose = FALSE)
HumanFetal_combined_IntegrationAnchors_Object_Exocrine <- RunUMAP(object = HumanFetal_combined_IntegrationAnchors_Object_Exocrine, 
                                                                  reduction = "pca", dims = 1:8)
res <- find_res(HumanFetal_combined_IntegrationAnchors_Object_Exocrine, initial_res = 0.1, jump = 0.1)
HumanFetal_combined_IntegrationAnchors_Object_Exocrine <- FindNeighbors(HumanFetal_combined_IntegrationAnchors_Object_Exocrine, dims = 1:8)
HumanFetal_combined_IntegrationAnchors_Object_Exocrine <- FindClusters(HumanFetal_combined_IntegrationAnchors_Object_Exocrine, 
                                                                       resolution = res)
output_folder <- paste("/your_directory/All Integrated Exocrine CellFindR")
dir.create(output_folder)
HumanFetal_combined_IntegrationAnchors_Object_Exocrine_labeled <- sub_clustering(HumanFetal_combined_IntegrationAnchors_Object_Exocrine, output_folder)
rm(HumanFetal_combined_IntegrationAnchors_Object_Exocrine)

#Neuronal CellFindR
HumanFetal_combined_IntegrationAnchors_Object_Neuronal <- RunPCA(object = HumanFetal_combined_IntegrationAnchors_Object_Neuronal, 
                                                                 npcs = 30, verbose = FALSE)
HumanFetal_combined_IntegrationAnchors_Object_Neuronal <- RunUMAP(object = HumanFetal_combined_IntegrationAnchors_Object_Neuronal, 
                                                                  reduction = "pca", dims = 1:5)
res <- find_res(HumanFetal_combined_IntegrationAnchors_Object_Neuronal, initial_res = 0.1, jump = 0.1)
HumanFetal_combined_IntegrationAnchors_Object_Neuronal <- FindNeighbors(HumanFetal_combined_IntegrationAnchors_Object_Neuronal)
HumanFetal_combined_IntegrationAnchors_Object_Neuronal <- FindClusters(HumanFetal_combined_IntegrationAnchors_Object_Neuronal, 
                                                                       resolution = res)
output_folder <- paste("/your_directory/All Integrated Neuronal CellFindR")
dir.create(output_folder)
HumanFetal_combined_IntegrationAnchors_Object_Neuronal_labeled <- sub_clustering(HumanFetal_combined_IntegrationAnchors_Object_Neuronal, 
                                                                                 output_folder)
rm(HumanFetal_combined_IntegrationAnchors_Object_Neuronal)

#Endothelial CellFindR 
HumanFetal_combined_IntegrationAnchors_Object_Endothelial <- RunPCA(object = HumanFetal_combined_IntegrationAnchors_Object_Endothelial, 
                                                                    npcs = 30, verbose = FALSE)
HumanFetal_combined_IntegrationAnchors_Object_Endothelial <- RunUMAP(object = HumanFetal_combined_IntegrationAnchors_Object_Endothelial, 
                                                                     reduction = "pca", dims = 1:5)
res <- find_res(HumanFetal_combined_IntegrationAnchors_Object_Endothelial, initial_res = 0.1, jump = 0.1)
HumanFetal_combined_IntegrationAnchors_Object_Endothelial <- FindNeighbors(HumanFetal_combined_IntegrationAnchors_Object_Endothelial)
HumanFetal_combined_IntegrationAnchors_Object_Endothelial <- FindClusters(HumanFetal_combined_IntegrationAnchors_Object_Endothelial, 
                                                                          resolution = res)
output_folder <- paste("/your_directory/All Integrated Endothelial CellFindR")
dir.create(output_folder)
HumanFetal_combined_IntegrationAnchors_Object_Endothelial_labeled <- sub_clustering(HumanFetal_combined_IntegrationAnchors_Object_Endothelial, 
                                                                                    output_folder)
rm(HumanFetal_combined_IntegrationAnchors_Object_Endothelial)

#Proliferating CellFindR 
HumanFetal_combined_IntegrationAnchors_Object_Proliferating <- RunPCA(object = HumanFetal_combined_IntegrationAnchors_Object_Proliferating, 
                                                                      npcs = 30, verbose = FALSE)
HumanFetal_combined_IntegrationAnchors_Object_Proliferating <- RunUMAP(object = HumanFetal_combined_IntegrationAnchors_Object_Proliferating, 
                                                                       reduction = "pca", dims = 1:20)
res <- find_res(HumanFetal_combined_IntegrationAnchors_Object_Proliferating, initial_res = 0.1, jump = 0.1)
HumanFetal_combined_IntegrationAnchors_Object_Proliferating <- FindNeighbors(HumanFetal_combined_IntegrationAnchors_Object_Proliferating)
HumanFetal_combined_IntegrationAnchors_Object_Proliferating <- FindClusters(HumanFetal_combined_IntegrationAnchors_Object_Proliferating, 
                                                                            resolution = res)
DimPlot(HumanFetal_combined_IntegrationAnchors_Object_Proliferating, reduction = "umap",  label = TRUE)
output_folder <- paste("/your_directory/All Integrated Proliferating CellFindR")
dir.create(output_folder)
HumanFetal_combined_IntegrationAnchors_Object_Proliferating_labeled <- sub_clustering(HumanFetal_combined_IntegrationAnchors_Object_Proliferating, 
                                                                                      output_folder)
rm(HumanFetal_combined_IntegrationAnchors_Object_Proliferating)

#Immune CellFindR 
HumanFetal_combined_IntegrationAnchors_Object_Immune <- RunPCA(object = HumanFetal_combined_IntegrationAnchors_Object_Immune, 
                                                               npcs = 30, verbose = FALSE)
HumanFetal_combined_IntegrationAnchors_Object_Immune <- RunUMAP(object = HumanFetal_combined_IntegrationAnchors_Object_Immune, 
                                                                reduction = "pca", dims = 1:10)
res <- find_res(HumanFetal_combined_IntegrationAnchors_Object_Immune, initial_res = 0.1, jump = 0.1)
HumanFetal_combined_IntegrationAnchors_Object_Immune <- FindNeighbors(HumanFetal_combined_IntegrationAnchors_Object_Immune)
HumanFetal_combined_IntegrationAnchors_Object_Immune <- FindClusters(HumanFetal_combined_IntegrationAnchors_Object_Immune, resolution = res)
output_folder <- paste("/your_directory/All Integrated Immune CellFindR")
dir.create(output_folder)
HumanFetal_combined_IntegrationAnchors_Object_Immune_labeled <- sub_clustering(HumanFetal_combined_IntegrationAnchors_Object_Immune, 
                                                                               output_folder)
rm(HumanFetal_combined_IntegrationAnchors_Object_Immune)


#### Map back the CellFindR clustering onto the merged (non-subsetted) object ####
#Attach the Broad Group cell type before each cluster 
endocrine <- data.frame('cellname' = rownames(HumanFetal_combined_IntegrationAnchors_Object_Endocrine_labeled@meta.data), 
                        'CellfindR' = paste('Endocrine', HumanFetal_combined_IntegrationAnchors_Object_Endocrine_labeled@meta.data$CellfindR, sep = '_'))
endothelial <- data.frame('cellname' = rownames(HumanFetal_combined_IntegrationAnchors_Object_Endothelial_labeled@meta.data), 
                          'CellfindR' = paste('Endothelial', HumanFetal_combined_IntegrationAnchors_Object_Endothelial_labeled@meta.data$CellfindR, sep = '_'))
combined_labels <- rbind(endocrine, endothelial)

Exocrine <- data.frame('cellname' = rownames(HumanFetal_combined_IntegrationAnchors_Object_Exocrine_labeled@meta.data), 
                       'CellfindR' = paste('Exocrine', HumanFetal_combined_IntegrationAnchors_Object_Exocrine_labeled@meta.data$CellfindR, sep = '_'))
combined_labels <- rbind(combined_labels, Exocrine)

Immune <- data.frame('cellname' = rownames(HumanFetal_combined_IntegrationAnchors_Object_Immune_labeled@meta.data), 
                     'CellfindR' = paste('Immune', HumanFetal_combined_IntegrationAnchors_Object_Immune_labeled@meta.data$CellfindR, sep = '_'))
combined_labels <- rbind(combined_labels, Immune)

Mesenchymal <- data.frame('cellname' = rownames(HumanFetal_combined_IntegrationAnchors_Object_Mesenchymal_labeled@meta.data), 
                          'CellfindR' = paste('Mesenchymal', HumanFetal_combined_IntegrationAnchors_Object_Mesenchymal_labeled@meta.data$CellfindR, sep = '_'))
combined_labels <- rbind(combined_labels, Mesenchymal)

Neuronal <- data.frame('cellname' = rownames(HumanFetal_combined_IntegrationAnchors_Object_Neuronal_labeled@meta.data), 
                       'CellfindR' = paste('Neuronal', HumanFetal_combined_IntegrationAnchors_Object_Neuronal_labeled@meta.data$CellfindR, sep = '_'))
combined_labels <- rbind(combined_labels, Neuronal)

Proliferating <- data.frame('cellname' = rownames(HumanFetal_combined_IntegrationAnchors_Object_Proliferating_labeled@meta.data), 
                            'CellfindR' = paste('Proliferating', HumanFetal_combined_IntegrationAnchors_Object_Proliferating_labeled@meta.data$CellfindR, sep = '_'))
combined_labels <- rbind(combined_labels, Proliferating)

combined_labels <- combined_labels[order(combined_labels$cellname),]
write.csv(combined_labels, '/your_directory/Desktop/labels.csv')
combined_labels <- read.csv(file = '/your_directory/labels.csv', header = TRUE)

list_label = c()
count = 0
# loop through all of the entire cellnames:
for (cell in rownames(HumanFetal_combined_IntegrationAnchors_Object@meta.data)){
  # if in these groups
  if (cell %in% combined_labels$cellname){
    list_label= c(list_label, toString(combined_labels[combined_labels$cellname == cell,]$CellfindR))
  }
  else{
    list_label = c(list_label, NA)
  }
  count = count + 1
  print(count)
}

HumanFetal_combined_IntegrationAnchors_Object@meta.data$CellFindR_update = list_label
HumanFetal_combined_IntegrationAnchors_Object_CellFindR <- SetIdent(HumanFetal_combined_IntegrationAnchors_Object, value = 'CellFindR_update')
DimPlot(HumanFetal_combined_IntegrationAnchors_Object_CellFindR, reduction = "umap", label = FALSE, pt.size = 1.5, 
        label.size = 6) +NoLegend()
rm(HumanFetal_combined_IntegrationAnchors_Object)

#Re-order the cluster to to group Broad Groups together
levels(x = HumanFetal_combined_IntegrationAnchors_Object_CellFindR) <- c("Mesenchymal_0", "Mesenchymal_1.0", "Mesenchymal_1.1", "Mesenchymal_2", "Mesenchymal_3.0", 
                                                                         "Mesenchymal_3.1", "Mesenchymal_3.2", "Mesenchymal_3.3", "Mesenchymal_4.0.0", "Mesenchymal_4.0.1", 
                                                                         "Mesenchymal_4.1", "Mesenchymal_4.2", "Mesenchymal_5", "Mesenchymal_6.0", "Mesenchymal_6.1", "Mesenchymal_6.2.0",
                                                                         "Mesenchymal_6.2.1", "Mesenchymal_7",
                                                                         
                                                                         "Exocrine_0.0", "Exocrine_0.1", "Exocrine_1", "Exocrine_2.0", "Exocrine_2.1.0", "Exocrine_2.1.1", "Exocrine_3.0", 
                                                                         "Exocrine_3.1", "Exocrine_3.2.0", "Exocrine_3.2.1", "Exocrine_3.2.2", "Exocrine_3.3.0", "Exocrine_3.3.1", "Exocrine_3.3.2", 
                                                                         "Exocrine_4.0", "Exocrine_4.1",
                                                                         
                                                                         "Endocrine_0", "Endocrine_1","Endocrine_2", "Endocrine_3.0.0", "Endocrine_3.0.1", "Endocrine_3.1.0", "Endocrine_3.1.1",
                                                                         "Endocrine_3.2", "Endocrine_4", "Endocrine_5.0", "Endocrine_5.1", "Endocrine_6", 
                                                                         
                                                                         "Endothelial_0", "Endothelial_1", "Endothelial_2", "Endothelial_3.0", "Endothelial_3.1.0", "Endothelial_3.1.1", 
                                                                         "Endothelial_3.2.0", "Endothelial_3.2.1", "Endothelial_3.3.0", "Endothelial_3.3.1", "Endothelial_4.0", 
                                                                         "Endothelial_4.1.0", "Endothelial_4.1.1", "Endothelial_4.2", "Endothelial_5.0", "Endothelial_5.1",
                                                                         
                                                                         "Immune_0", "Immune_1", "Immune_2", "Immune_3", "Immune_4",
                                                                         "Immune_5", "Immune_6", "Immune_7", "Immune_8", "Immune_9", "Immune_10.0.0", "Immune_10.0.1", 
                                                                         "Immune_10.1", "Immune_11", "Immune_12.0", "Immune_12.1", "Immune_13", "Immune_14.0.0", 
                                                                         "Immune_14.0.1", "Immune_14.1.0", "Immune_14.1.1", "Immune_14.2", "Immune_15.0", "Immune_15.1", 
                                                                         "Immune_15.2", "Immune_15.3", "Immune_15.4", "Immune_16", "Immune_17.0", "Immune_17.1", 
                                                                         "Immune_17.2", "Immune_17.3", "Immune_18.0.0", "Immune_18.0.1", "Immune_18.1", "Immune_19", "Immune_20.0", 
                                                                         "Immune_20.1", "Immune_21.0", "Immune_21.1", "Immune_22.0", "Immune_22.1", "Immune_23", "Immune_24", "Immune_25", "Immune_26",
                                                                         
                                                                         "Proliferating_0.0", "Proliferating_0.1.0", 
                                                                         "Proliferating_0.1.1", "Proliferating_0.2.0", "Proliferating_0.2.1", "Proliferating_0.2.2",
                                                                         "Proliferating_0.3", "Proliferating_1", "Proliferating_2.0.0", "Proliferating_2.0.1", 
                                                                         "Proliferating_2.1.0", "Proliferating_2.1.1", "Proliferating_2.1.2", "Proliferating_2.2.0", 
                                                                         "Proliferating_2.2.1", "Proliferating_3", "Proliferating_4.0", "Proliferating_4.1",
                                                                         "Proliferating_5.0", "Proliferating_5.1",
                                                                         
                                                                         "Neuronal_0.0", "Neuronal_0.1", "Neuronal_0.2", "Neuronal_1.0.0", "Neuronal_1.0.1", "Neuronal_1.1.0",
                                                                         "Neuronal_1.1.1", "Neuronal_1.2", "Neuronal_2.0.0", "Neuronal_2.0.1", "Neuronal_2.1.0", "Neuronal_2.1.1", "Neuronal_2.2")


#### Remove endocrine clusters that we deem as non-biological (doublets, empty droplets) to get final dataset ####
HumanFetal_combined_IntegrationAnchors_Object_CellFindR_Final <- subset(HumanFetal_combined_IntegrationAnchors_Object_CellFindR, 
                                                                        idents = c("Endocrine_5.1", "Endocrine_6", "Endocrine_3.2", "Endocrine_2.0"), 
                                                                        invert = TRUE)
HumanFetal_combined_IntegrationAnchors_Object_CellFindR_Final$CellFindR_update <- plyr::mapvalues(
  x = HumanFetal_combined_IntegrationAnchors_Object_CellFindR_Final$CellFindR_update, 
  from = c("Mesenchymal_0", "Mesenchymal_1.0", "Mesenchymal_1.1", "Mesenchymal_2", "Mesenchymal_3.0", 
           "Mesenchymal_3.1", "Mesenchymal_3.2", "Mesenchymal_3.3", "Mesenchymal_4.0.0", "Mesenchymal_4.0.1", 
           "Mesenchymal_4.1", "Mesenchymal_4.2", "Mesenchymal_5", "Mesenchymal_6.0", "Mesenchymal_6.1", "Mesenchymal_6.2.0",
           "Mesenchymal_6.2.1", "Mesenchymal_7",
           
           "Exocrine_0.0", "Exocrine_0.1", "Exocrine_1", "Exocrine_2.0", "Exocrine_2.1.0", "Exocrine_2.1.1", "Exocrine_3.0", 
           "Exocrine_3.1", "Exocrine_3.2.0", "Exocrine_3.2.1", "Exocrine_3.2.2", "Exocrine_3.3.0", "Exocrine_3.3.1", "Exocrine_3.3.2", 
           "Exocrine_4.0", "Exocrine_4.1",
           
           "Endocrine_0", "Endocrine_1", "Endocrine_3.0.0", "Endocrine_3.0.1", "Endocrine_3.1.0", "Endocrine_3.1.1",
           "Endocrine_4", "Endocrine_5.0", 
           
           "Endothelial_0", "Endothelial_1", "Endothelial_2", "Endothelial_3.0", "Endothelial_3.1.0", "Endothelial_3.1.1", 
           "Endothelial_3.2.0", "Endothelial_3.2.1", "Endothelial_3.3.0", "Endothelial_3.3.1", "Endothelial_4.0", 
           "Endothelial_4.1.0", "Endothelial_4.1.1", "Endothelial_4.2", "Endothelial_5.0", "Endothelial_5.1",
           
           "Immune_0", "Immune_1", "Immune_2", "Immune_3", "Immune_4",
           "Immune_5", "Immune_6", "Immune_7", "Immune_8", "Immune_9", "Immune_10.0.0", "Immune_10.0.1", 
           "Immune_10.1", "Immune_11", "Immune_12.0", "Immune_12.1", "Immune_13", "Immune_14.0.0", 
           "Immune_14.0.1", "Immune_14.1.0", "Immune_14.1.1", "Immune_14.2", "Immune_15.0", "Immune_15.1", 
           "Immune_15.2", "Immune_15.3", "Immune_15.4", "Immune_16", "Immune_17.0", "Immune_17.1", 
           "Immune_17.2", "Immune_17.3", "Immune_18.0.0", "Immune_18.0.1", "Immune_18.1", "Immune_19", "Immune_20.0", 
           "Immune_20.1", "Immune_21.0", "Immune_21.1", "Immune_22.0", "Immune_22.1", "Immune_23", "Immune_24", "Immune_25", "Immune_26",
           
           "Proliferating_0.0", "Proliferating_0.1.0", 
           "Proliferating_0.1.1", "Proliferating_0.2.0", "Proliferating_0.2.1", "Proliferating_0.2.2",
           "Proliferating_0.3", "Proliferating_1", "Proliferating_2.0.0", "Proliferating_2.0.1", 
           "Proliferating_2.1.0", "Proliferating_2.1.1", "Proliferating_2.1.2", "Proliferating_2.2.0", 
           "Proliferating_2.2.1", "Proliferating_3", "Proliferating_4.0", "Proliferating_4.1",
           "Proliferating_5.0", "Proliferating_5.1",
           
           "Neuronal_0.0", "Neuronal_0.1", "Neuronal_0.2", "Neuronal_1.0.0", "Neuronal_1.0.1", "Neuronal_1.1.0",
           "Neuronal_1.1.1", "Neuronal_1.2", "Neuronal_2.0.0", "Neuronal_2.0.1", "Neuronal_2.1.0", "Neuronal_2.1.1", "Neuronal_2.2"), 
  
  to = c("Mesenchymal_0", "Mesenchymal_1.0", "Mesenchymal_1.1", "Mesenchymal_2", "Mesenchymal_3.0", 
         "Mesenchymal_3.1", "Mesenchymal_3.2", "Mesenchymal_3.3", "Mesenchymal_4.0.0", "Mesenchymal_4.0.1", 
         "Mesenchymal_4.1", "Mesenchymal_4.2", "Mesenchymal_5", "Mesenchymal_6.0", "Mesenchymal_6.1", "Mesenchymal_6.2.0",
         "Mesenchymal_6.2.1", "Mesenchymal_7",
         
         "Exocrine_0.0", "Exocrine_0.1", "Exocrine_1", "Exocrine_2.0", "Exocrine_2.1.0", "Exocrine_2.1.1", "Exocrine_3.0", 
         "Exocrine_3.1", "Exocrine_3.2.0", "Exocrine_3.2.1", "Exocrine_3.2.2", "Exocrine_3.3.0", "Exocrine_3.3.1", "Exocrine_3.3.2", 
         "Exocrine_4.0", "Exocrine_4.1",
         
         "Endocrine_0", "Endocrine_1","Endocrine_2", "Endocrine_3.0.0", "Endocrine_3.0.1", "Endocrine_3.1.0", "Endocrine_3.1.1",
         "Endocrine_4", "Endocrine_5",
         
         "Endothelial_0", "Endothelial_1", "Endothelial_2", "Endothelial_3.0", "Endothelial_3.1.0", "Endothelial_3.1.1", 
         "Endothelial_3.2.0", "Endothelial_3.2.1", "Endothelial_3.3.0", "Endothelial_3.3.1", "Endothelial_4.0", 
         "Endothelial_4.1.0", "Endothelial_4.1.1", "Endothelial_4.2", "Endothelial_5.0", "Endothelial_5.1",
         
         "Immune_0", "Immune_1", "Immune_2", "Immune_3", "Immune_4",
         "Immune_5", "Immune_6", "Immune_7", "Immune_8", "Immune_9", "Immune_10.0.0", "Immune_10.0.1", 
         "Immune_10.1", "Immune_11", "Immune_12.0", "Immune_12.1", "Immune_13", "Immune_14.0.0", 
         "Immune_14.0.1", "Immune_14.1.0", "Immune_14.1.1", "Immune_14.2", "Immune_15.0", "Immune_15.1", 
         "Immune_15.2", "Immune_15.3", "Immune_15.4", "Immune_16", "Immune_17.0", "Immune_17.1", 
         "Immune_17.2", "Immune_17.3", "Immune_18.0.0", "Immune_18.0.1", "Immune_18.1", "Immune_19", "Immune_20.0", 
         "Immune_20.1", "Immune_21.0", "Immune_21.1", "Immune_22.0", "Immune_22.1", "Immune_23", "Immune_24", "Immune_25", "Immune_26",
         
         "Proliferating_0.0", "Proliferating_0.1.0", 
         "Proliferating_0.1.1", "Proliferating_0.2.0", "Proliferating_0.2.1", "Proliferating_0.2.2",
         "Proliferating_0.3", "Proliferating_1", "Proliferating_2.0.0", "Proliferating_2.0.1", 
         "Proliferating_2.1.0", "Proliferating_2.1.1", "Proliferating_2.1.2", "Proliferating_2.2.0", 
         "Proliferating_2.2.1", "Proliferating_3", "Proliferating_4.0", "Proliferating_4.1",
         "Proliferating_5.0", "Proliferating_5.1",
         
         "Neuronal_0.0", "Neuronal_0.1", "Neuronal_0.2", "Neuronal_1.0.0", "Neuronal_1.0.1", "Neuronal_1.1.0",
         "Neuronal_1.1.1", "Neuronal_1.2", "Neuronal_2.0.0", "Neuronal_2.0.1", "Neuronal_2.1.0", "Neuronal_2.1.1", "Neuronal_2.2"))


