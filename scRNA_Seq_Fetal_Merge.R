
###########################################################################################################
### This script details the initial loading, filtering and merging of all human fetal pancreas scRNA Seq ## 
### datasets as shown in figure 1 of the paper.                                                          ##
###########################################################################################################

library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(cowplot)

#### Load individual data sets and merge 10x lanes #### 
#8w
Eight_w_lane1 <- Read10X(data.dir = directory_8w_1)
Eight_w_lane2 <- Read10X(data.dir = directory_8w_2)
Eight_w_lane1 <- CreateSeuratObject(counts = Eight_w_lane1, project = "Eight_w_lane1")
Eight_w_lane2 <- CreateSeuratObject(counts = Eight_w_lane2, project = "Eight_w_lane2")
Eight_w <- merge(Eight_w_lane1, y = Eight_w_lane2, add.cell.ids = c("lane1", "lane2"), project = "Eight_w")
rm(Eight_w_lane1, Eight_w_lane2)

#10w
Ten_w_lane1 <- Read10X(data.dir = directory_10w_1)
Ten_w_lane2 <- Read10X(data.dir = directory_10w_2)
Ten_w_lane1 <- CreateSeuratObject(counts = Ten_w_lane1, project = "Ten_w_lane1")
Ten_w_lane2 <- CreateSeuratObject(counts = Ten_w_lane2, project = "Ten_w_lane2")
Ten_w <- merge(Ten_w_lane1, y = Ten_w_lane2, add.cell.ids = c("lane1", "lane2"), project = "Ten_w")
rm(Ten_w_lane1, Ten_w_lane2)

#12w1
Twelve_w1_lane1 <- Read10X(data.dir = directory_Twelve_w1_lane1)
Twelve_w1_lane2 <- Read10X(data.dir = directory_Twelve_w1_lane2)
Twelve_w1_lane1 <- CreateSeuratObject(counts = Twelve_w1_lane1, project = "Twelve_w1_lane1")
Twelve_w1_lane2 <- CreateSeuratObject(counts = Twelve_w1_lane2, project = "Twelve_w1_lane2")
Twelve_w1 <- merge(Twelve_w1_lane1, y = Twelve_w1_lane2, add.cell.ids = c("lane1", "lane2"), project = "Twelve_w1")
rm(Twelve_w1_lane1, Twelve_w1_lane2)

#12w2
Twelve_w2_lane1 <- Read10X(data.dir = directory_Twelve_w2_lane1)
Twelve_w2_lane2 <- Read10X(data.dir = directory_Twelve_w2_lane2)
Twelve_w2_lane1 <- CreateSeuratObject(counts = Twelve_w2_lane1, project = "Twelve_w2_lane1")
Twelve_w2_lane2 <- CreateSeuratObject(counts = Twelve_w2_lane2, project = "Twelve_w2_lane2")
Twelve_w2 <- merge(Twelve_w2_lane1, y = Twelve_w2_lane2, add.cell.ids = c("lane1", "lane2"), project = "Twelve_w2")
rm(Twelve_w2_lane1, Twelve_w2_lane2)


#12w3_Epcam
Twelve_w3_Epcam_lane1 <- Read10X(data.dir = directory_Twelve_w3_Epcam_lane1)
Twelve_w3_Epcam_lane2 <- Read10X(data.dir = directory_Twelve_w3_Epcam_lane2)
Twelve_w3_Epcam_lane1 <- CreateSeuratObject(counts = Twelve_w3_Epcam_lane1, project = "Twelve_w3_Epcam_lane1")
Twelve_w3_Epcam_lane2 <- CreateSeuratObject(counts = Twelve_w3_Epcam_lane2, project = "Twelve_w3_Epcam_lane2")
Twelve_w3_Epcam <- merge(Twelve_w3_Epcam_lane1, y = Twelve_w3_Epcam_lane2, add.cell.ids = c("lane1", "lane2"), project = "Twelve_w3_Epcam")
rm(Twelve_w3_Epcam_lane1, Twelve_w3_Epcam_lane2)

#16w
Sixteen_w_lane1 <- Read10X(data.dir = directory_16w_lane1)
Sixteen_w_lane2 <- Read10X(data.dir = directory_16w_lane2)
Sixteen_w_lane1 <- CreateSeuratObject(counts = Sixteen_w_lane1, project = "Sixteen_w_lane1")
Sixteen_w_lane2 <- CreateSeuratObject(counts = Sixteen_w_lane2, project = "Sixteen_w_lane2")
Sixteen_w <- merge(Sixteen_w_lane1, y = Sixteen_w_lane2, add.cell.ids = c("lane1", "lane2"), project = "Sixteen_w")
rm(Sixteen_w_lane1, Sixteen_w_lane2)

#19w_1_CD45_neg
Ninteen_w_1_CD45_neg <- Read10X(data.dir = directory_19w_1_CD45_neg)
Ninteen_w_1_CD45_neg <- CreateSeuratObject(counts = Ninteen_w_1_CD45_neg, project = "Ninteen_w_1_CD45_neg")

#19w_1_Epcam
Ninteen_w_1_Epcam <- Read10X(data.dir = directory_19w_1_Epcam)
Ninteen_w_1_Epcam <- CreateSeuratObject(counts = Ninteen_w_1_Epcam, project = "Ninteen_w_1_Epcam")

#20w_Epcam
Twenty_w_Epcam <- Read10X(data.dir = directory_20w_Epcam)
Twenty_w_Epcam <- CreateSeuratObject(counts = Twenty_w_Epcam, project = "Twenty_w_Epcam")

#### QC filter all datasets ####
Eight_w <- subset(Eight_w, subset = nFeature_RNA > 200 & nFeature_RNA < 5500 & percent.mt < 10)
Ten_w <- subset(Ten_w, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
Twelve_w1 <- subset(Twelve_w1, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
Twelve_w2 <- subset(Twelve_w2, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 15)
Twelve_w3_Epcam <- subset(Twelve_w3_Epcam, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 & percent.mt < 15)
Sixteen_w <- subset(Sixteen_w, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
Ninteen_w_1_CD45_neg <- subset(Ninteen_w_1_CD45_neg, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 15)
Ninteen_w_1_Epcam <- subset(Ninteen_w_1_Epcam, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 15)
Twenty_w_Epcam <- subset(Twenty_w_Epcam, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 15)

#### Normalization and HVG ####
Eight_w <- NormalizeData(Eight_w)
Ten_w <- NormalizeData(Ten_w)
Twelve_w1 <- NormalizeData(Twelve_w1)
Twelve_w2 <- NormalizeData(Twelve_w2)
Twelve_w3_Epcam <- NormalizeData(Twelve_w3_Epcam)
Sixteen_w <- NormalizeData(Sixteen_w)
Ninteen_w_1_CD45_neg <- NormalizeData(Ninteen_w_1_CD45_neg)
Ninteen_w_1_Epcam <- NormalizeData(Ninteen_w_1_Epcam)
Twenty_w_Epcam <- NormalizeData(Twenty_w_Epcam)

Eight_w <- FindVariableFeatures(Eight_w, selection.method = "vst", nfeatures = 2000)
Ten_w <- FindVariableFeatures(Ten_w, selection.method = "vst", nfeatures = 2000)
Twelve_w1 <- FindVariableFeatures(Twelve_w1, selection.method = "vst", nfeatures = 2000)
Twelve_w2 <- FindVariableFeatures(Twelve_w2, selection.method = "vst", nfeatures = 2000)
Twelve_w3_Epcam <- FindVariableFeatures(Twelve_w3_Epcam, selection.method = "vst", nfeatures = 2000)
Sixteen_w <- FindVariableFeatures(Sixteen_w, selection.method = "vst", nfeatures = 2000)
Ninteen_w_1_Epcam <- FindVariableFeatures(Ninteen_w_1_Epcam, selection.method = "vst", nfeatures = 2000)
Twenty_w_Epcam <- FindVariableFeatures(Twenty_w_Epcam, selection.method = "vst", nfeatures = 2000)

#### Integrate data with Standard Integration #### 
HumanFetal_combined_IntegrationAnchors <- FindIntegrationAnchors(object.list = list(
  Eight_w, 
  Ten_w, 
  Twelve_w1, 
  seur_ob_10wpc_1_merged_QCfiltered, 
  Twelve_w2,
  Twelve_w3_Epcam,
  Sixteen_w,
  Ninteen_w_1_Epcam,
  Twenty_w_Epcam), 
  dims = 1:30)
HumanFetal_combined_IntegrationAnchors_Object <- IntegrateData(
  anchorset = HumanFetal_combined_IntegrationAnchors, 
  dims = 1:30)

#### Return to standard processing ####
HumanFetal_combined_IntegrationAnchors_Object <- ScaleData(object = HumanFetal_combined_IntegrationAnchors_Object)
HumanFetal_combined_IntegrationAnchors_Object <- RunPCA(object = HumanFetal_combined_IntegrationAnchors_Object, 
                                                        npcs = 30, 
                                                        verbose = FALSE)
HumanFetal_combined_IntegrationAnchors_Object <- FindNeighbors(object = HumanFetal_combined_IntegrationAnchors_Object, 
                                                               reduction = "pca", 
                                                               dims = 1:24)
HumanFetal_combined_IntegrationAnchors_Object <- FindClusters(HumanFetal_combined_IntegrationAnchors_Object, 
                                                              resolution = 2.0)
HumanFetal_combined_IntegrationAnchors_Object <- RunUMAP(HumanFetal_combined_IntegrationAnchors_Object)
