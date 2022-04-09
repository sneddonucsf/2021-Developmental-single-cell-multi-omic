# 2021-Developmental-single-cell-multi-omic

Scripts used for analysis of data in manuscript @ https://www.biorxiv.org/content/10.1101/2022.02.17.480942v1. 

SCRIPTS DESCRIPTION: 

scRNA_Seq_Fetal_Merge.R - Initialize Seurat objects from 10x Genomics CellRanger outs and computationally merge all datasets together using Seurat v3.2 Standard Integration. Script details loading the datasets, filtering, normalization, variable gene identification, integration, PCA analysis, clustering and UMAP visualization. Used for all human fetal scRNA Seq datasets. 

Broad_Group_CellFindR_Subcluster.R - Subset merged dataset for each Broad Group as defined by marker expression and run CellFindR sub-clustering. Script details sub-setting each Broad Group and running CellFindR, then mapping the data back onto the merged dataset.

CellFindR_Functions.R - Script for creating the CellFindR functions ran on the various datasets in the paper. 

Human_Fetal_Endocrine_scRNA_Seq_Slingshot_TradeSeq.R - Perform pseudotime reconstruction on human fetal endocrine dataset and differential gene expression analysis on a per-lineage basis. Script details performing pseudotemporal ordering with Slingshot and differential gene expression analysis with TradeSeq.

Human_Fetal_Endocrine_SCENIC_analysis.R - Perform Gene Regulatory Network (GRN) analysis on human fetal endocrine dataset using R package SCENIC. Script details building TF/RNA expression modules, pruning off target targets of TFs and area-under-the-curve (AUC) calculations. 

snATAC_Seq_Analysis.R - Initialize snATAC Seq analysis from 10x Genomics CellRanger outs with ArchR. Script details creating Arrow Files, initializing ArchR project, filtering, dimensional reduction, clustering, UMAP visualization, and snATAC-Seq scRNA Seq unconstrained integration. 

GWAS_SNP_analysis - Folder containing the relevant scripts and files in order to perform GWAS enrichment scoring on the snATAC Seq endocrine dataset. The folder contains all of the necessary files and scripts to calculate T1D and T2D SNP enrichment on a cell-type specific manner.

In_Vitro_FEV_KO_DESeq2.R - Generate differentially expressed genes between FEV WT and KO conditions from bulk RNA Sequencing datasets. Script details loading in count files from 6 samples (3 pairs of WT and KO differentiations) and performing differential gene expression analysis with DESeq2. 





