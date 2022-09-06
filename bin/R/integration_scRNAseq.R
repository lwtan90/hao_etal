## integration_scRNAseq.R ##
#
# Author: Wilson Tan
# Date: 9/6/2022
#
# Purpose: Main body of script for integration of scRNA-seq data
# 
# Section 1: Input reading from individual replicates
# Section 2: Data Integration
# Section 3: DimRed and FindNeightbous
# Section 4: Output RDS for downstream analysis
#############################

## Loaind required packages
library(Seurat)
library(ggplot2)
library(patchwork)
require(data.table)
require(sctransform)
library(glmGamPoi)

## For PNG file
options(bitmapType="cairo")

## Setting random number
set.seed(1234)

# Section 1: Input reading from indiidual replicates
## The seurat object should be SCTransformed
C1 = readRDS("../S1/S1/seurat/heart.rds") # Sham
C2 = readRDS("../T1/T1/seurat/heart.rds") # TAC
C3 = readRDS("../A1/A1/seurat/heart.rds") # ART

heart.list <- list(C1=C1,C2=C2,C3=C3)

## Feature selection
### Here I pick 3000 as default (recommended)
features <- SelectIntegrationFeatures(object.list=heart.list,nfeatures=3000)

## Preparing the SCTrasformed data for integeration
heart.list = PrepSCTIntegration(object.list=heart.list, anchor.features=features)

## Find a set of anchors for data integration
### You must set normalization.method = "SCT" to enable SCTrasnform
heart.anchors <- FindIntegrationAnchors(object.list=heart.list, normalization.method="SCT", anchor.features=features)

## Key steps in Data Integration
heart.combined <- IntegrateData(anchorset=heart.anchors,normalization.method="SCT")

# Section 3: DimReduction
## Key step: resolution = 0.5. I find resolution = 1 results in infinitely lots of clusters.
heart.combined <- RunPCA(heart.combined)
heart.combined <- RunUMAP(heart.combined,reduction="pca",dims=1:30)
heart.combined <- FindNeighbors(heart.combined,reduction="pca",dims=1:30)
heart.combined <- FindClusters(heart.combined,resolution=0.5)
save(heart.combined,file="step3.merged.RData")

## Optional: Depending on your sample, you can rename the sample
heart.combined$sample = rep("Sham")
heart.combined$sample[grep("_2",rownames(heart.combined@meta.data))] = "TAC"
heart.combined$sample[grep("_3",rownames(heart.combined@meta.data))] = "ART"

# Section 4: Output RDS
### This should be the final RDS
saveRDS(heart.combined,file="final.integrated.RDS")

## Output UMAP
png("UMAP_RNA_integration.png",width=3000,height=1500,res=300)
## Map colored by condition(Sham,TAC,ART)
p1 <- DimPlot(heart.combined, reduction = "umap", group.by = "sample")
## Map colored by clusters (resolution = 0.5)
p2 <- DimPlot(heart.combined, reduction = "umap", label = TRUE, repel = TRUE,raster = FALSE)
p1 + p2
dev.off()
