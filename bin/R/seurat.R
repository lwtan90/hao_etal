## seurat.R ##
# Author: Wilson Tan
# Date: 9/6/2022
# This script takes output from cellranger for downstream analysis:
# S1: Quality control per sample
# S2: Filter low-quality cells based on "aggregated QC" with other replictaes.
# S3: Normalization/Sctransform 
# S4: DimReduction, Clustering, and findNeighbours
# S5: Output Data for Integrative Analysis
#############################################################

## Loading libraries
## General (single-cell)
library(Seurat)
## SCTransform-specific
library(sctransform)
library(glmGamPoi)
# General packages
library(scales)
library(ggplot2)
library(patchwork)
require(data.table)
require(stringr)

## Setting seed for random number generation
set.seed(1234)

## For PNG file
## Stanford cluster doesnt have PNG driver
## Optional for other servers
options(bitmapType="cairo")


# Section 0: Reading input from cellranger
## This is the filtered compressed data produced by cellranger
## Seurat has other formats options as well
counts = Read10X("filtered_feature_bc_matrix")
## Reading in the count matrix to form seurat object
### min.cells = 10 and min.features = 300 are just default values (can be tuned again if cell numbers too little)
heart <- CreateSeuratObject(counts = counts, project = "FB", min.cells = 10, min.features = 300)


# Section 1: Quality control per sample
## Calculate percentage of reads contributed by mitochondrial gene, and stored in "percent.mt" in meta.data
## Depending on the genome, hg tends to be "MT-". Other genome might be "Mt-"
heart[["percent.mt"]] <- PercentageFeatureSet(heart, pattern = "mt-")
## pass the metadata as df, and write into text file for integrative QC process.
qc = heart@meta.data
write.table(qc,file="QC.txt",sep="\t",quote=FALSE)
## Checkpoint 1: save as RDS for loading after we have determined the threshold for removal of low quality cells.
saveRDS(heart,file="heart.preQC.RDS")


# Section 2: Filter low-quality cells based on "aggregated QC" with other replictaes.
heart = readRDS("heart.preQC.RDS")
## various thresholding was set up based on the threshold established with the other replicates.
## See qc.R
heart <- subset(heart, subset = nFeature_RNA > 900 & nFeature_RNA<15000 & nCount_RNA>1000 & nCount_RNA<10000 & percent.mt < 10)
## Checkpoint 2: save as RDS for the next step of analysis (assumed cleaned data, removed low quality cells)
saveRDS(heart, file="heart.postQC.RDS")


# Section 3: Normalization / SCTransform
## Here I normalize the data based on percentage mitochondrial RNA content. I have yet to test if the effect of mitochondrial is huge.
heart <- SCTransform(heart, method="glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
## Checkpoint 3: save as RDS for the next step of analysis.
saveRDS(heart, file = "heart.postQC.2.RDS")

# Section 4: DimReduction / Clustering / FindNeighbours
## Pretty Standard, not very important at this stage. More important at integration stage.
heart <- RunPCA(heart)
heart <- FindNeighbors(heart, dims=1:30)
heart <- FindClusters(heart, resolution = 0.5)
heart <- RunUMAP(heart,dims=1:30)
## Checkpoint 3: save as RDS, and ready for the integration analysis
saveRDS(heart, file = "heart.RDS")

