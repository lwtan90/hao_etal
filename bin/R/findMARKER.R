## findMARKER.R ##
#
# Author: Wilson Tan
# Date: 9/6/2022
# Purpose: perform marker discovery for each cluster / cell type, depending on the Ident of the seurat object
# Note: Set DefaultAssay to "SCT" for RNA
# Note: Assuming tha the data has been SCTrasnformed and Prepped (scaled using PrepSCTFindMarkers)
###################

### Loading packages
require(Seurat)

### Reading input
seuOBJ = readRDS("SCT.prepped.RDS")

### This is the user argument for which cluster is compared here
args = commandArgs(TRUE)
## a quick check to see if the Ident of the cells have been set correctly
table(Idents(seuOBJ))

## Performing differential gene expression analysis between the Ident of interest with the rest. (One vs All comparison)
## Default test is Wilcox. Only showing only.pos (upregulated in the cluster of interest)
## To save time, only focus on those gene expressed at leats 20% of the cluster of interest with logfc of 0.2 (edit table)
seuOBJ.markers <- FindMarkers(seuOBJ, assay="SCT", ident.1 = args[1], min.pct = 0.2, logfc.threshold = 0.2, only.pos=TRUE)

## output the de output as a text file
write.table(seuOBJ.markers,file=paste("marker_",args[1],".txt",sep=""),sep="\t",quote=FALSE)

