## featurePLOT.R ##
#
# Author: Wilson Tan
# Date: 9/6/2022
# Purpose: This script plots featureplots based on gene of interest
#
###################

## Library Loading
library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr)

## PNG driver
options(bitmapType="cairo")

## User input for gene of interest
args = commandArgs(TRUE)
genetarget = args[1]

## read seurat object
seuOBJ = readRDS("SCT.prepped.RDS")
DefaultAssay(seuOBJ) = "SCT"
## added a column named "sample" to represent condition
## here we also rank the condition for visualization purpose
seuOBJ$sample = as.factor(seuOBJ$sample)
seuOBJ$sample = factor(seuOBJ$sample,levels=c("Sham","TAC","ART"))

## Function body to generate multi-plot panel for feature plots
FeaturePlotSingle<- function(obj, feature, metadata_column, ...){
	all_cells<- colnames(obj)
	groups<- levels(obj@meta.data[, metadata_column]) 
	minimal<- min(obj[['SCT']]@data[feature, ])
	maximal<- max(obj[['SCT']]@data[feature, ])
	ps<- list()
	for (group in groups) {
	  subset_indx<- obj@meta.data[, metadata_column] == group
	  subset_cells<- all_cells[subset_indx]
	  p<- FeaturePlot(obj, features = feature, cells= subset_cells, ...) +
	    scale_color_viridis_c(limits=c(minimal, maximal), direction = 1) +
	    ggtitle(paste(feature,"(",group,")")) +
	    theme(plot.title = element_text(size = 10, face = "bold"))
	  ps[[group]]<- p
	}
	
	
	return(ps)
}

layout = "ABC"

p_list<- FeaturePlotSingle(seuOBJ, feature=genetarget, metadata_column = "sample", pt.size = 0.5, order =TRUE)
png(paste(args[1],"_splitsample_featureplot.png",sep=""),width=2300,height=1000,res=300)
wrap_plots(p_list, guides='collect',design = layout)
dev.off()
