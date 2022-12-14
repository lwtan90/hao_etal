# Code Repo  
Author: Wilson Tan  
Date: 8/31/2022  
Description: Code repository for single-cell analysis  
  
  
## Required Packages:  
```
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
library(data.table)
library(stringr)

```  
  
## Suggested Workflow for Running the Scripts  
This is the flow that I suggest:  
```
# for individual replicate, run these:
## Untill you have the QC then stop
Rscript seurat.R
## Unless you have all replicates
## edit this script to the path of each replicates
Rscript qc.R

# Once you have all sctransformed data
## This is the step for data integration
Rscript integration_scRNAseq.R

# Perform Gene Marker Identitifcation
## Eg if you are interested in cluster "CM"
## Note: if running on RStudio, you can change celltype = args[1] into celltype = "CM"
Rscript findMARKER.R CM

# Utility Function
## Plot dotplot and heatmap for top 5 markers
Rscript dotplot.R  


```  
  


