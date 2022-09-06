## dotplot.R ##
# Author: Wilson Tan
# Date: 9/6/2022
# Purpose: Generate dotplot and heatmap for the top 5/10 gene markers
################

## Library loading
require(Seurat)
require(ggplot2)

## PNG driver
options(bitmapType="cairo")

## global variables
features = c() # vector containing top 5 genes

## list directory
# Function that reads the marker differential file and return a list of top 5/10 genes
# Input: filename
# Output: a vector of gene names (top 5/10)
# Note: to change 5 to 10, just edit the return statement
top10marker <- function(filename)
{
	data = read.table(filename,header=TRUE)
	data = data[ order(-data$avg_log2FC),]
	return(head(rownames(data),5))
}

# The filenames (assume they started marker_)
# Loop through the list, and append the list to the features vector.
file.list = list.files(pattern="^marker_")
for( f in file.list)
{
	print(f)
	features = c(features, top10marker(f))
}

## Only allow unique gene elements
features = unique(features)
## Exclude mitochondrial
features = features[ grep("mt-",features,invert=TRUE) ]

## Output the marker file (for future uses: not sure)
write.table(features,file="markers.txt",sep="\t",quote=FALSE,row.names=F)

## extract top 5 markers per cluster
## plot the dotplot
## Read in Seurat object containing the integrated object
## Can be any. No restriction.
## If your data is .RData format, use load instead of readRDS
heart.combined = readRDS("SCT.prepped.RDS")

## In the current project, we have already renamed the cell types.
## This step is to arrange the columns of the dotplot
levels(heart.combined) <- c('CM','EC','Epicardial','FB','Lymphatic_Endo','Lymphocyte','MC','Neuron','Pericyte','Prolif_Macro','UNK1','UNK2')

## Use Seurat utility function to put the dotplot
p1 <- DotPlot(heart.combined, features = features) & coord_flip()
png("top5marker_dotplot.png",width=2500,height=5000,res=300)
p1
dev.off()
## Use Seurat utility function to plot the heatmap
png("top5marker_heatmap.png",width=6000,height=4000,res=300)
DoHeatmap(heart.combined, features = features, cells = 1:500, size = 4,angle = 90) + NoLegend()
dev.off()


