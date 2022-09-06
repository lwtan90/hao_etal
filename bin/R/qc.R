## qc.R ##
#
# Purpose: perform the quality control across all replicates 
# Optional if you have only n=1
#
###########

# load packages
## packages
require(ggplot2)

## required for PNG driver
## Optional for other servers
options(bitmapType="cairo")

## Read the QC reports
meta1 = read.table("../A1/A1/seurat/QC.txt",header=TRUE)
meta1$sample = rep("A1")

meta2 = read.table("../T1/T1/seurat/QC.txt",header=TRUE)
meta2$sample = rep("T1")

meta3 = read.table("../S1/S1/seurat/QC.txt",header=TRUE)
meta3$sample = rep("S1")

meta = rbind(meta1,meta2,meta3)
summary(meta)

## plot nCount_RNA
p1 <- ggplot(meta,aes(x=nCount_RNA,fill=sample)) + geom_density(alpha=0.6) + theme_bw() + theme(panel.grid=element_blank())
p1 <- p1 + scale_fill_manual(values=c("red","orange","yellow","green","blue","turquoise")) + scale_x_log10() + geom_vline(xintercept=c(1000,100000))
png("nCount_RNA.png",width=1500,height=1300,res=300)
print(p1)
dev.off()

## plot nFeature_RNA
p1 <- ggplot(meta,aes(x=nFeature_RNA,fill=sample)) + geom_density(alpha=0.6) + theme_bw() + theme(panel.grid=element_blank())
p1 <- p1 + scale_fill_manual(values=c("red","orange","yellow","green","blue","turquoise")) + scale_x_log10() + geom_vline(xintercept=c(1000,15000))
png("nFeature_RNA.png",width=1500,height=1300,res=300)
print(p1)
dev.off() 

## plot percent.mt
p1 <- ggplot(meta,aes(x=percent.mt,fill=sample)) + geom_density(alpha=0.6) + theme_bw() + theme(panel.grid=element_blank())
p1 <- p1 + scale_fill_manual(values=c("red","orange","yellow","green","blue","turquoise"))  + geom_vline(xintercept=40)
png("percent.mt.png",width=1500,height=1300,res=300)
print(p1)
dev.off()

## plot correlation between nCount_RNA and nFeature_RNA
p1 <- ggplot(meta,aes(x=nCount_RNA,y=nFeature_RNA)) + geom_bin_2d(bins=1000) +scale_fill_continuous(type = "viridis")+ theme_bw() + theme(panel.grid=element_blank())+scale_x_log10()
p1 <- p1 + facet_grid(.~sample)
png("nCount_RNA_nFeature_RNA.png",width=2800,height=1300,res=300)
print(p1)
dev.off()


