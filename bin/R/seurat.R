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
## For scATAC
###library(GenomeInfoDb)
###library(EnsDb.Hsapiens.v79)
## General packages
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


if(1){
heart = readRDS("heart.preQC.RDS")
##heart = readRDS("heart.preQC.RDS")
heart <- subset(heart, subset = nFeature_RNA > 900 & nFeature_RNA<15000 & nCount_RNA>1000 & nCount_RNA<10000 & percent.mt < 10)
saveRDS(heart, file="heart.postQC.RDS")

heart <- SCTransform(heart, method="glmGamPoi", vars.to.regress = "percent.mt", verbose = FALSE)
saveRDS(heart, file = "heart.postQC.2.RDS")
heart <- RunPCA(heart)

## diagnostic
##pc = as.data.frame(heart@reductions$pca@cell.embeddings)

##meta.data = heart@meta.data
##meta.data$nCount_RNA = log10(meta.data$nCount_RNA)
##meta.data$nFeature_RNA = log10(meta.data$nFeature_RNA)
##meta.data$nCount_SCT = log10(meta.data$nCount_SCT)
##meta.data$nFeature_SCT = log10(meta.data$nFeature_SCT)
## remove orig.ident
##meta.data = meta.data[,-1]

##require(reshape)
##require(ggplot2)
##require(scales)

##cor.pc.feature = cor(pc,meta.data)
##subset.cor.pc.feature = head(cor.pc.feature,10)
##cordata = melt(as.matrix(subset.cor.pc.feature))
##names(cordata) = c("PC","Feature","r")
##
##cordata = cordata[ cordata$PC!="PC_10", ]
##png("PC_correlation_heatmap.png",width=2500,height=2300,res=300)
##p1 <- ggplot(cordata,aes(x=PC,y=Feature)) + geom_tile(aes(fill=r)) + theme_bw() + theme(panel.grid=element_blank()) + geom_text(aes(x=PC,y=Feature,label=format(round(r,2),2))) + scale_fill_gradientn(values=rescale(c(-0.5,0.5)),colors=c("red","white","blue"),limits=c(-1,1))
##print(p1)
##dev.off()
##
##png("top_PC.png",width=2000,height=2000,res=300)
##VizDimLoadings(heart, dims = 1:2, reduction = "pca")
##dev.off()

png("PCA.png",width=5000,height=5000,res=300)
DimPlot(heart, reduction = "pca")
dev.off()

png("topgene_pca_heatmap.png",width=5000,height=5000,res=300)
DimHeatmap(heart, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()

png("elbow.png",width=2000,height=2000,res=300)
ElbowPlot(heart,ndims=50)
dev.off()

heart <- FindNeighbors(heart, dims=1:30)
heart <- FindClusters(heart, resolution = 0.5)
heart <- RunUMAP(heart,dims=1:30)

saveRDS(heart, file = "heart.rds")
umap.coord = as.data.frame(heart@reductions$umap@cell.embeddings)
umap.coord$group = as.factor(Idents(heart))
png("UMAP_SCTransform.png",width=3000,height=3000,res=300)
p1 <- ggplot(umap.coord,aes(x=UMAP_1,y=UMAP_2))+geom_point(aes(color=group),size=1.2)+theme_bw()+theme(panel.grid=element_blank())
##p1 <- p1 + scale_color_manual(values=c("darkred","red","orange","yellow","turquoise","blue","midnightblue","black","brown","purple","pink","salmon","bisque","grey30","grey60","green","darkgreen","yellowgreen"))
print(p1)
dev.off()
quit()
}



heart.markers <- FindAllMarkers(heart, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(heart.markers,file="heart.DE.markers.txt",sep="\t",quote=FALSE)

top10 <- heart.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)

saveRDS(heart.markers,file="heart.markers.RDS")

scale.data = GetAssayData(heart,slot = "scale.data")
target.scaledata = scale.data[ rownames(scale.data) %in% top10$gene, ]

require(reshape)
require(ggplot2)
require(scales)

metadata = heart@meta.data
metadata$cell =rownames(metadata)


scaled.data= as.matrix(target.scaledata)
scaled.data = (scaled.data-rowMeans(scaled.data))/apply(scaled.data,1,sd)
head(scaled.data[,1:20])

plotdata = melt(as.matrix(scaled.data))
names(plotdata)=c("gene","sample","z")
plotdata$cluster = metadata$seurat_clusters[match(plotdata$sample,metadata$cell)]
plotdata$de.group = top10$cluster[match(plotdata$gene,top10$gene)]


umap.agg = aggregate(plotdata$z,by=list(plotdata$sample,plotdata$de.group),FUN=mean)
names(umap.agg) = c("sample","de.group","mean")

#### UMAP
for(i in c(0:max(as.numeric(umap.agg$de.group))))
{
        test = umap.agg[ umap.agg$de.group==i, ]
        umap.coord$group = test$mean[match(rownames(umap.coord),test$sample)]
        umap.coord$group[ umap.coord$group>3]=3
        umap.coord$group[ umap.coord$group<(-3)]=-3
        print(head(umap.coord))
        png(paste("UMAP_DE_",i,".png",sep=""),width=3000,height=3000,res=300)
        p1 <- ggplot(umap.coord,aes(x=UMAP_1,y=UMAP_2))+geom_point(aes(color=group),size=1.2)+theme_bw()+theme(panel.grid=element_blank())
        p1 <- p1 + scale_color_gradientn(values=rescale(c(-2,0,2)),colors=c("black","red","orange","yellow"))
        p1 <- p1 + ggtitle(paste("UMAP for cluster marker",i))
        print(p1)
        dev.off()

}
####


plotdata$z[ plotdata$z>3]=3
plotdata$z[ plotdata$z<(-3)]=-3

options(bitmapType='cairo')
png("DE_heatmap.png",width=4000,height=4000,res=300)
p1<-ggplot(plotdata,aes(x=sample,y=gene))+geom_tile(aes(fill=z))+theme_bw()+theme(panel.grid=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.spacing=unit(0,"lines"))+facet_grid(de.group~cluster,space="free",scale="free")
p1<-p1+scale_fill_gradientn(values=rescale(c(-2,-1,-0.5,0,1,2)),colors=c("midnightblue","blue3","blue","turquoise","white","orange","red"))+xlab("")+ylab("")
print(p1)
dev.off()



#### Pathway
m_df<- msigdbr(species = "Rattus norvegicus", category = "C5")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

### pathway analysis
for(i in c(1:max(as.numeric(heart.markers$cluster))))
{
	cluster0.genes <- heart.markers %>% dplyr::filter(cluster == i) %>% arrange(desc(avg_log2FC)) %>% dplyr::select(gene,avg_log2FC)
	ranks <- cluster0.genes$avg_log2FC
	names(ranks)= cluster0.genes$gene
	fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000, scoreType="pos")
	fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
	print(head(fgseaResTidy))
	temp2 = as.data.frame(fgseaResTidy)
	temp = fgseaResTidy %>% filter(pval < 0.01) %>% head(n= 20)
	temp = as.data.frame(temp)
	print(temp)
	png(paste("FGSEA_",i,".png",sep=""),width=3000,height=3000,res=300)
	p1 <- ggplot(temp, aes(reorder(pathway, NES), NES)) 
  	p1 <- p1 + geom_col(aes(fill= NES < 7.5)) + coord_flip() + labs(x="Pathway", y="Normalized Enrichment Score", title=paste("GO NES from GSEA Cluster",i))
	p1 <- p1 + theme_minimal()
	print(p1)
	dev.off()
	temp2 = temp2[,c(1:5)]
	write.table(temp2,file=paste("Pathway_",i,".txt",sep=""),row.names=F,sep="\t",quote=FALSE)
}



if(0){
png("MARKER_featureplot.png",width=6000,height=6000,res=300)
FeaturePlot(heart, features=c("Myl2","Tnnt2","Myh6","Nppa","Dach1","Emcn","Egfl7","Vwf","Cdh5","Tie1","Postn","Dcn","Fn1","Rgs5","Abcc9","Pcdh7","Cd163","Mrc1","Ikzf1","Fyb1"),raster=FALSE)
dev.off()

options(bitmapType="cairo")
heart = readRDS("heart.rds")

heart <- RenameIdents(
        object = heart,
        "0"="CM",
        "1"="CM",
        "2"="CM",
        "3"="CM",
        "4"="FB",
        "5"="ENDO",
        "6"="dCM",
        "7"="FB",
        "8"="dCM",
        "9"="PERI",
        "10"="CM",
        "11"="FB",
        "12"="MACRO2",
        "13"="CM2",
        "14"="ENDO",
        "15"="ENDO",
        "16"="ENDO",
        "17"="FB"
)

umap.coord = as.data.frame(heart@reductions$umap@cell.embeddings)
umap.coord$group = as.factor(Idents(heart))
png("UMAP_namedTYPE.png",width=3000,height=3000,res=300)
p1 <- ggplot(umap.coord,aes(x=UMAP_1,y=UMAP_2))+geom_point(aes(color=group),size=1.2)+theme_bw()+theme(panel.grid=element_blank())
p1 <- p1 + scale_color_manual(values=c("darkred","red","orange","yellow","turquoise","blue","pink","black","brown","purple","pink","salmon","bisque","grey30","grey60","green","darkgreen","yellowgreen"))
print(p1)
dev.off()

saveRDS(heart,file = "majorHEART.rds")
}


heart.markers <- FindAllMarkers(heart, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.table(heart.markers,file="heart.DE.celltype.markers.txt",sep="\t",quote=FALSE)
saveRDS(heart.markers,file="heart.celltype.markers.RDS")


heart.markers = readRDS("heart.celltype.markers.RDS")
umap.coord = as.data.frame(heart@reductions$umap@cell.embeddings)
umap.coord$group = as.factor(Idents(heart))

top10 <- heart.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)


scale.data = GetAssayData(heart,slot = "scale.data")
target.scaledata = scale.data[ rownames(scale.data) %in% top10$gene, ]


metadata = heart@meta.data
metadata$cell =rownames(metadata)


scaled.data= as.matrix(target.scaledata)
scaled.data = (scaled.data-rowMeans(scaled.data))/apply(scaled.data,1,sd)
head(scaled.data[,1:20])

plotdata = melt(as.matrix(scaled.data))
names(plotdata)=c("gene","sample","z")
plotdata$cluster = Idents(heart)[match(plotdata$sample,names(Idents(heart)))]
plotdata$de.group = top10$cluster[match(plotdata$gene,top10$gene)]
head(plotdata)

umap.agg = aggregate(plotdata$z,by=list(plotdata$sample,plotdata$de.group),FUN=mean)
names(umap.agg) = c("sample","de.group","mean")

#### UMAP
ct = unique(Idents(heart))
for(i in c(1:length(ct)))
{
        test = umap.agg[ umap.agg$de.group==ct[i], ]
        umap.coord$group = test$mean[match(rownames(umap.coord),test$sample)]
        umap.coord$group[ umap.coord$group>3]=3
        umap.coord$group[ umap.coord$group<(-3)]=-3
        print(head(umap.coord))
        png(paste("UMAP_celltypeDE_",ct[i],".png",sep=""),width=3000,height=3000,res=300)
        p1 <- ggplot(umap.coord,aes(x=UMAP_1,y=UMAP_2))+geom_point(aes(color=group),size=1.2)+theme_bw()+theme(panel.grid=element_blank())
        p1 <- p1 + scale_color_gradientn(values=rescale(c(-2,0,2)),colors=c("black","red","orange","yellow"))
        p1 <- p1 + ggtitle(paste("UMAP for cluster marker",ct[i]))
        print(p1)
        dev.off()

}
####

plotdata$z[ plotdata$z>3]=3
plotdata$z[ plotdata$z<(-3)]=-3

options(bitmapType='cairo')
png("celltypeDE_heatmap.png",width=4000,height=4000,res=300)
p1<-ggplot(plotdata,aes(x=sample,y=gene))+geom_tile(aes(fill=z))+theme_bw()+theme(panel.grid=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),panel.spacing=unit(0,"lines"))+facet_grid(de.group~cluster,space="free",scale="free")
p1<-p1+scale_fill_gradientn(values=rescale(c(-2,-1,-0.5,0,1,2)),colors=c("midnightblue","blue3","blue","turquoise","white","orange","red"))+xlab("")+ylab("")
print(p1)
dev.off()

quit()

#### Pathway
m_df<- msigdbr(species = "Rattus norvegicus", category = "C5")
fgsea_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

### pathway analysis
for(i in c(1:length(ct)))
{
        cluster0.genes <- heart.markers %>% dplyr::filter(cluster == ct[i]) %>% arrange(desc(avg_log2FC)) %>% dplyr::select(gene,avg_log2FC)
        ranks <- cluster0.genes$avg_log2FC
        names(ranks)= cluster0.genes$gene
        fgseaRes<- fgsea(fgsea_sets, stats = ranks, nperm = 1000, scoreType="pos")
        fgseaResTidy <- fgseaRes %>% as_tibble() %>% arrange(desc(NES))
        print(head(fgseaResTidy))
        temp2 = as.data.frame(fgseaResTidy)
        temp = fgseaResTidy %>% filter(pval < 0.01) %>% head(n= 20)
        temp = as.data.frame(temp)
        print(temp)
        png(paste("FGSEA_",i,".png",sep=""),width=3000,height=3000,res=300)
        p1 <- ggplot(temp, aes(reorder(pathway, NES), NES))
        p1 <- p1 + geom_col(aes(fill= NES < 7.5)) + coord_flip() + labs(x="Pathway", y="Normalized Enrichment Score", title=paste("GO NES from GSEA Cluster",ct[i]))
        p1 <- p1 + theme_minimal()
        print(p1)
        dev.off()
        temp2 = temp2[,c(1:5)]
        write.table(temp2,file=paste("Pathway_",ct[i],".txt",sep=""),row.names=F,sep="\t",quote=FALSE)
}
}


