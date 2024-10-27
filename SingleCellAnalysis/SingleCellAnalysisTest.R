library(dplyr)
library(Seurat)
library(ggplot2)
library(patchwork)
library(celldex)
library(SingleR)

hpsc.data <- Read10X(data.dir='./Single Cell Analysis/data/GSM2572205')
hpsc <- CreateSeuratObject(counts = pbmc.data, project='hPSC', min.cells=3, min.features=200)
hpsc

hpsc[['percent.mt']] <- PercentageFeatureSet(hpsc, pattern='^MT-')
VlnPlot(hpsc, features=c('nFeature_RNA', 'nCount_RNA', 'percent.mt'), ncol=3)
plot1 <- FeatureScatter(hpsc, feature1='nCount_RNA', feature2='percent.mt')
plot2 <- FeatureScatter(hpsc, feature1='nCount_RNA', feature2='nFeature_RNA')
plot1 + plot2

#Normalize & Find Variable Features + Plotting
hpsc <- NormalizeData(hpsc, normalization.method = 'LogNormalize')
hpsc <- FindVariableFeatures(hpsc, selection.method = 'vst', nfeatures=2000)
top10 <- head(VariableFeatures(hpsc),10)
plot1 <- VariableFeaturePlot(hpsc)
plot2 <- LabelPoints(plot=plot1, points=top10, repel=T)
plot2

#Scaling data to makes all the data contribute evenly and do PCA before UMAP to get a better result
all.genes <- rownames(hpsc)
hpsc <- ScaleData(hpsc, features = all.genes)
hpsc <- RunPCA(hpsc, features = VariableFeatures(object=hpsc))
print(hpsc[['pca']], dims=1:2, nfeatures=2)

VizDimLoadings(hpsc, dims=1:2, nfeatures=15, reduction='pca')

DimPlot(hpsc, reduction='pca')
DimHeatmap(hpsc, dims = 1, cells = 500, balanced = T)

#Find Neighbors Clusters and Run UMAP
ElbowPlot(hpsc)
hpsc <- FindNeighbors(hpsc, dims=1:10)
hpsc <- FindClusters(hpsc, resolution = 0.5)

head(Idents(hpsc),5)
hpsc <- RunUMAP(hpsc, dims=1:10)
DimPlot(hpsc, reduction='umap')
DimHeatmap(hpsc, dims = 1, cells = 500, balanced = T)

#Look into the cluster1 data
cluster1.markers <- FindMarkers(hpsc, ident.1 = 1, min.pct=0.25)
head(cluster1.markers, n=5)
VlnPlot(hpsc, features = c(row.names(cluster1.markers)[1], row.names(cluster1.markers)[2]))

#Look into all markers distinguishing cluster 5 from cluster 0 and 3
cluster5.markers <- FindMarkers(hpsc, ident.1 = 5, ident.2 = c(0,3), min.pct=0.25)
head(cluster1.markers, n=5)
VlnPlot(hpsc, features = c(row.names(cluster5.markers)[1], row.names(cluster5.markers)[2]))

#Look into all markers
hpsc.markers <- FindAllMarkers(hpsc, only.pos=T, min.pct=0.25, logfc.threshold = 0.25)
x <- hpsc.markers %>% group_by(cluster) %>% top_n(n=2, wt=avg_log2FC)   #get top 2 markers for each cluster
FeaturePlot(hpsc, features = x$gene[1:8])
FeaturePlot(hpsc, features = x$gene[5:8])