library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)

#raw data input and plot
pbmc.data <- Read10X(data.dir = './R/Single Cell Analysis/data/filtered_gene_bc_matrices/hg19')
pbmc <- CreateSeuratObject(counts=pbmc.data, project = 'pbmc3k', min.cells=3, min.features=200)
View(pbmc)

#add new columns to object metadata using [[]]
pbmc[['percent.mt']] <- PercentageFeatureSet(pbmc, pattern='^MT-')
VlnPlot(pbmc, features=c('nFeature_RNA','nCount_RNA','percent.mt'), ncol=3)

#Log transformation and variance stabilizing transformation
pbmc <- NormalizeData(pbmc, normalization.method = 'LogNormalize')
pbmc <- FindVariableFeatures(pbmc, selection.method = 'vst', nfeatures = 2000)

#plot variable features with and without labels
top10 <- head(VariableFeatures(pbmc),20)
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot=plot1, points=top10, repel=T)
plot1
plot2

#scaling
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features=all.genes)

#run pca & visualize it
pbmc <- RunPCA(pbmc, features=VariableFeatures(object=pbmc))
print(pbmc[['pca']], dims=1:2, nfeatures=5)
VizDimLoadings(pbmc, dims=1:2, nfeatures=15, reduction = 'pca')
DimPlot(pbmc, reduction = 'pca')
DimHeatmap(pbmc, dims=1, cells=500, balanced=T)

#or another alternative is using jackstraw
pbmc <- JackStraw(pbmc, num.replicate=100)
pbmc <- ScoreJackStraw(pbmc, dims=1:20)
JackStrawPlot(pbmc, dims=1:15)
#but can be substitute by ElbowPlot
ElbowPlot(pbmc)

