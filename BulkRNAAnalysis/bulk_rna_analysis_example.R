library(DESeq2)
library(vsn)
library(EnhancedVolcano)
library(pheatmap)
library(tidyverse)

ge.matrix <- read.table('bulk_RNAseq_raw_counts.txt.gz',
                        header=T, sep='\t')
dim(ge.matrix)
ge.matrix[1:4,1:4]

pheno.matrix <- read.table('bulk_RNAseq_metadata.txt.gz',
                           header=T, sep='\t', stringsAsFactors=T)
pheno.matrix[1:4,]

#make the first column to become sample id
rownames(pheno.matrix) <- pheno.matrix$sample_id

#Confirming whether the ge and pheno data has the same values
all(rownames(pheno.matrix)==colnames(ge.matrix))

#select data that we want to focus
stimTime <- '5d'
conditions <- c('Th2', 'Th0')
celltype <- 'CD4_Memory'

toSelect <- pheno.matrix$stimulation_time == stimTime &
  pheno.matrix$cytokine_condition %in% conditions &
  pheno.matrix$cell_type == celltype

pheno.matrix.subset <- pheno.matrix[toSelect,]
ge.matrix.subset <- ge.matrix[,toSelect]

#Converting data into DESeq that contains numerical data and the meta data
dds <- DESeqDataSetFromMatrix(countData=ge.matrix.subset,
                              colData= ,pheno.matrix.subset,
                              design= ~cytokine_condition)

#only data that has 10 or more infos are processed
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds


#apply a pseudocount of 1 and apply log2
normtfd <- normTransform(dds)
#Compare mean to sd
meanSdPlot(assay(normtfd))

#calculate rlog values to remove variance on the mean(variance stabilizing transformations VST) and regularized logarithm (rlog)
rltfd <- rlog(dds, blind=F)
meanSdPlot(assay(rltfd))

#Normalization
dds <- estimateSizeFactors(dds)
sizeFactors(dds)

plot(sizeFactors(dds),
     colSums(counts(dds, normalized=F)),
     xlab='Size Factor',
     ylab='Total Number of Reads',
     pch=19)

rltfd.pca <-prcomp(t(assay(rltfd)), scale=T)

require(factoextra)
fviz_eig(rltfd.pca, addlabels=T)

fviz_pca_ind(rltfd.pca)

#plot with conditions applied (using normal data not pca data)
plotPCA(rltfd, intgroup='sequencing_batch', ntop=26656) #sequencing batch isnot biased between 2 groups -> well done
#scale is not applied in plotPCA thus PC1 and PC2 variance number is different in plotpCA and prcomp

plotPCA(rltfd, intgroup='cytokine_condition')

dds <- DESeq(dds)
res <- results(dds)
dim(res)
res

#using padj because our test data is too large that even 5% of error will become many 
sum(res$padj <= 0.01 &
      abs(res$log2FoldChange)>1, na.rm=T)

#Visualization using Volcano Plot
EnhancedVolcano(res, lab=rownames(res),
                 x ='log2FoldChange', y='padj',
                 subtitle = 'Th2 vs Th0', labSize=3,
                 pCutoff = 0.01,
                 FCcutoff = 1,
                 drawConnectors = T)

#Visualization using heatmap (first we select our data, apply it to our results and make it into dataframe)
DEG.idx <- which(res$padj <= 0.01 &
                   abs(res$log2FoldChange) > 1)
res[DEG.idx,]
df <- as.data.frame(colData(dds)[,c('cytokine_condition','donor_id','sequencing_batch')])

pheatmap(assay(rltfd)[DEG.idx,], annotation_col=df,
         treeheight_row=0, treeheight_col=0, scale='row')
