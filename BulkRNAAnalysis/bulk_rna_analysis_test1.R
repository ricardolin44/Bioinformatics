library(DESeq2)
library(EnhancedVolcano)
library(vsn)

ge.matrix <- read.table('bulk_RNAseq_raw_counts.txt.gz',
                        header=T, sep='\t')
pheno.matrix <- read.table('bulk_RNAseq_metadata.txt.gz',
                           header=T, sep='\t', stringsAsFactors=T)

rownames(pheno.matrix) <- pheno.matrix$sample_id
all(rownames(pheno.matrix) == colnames(ge.matrix))

stimTime2 <- c('16h', '5d')
conditions2 <- 'Th0'
celltype2 <- 'CD4_Memory'

toSelect2 <- pheno.matrix$stimulation_time == stimTime2 &
  pheno.matrix$cytokine_condition == conditions2 &
  pheno.matrix$cell_type == celltype2

pheno.matrix.subset2 <- pheno.matrix[toSelect2,]
ge.matrix.subset2 <- ge.matrix[,toSelect2]

dds2 <- DESeqDataSetFromMatrix(countData = ge.matrix.subset2,
                               colData = pheno.matrix.subset2,
                               design = ~ stimulation_time)

dds2 <- DESeq(dds2)
res2 <- results(dds2)
res2

rltfd2 <- rlog(dds2, blind=F)
meanSDPlot(assay(rltfd2))



#plotting enhanced volcano
EnhancedVolcano(res2, lab=rownames(res2),
                x ='log2FoldChange', y='padj',
                subtitle = '16h vs 5d', labSize=3,
                pCutoff = 0.01,
                FCcutoff = 1,
                drawConnectors = T)

#plotting heatmap <- need to filter the data
filter.idx <- which(res2$padj<0.01 &
                      res2$log2FoldChange>1)
df2 <- as.data.frame(colData(dds2)[,c('stimulation_time', 'donor_id', 'sequencing_batch')]) #for annotating in graph
                     
pheatmap(assay(rltfd2)[filter.idx,], annotation_col = df2,
         treeheight_row=0, treeheight_col=0, scale='row')

