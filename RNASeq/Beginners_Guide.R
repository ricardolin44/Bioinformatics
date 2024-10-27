# library(renv)
# renv::init()
library(edgeR)
library(readr)
library(DESeq2)
library(Seurat)
library(ggplot2)
library(pheatmap)
library(patchwork)

system.time({
    rna_seq_data <- read_tsv("./RNA_Seq/data/GSE116583_transplant.am.htseq.all.rpkm.txt.gz", col_names = TRUE)
    Sys.sleep(2) # Example command that pauses for 2 seconds
})

View(rna_seq_data)
rna_seq_data <- as.data.frame(rna_seq_data) # need to convert it into data frame first
rownames(rna_seq_data) <- rna_seq_data$Symbol
rna_seq_data <- rna_seq_data[, -1]
rna_seq_data <- as.data.frame(lapply(rna_seq_data, as.integer))

col_data <- data.frame(
    row.names = colnames(rna_seq_data),
    condition = factor(c(rep("Naive", 4), rep("Allo_2H", 4), rep("Allo_24H", 4)))
)
col_data$Sample <- rownames(col_data)
col_data$num <- seq(1,12)

dds <- DESeqDataSetFromMatrix(countData = as.matrix(rna_seq_data), colData = col_data, design = ~condition)

#Normalize the data
dds <- DESeq(dds)
dds_n <- scale(rna_seq_data)

# regularized log transformation
rld <- rlog(dds) 
rld_n <- log2 

#extract the transformed data matrix
rld_mat <- assay(rld) 
pca_res <- prcomp(t(rld_mat))
pca_res_n <- prcomp(t(rld_n))

# Create a data frame for ggplot
pca_df <- as.data.frame(pca_res$x)
pca_df_n <- as.data.frame(pca_res_n$x)

#can be skip if it already contain the data we want e.g. condition
pca_df_n$Sample <- rownames(pca_df)

pca_df_n <- merge(pca_df_n, col_data, by = "Sample")
pca_df$Sample <- rownames(pca_df)
pca_df <- merge(pca_df, col_data, by = "Sample")
boxplot(rna_seq_data)

# dev.off()
# Plot the PCA results
plot1 <- ggplot(pca_df, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4) +
  labs(title = "PCA of RNA-Seq Data",
       x = paste0("PC1: ", round(summary(pca_res)$importance[2, 1] * 100, 1), "% variance"),
       y = paste0("PC2: ", round(summary(pca_res)$importance[2, 2] * 100, 1), "% variance")) +
  theme_minimal()
plot2 <- ggplot(pca_df_n, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4) +
  labs(title = "PCA of RNA-Seq Data",
       x = paste0("PC1: ", round(summary(pca_res)$importance[2, 1] * 100, 1), "% variance"),
       y = paste0("PC2: ", round(summary(pca_res)$importance[2, 2] * 100, 1), "% variance")) +
  theme_minimal()
plot1 + plot2

pc1 <- pca_df[, 2]

#pearson correlation method
corr_matrix <- cor(t(pca_res$x), method = "pearson")
corr_matrix_n <- cor(t(pca_res_n$x), method = 'pearson')
corr_test <- cor()

View(corr_matrix)
heatmap(corr_matrix, symm = TRUE, col = colorRampPalette(c("blue", "white", "red"))(100))
heatmap(corr_matrix_n, symm = TRUE, col = colorRampPalette(c("blue", "white", "red"))(100))



#Using edgeR to find DEG
# Example data: counts matrix and sample information
counts <- matrix(rpois(20000, lambda = 10), nrow = 1000, ncol = 20)  # 1000 genes, 20 samples
group <- factor(c(rep("Control", 10), rep("Treatment", 10)))

# Create a DGEList object
col_data_group <- factor(c(rep("Naive", 4), rep("Allo_2H", 4), rep("Allo_24H", 4)))
y <- DGEList(counts = rna_seq_data, group = col_data_group)

#Filter Lowly Expressed Genes
keep <- filterByExpr(y)
y <- y[keep, , keep.lib.sizes = FALSE]

normalized_counts <- cpm(y, normalized.lib.sizes = TRUE)
# Apply TMM (trimmed mean of M-values) normalization
y <- calcNormFactors(y)

# Estimate the common, trended, and tagwise dispersions
y <- estimateDisp(y)

# Differential Expression Testing, Fit a generalized linear model (GLM) and perform likelihood ratio tests
fit <- glmFit(y)
lrt <- glmLRT(fit, coef = 2)  # coef=2 compares the second group to the first

# Extract DEG
topTags(lrt)

# Visualization
plotMD(lrt)  # MA plot

correlation_matrix_pearson <- cor(y$counts, method = "pearson")
correlation_matrix_pearson_n <- cor(normalized_counts, method = "pearson")
heatmap(correlation_matrix_pearson, symm = TRUE, col = colorRampPalette(c("blue", "white", "red"))(100))
heatmap(correlation_matrix_pearson_n, symm = TRUE, col = colorRampPalette(c("blue", "white", "red"))(100))

