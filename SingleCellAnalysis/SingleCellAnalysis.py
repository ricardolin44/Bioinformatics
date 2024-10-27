import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

scdr = sc.read_h5ad('./output_data/scdr_preprocessed.h5ad')
#little unknown bug in scanpy
scdr.uns['log1p'] = {'base':None}

#PCA
sc.tl.pca(scdr, svd_solver='arpack')

#PCA Plot
sc.pl.pca_variance_ratio(scdr, log=True)
#PCA result is saved in obsm

# 2. Find neighbors and clusters
sc.pp.neighbors(scdr, n_pcs=10)  # Using the first 10 PCs
sc.tl.leiden(scdr, resolution=0.5)  # Clustering

# Display the first 5 cluster assignments
print(scdr.obs['leiden'].head())

# 3. Run UMAP for dimensionality reduction and visualization
sc.tl.umap(scdr)
sc.pl.umap(scdr, color='leiden')  # Coloring by cluster

# 4. Heatmap of genes for a specific component
sc.pl.heatmap(scdr, var_names=scdr.var_names[:10], groupby='leiden', use_raw=True)

# 5. Find markers for a specific cluster (e.g., cluster 1)
cluster1_markers = sc.tl.rank_genes_groups(scdr, groupby='leiden', groups=['1'], reference='rest', method='wilcoxon')
sc.pl.rank_genes_groups_violin(scdr, groups=['1'], n_genes=5)  # Violin plot for top markers in cluster 1

# 6. Find all markers for each cluster
sc.tl.rank_genes_groups(scdr, groupby='leiden', method='wilcoxon')
marker_genes = scdr.uns['rank_genes_groups']
groups = marker_genes['names'].dtype.names  # Cluster names
# Create a DataFrame from all ranked marker genes, including cluster labels
markers_df = pd.DataFrame({'cluster': [group for group in groups for _ in range(len(marker_genes['names'][group]))],'gene': [gene for group in groups for gene in marker_genes['names'][group]]})
# Get top 2 markers for each cluster
top_markers = markers_df.groupby('cluster').head(2)

# 7. Feature plots for specific genes
sc.pl.umap(scdr, color=top_markers['gene'][:4])  # First 4 markers
sc.pl.umap(scdr, color=top_markers['gene'][4:8])  # Next 4 markers

