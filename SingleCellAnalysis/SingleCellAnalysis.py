import numpy as np
import pandas as pd
import scanpy as sc

scdr = sc.read_h5ad('scdr_preprocessed.h5ad')
#little unknown bug in scanpy
scdr.uns['log1p'] = {'base':None}

#PCA
sc.tl.pca(scdr, svd_solver='arpack')

#PCA Plot
sc.pl.pca_variance_ratio(scdr, log=True)
#PCA result is saved in obsm

