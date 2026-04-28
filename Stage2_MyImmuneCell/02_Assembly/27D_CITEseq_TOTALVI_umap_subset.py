#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
import scipy.io as sio
import scanpy.external as sce
import matplotlib.pyplot as plt
import re
import anndata as ad
import statistics
import torch
import scvi
import tempfile
import sklearn
import mudata as md
import muon as mu
md.set_options(pull_on_update=False)
import muon
from datetime import datetime
from scib_metrics.benchmark import Benchmarker
from sklearn_ann.kneighbors.annoy import AnnoyTransformer
from multiprocessing import Pool
import uuid
print("Last run with scvi-tools version:", scvi.__version__)
scvi.settings.num_threads = 16
scvi.settings.seed = 0
sc._settings.n_jobs=16
sc.settings.verbosity = 4
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=400,
    facecolor = 'white', figsize=(8,8), format='png')
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 100)

Classification = sys.argv[1]
dataset = sys.argv[2]

obj_path = f'/home/liyanguo/MyImmuCell/06_Finnal_raw_count/{Classification}/{dataset}/'
sc.settings.figdir = obj_path

adata = sc.read_h5ad(f"{obj_path}{dataset}_scRNA_count_HVG.h5ad")

adt = sc.read_h5ad(f"{obj_path}{dataset}_scADT_count.h5ad")

indices = pd.read_csv(f"{obj_path}/{dataset}_indices_labels.csv",index_col=0)

indices = indices[['Classification_L4','Classification_L3','Classification_L2','Classification_L1','leiden_cluster']]

if (adata.obs.index == indices.index).all():
    adata.obs = adata.obs.join(indices,how='left')
    adt.obs = adt.obs.join(indices, how='left')
else:
    raise ValueError('Indices error.')


# # 1. Load latent representation and umap

adata.obsm["X_TOTALVI"]=np.load(f"{obj_path}{dataset}_X_TOTALVI.npy")

print(f"Do sc.pp.neighbors with AnnoyTransformer,{datetime.now()}")
sc.pp.neighbors(adata, transformer=AnnoyTransformer(20), use_rep='X_TOTALVI')

print(f"Run UMAP,{datetime.now()}")
sc.tl.umap(adata, min_dist=0.5)

indices['UMAP_1'] = adata.obsm['X_umap'][:, 0]  #
indices['UMAP_2'] = adata.obsm['X_umap'][:, 1]  #

indices.to_csv(f"{obj_path}/{dataset}_indices_labels_umap.csv")

# # 2. Plot leiden umap

groupby='Classification_L4'

print(f"Plot umap,{datetime.now()}")
sc.pl.umap(
    adata,
    color=groupby,
    legend_fontsize=6,legend_loc='on data',ncols=3,frameon=False,
    save=f'_{dataset}_{groupby}',show=False
)

# # 4. Process adt

adt.obsm = adata.obsm

adt.X.max()

adt.layers["clr"] = mu.prot.pp.clr(adt,inplace= False).X.copy()

adt.write(f"{obj_path}{dataset}_scADT_count.h5ad",compression="gzip")

import session_info
session_info.show()