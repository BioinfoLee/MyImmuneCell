#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import sys
import tempfile
import pandas as pd
import ray
import scanpy as sc
import scvi
import mudata as md
import muon
import seaborn as sns
import matplotlib.pyplot as plt
import torch
from ray import tune
from scvi import autotune
import celltypist
from sklearn_ann.kneighbors.annoy import AnnoyTransformer
print("Last run with scvi-tools version:", scvi.__version__)
os.environ["TUNE_DISABLE_STRICT_METRIC_CHECKING"] = "1"


# In[ ]:


scvi.settings.num_threads = 24
scvi.settings.seed = 0

sc._settings.ScanpyConfig.n_jobs=4
sc.settings.verbosity = 1
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=400,frameon=False,
    facecolor = 'white', figsize=(8,8), format='png')

sns.set_theme()
#torch.set_float32_matmul_precision("high")


# # 0. Loading and preprocessing the dataset

# In[ ]:


os.environ["CUDA_VISIBLE_DEVICES"] = "6"

obj_path = '/home/liyanguo/MyImmuCell/04_MyImmuCell_TOTALVI/'
sc.settings.figdir = obj_path
dataset = sys.argv[1] #high_nCount_RNA or low_nCount_RNA


adata = sc.read_h5ad(f"{obj_path}scRNA_MyImmuCell_{dataset}_HVG.h5ad",backed='r')
adt = sc.read_h5ad(f"{obj_path}scADT_MyImmuCell_{dataset}.h5ad",backed='r')

# downsample n cell in each cell type
ncells = int(sys.argv[2])
cell_index = celltypist.samples.downsample_adata(adata,mode = 'each', n_cells = ncells,
                                                 by = 'Reference_Atlas_L1L2_pl',
                                                 return_index = True,random_state=0)

adt = adt[cell_index].to_memory()

adata = adata[cell_index].to_memory()


# # 1. Model hyperparameter tuning with TOTALVI

# In[ ]:


adata.obsm["protein"] = adt.X.copy()
adata.X = adata.layers['counts'].copy()
del adata.layers


# In[ ]:


adata.obs['Batch'] = adt.obs['Batch']


# In[ ]:


adata


# In[ ]:


adata.write(f"{obj_path}Model_hyperparameter_TOTALVI_MyImmuCell_{dataset}.h5ad",compression="gzip")


# In[ ]:


search_space = {
    "model_params": {"n_layers_decoder": tune.choice([1, 2]),"n_hidden": tune.choice([128, 256]),"n_latent": tune.choice([20, 30])},
    "train_params": {"max_epochs": 150, 'batch_size':tune.choice([256, 512, 1024,4096]), "lr": tune.choice([1e-5,1e-4,5e-4,1e-3,4e-3,1e-2])},
}


# In[ ]:


ray.shutdown()


# In[ ]:


model_cls = scvi.model.TOTALVI

model_cls.setup_anndata(adata,
                        protein_expression_obsm_key='protein',
                        batch_key="Batch"
                       )

results = autotune.run_autotune(
    model_cls,
    data=adata,
    mode="min",
    seed=0,
    metrics=["validation_loss",'train_loss_step','train_loss_epoch','elbo_train','elbo_validation'],
    search_space=search_space,
    num_samples=90,
    resources={"cpu": 24, "gpu": 1},
)

print(results.result_grid)

model_cls.save(obj_path, overwrite=True, prefix=f'{dataset}_CITEseq_TOTALVI_')

# In[ ]:


ray.shutdown()