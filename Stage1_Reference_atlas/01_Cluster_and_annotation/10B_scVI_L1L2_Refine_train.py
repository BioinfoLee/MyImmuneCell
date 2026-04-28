#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import sys
import scanpy as sc
import scvi
import mudata as md
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
print("Last run with scvi-tools version:", scvi.__version__)
scvi.settings.num_threads = 72
scvi.settings.seed = 0
sc._settings.ScanpyConfig.n_jobs=8
sc.settings.verbosity = 1


# In[2]:


sc.logging.print_header()


# In[3]:


celltype = sys.argv[1]

# In[5]:


obj_path = f'/home/liyanguo/MyImmuCell/05_Ref_Atlas_subpopulation/Level2_Refine_R1/{celltype}/'


# # 1. Run on GPU or CPU

# In[6]:


adata = sc.read_h5ad(f"{obj_path}{celltype}_preprocess_scRNA.h5ad")


# In[7]:


scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
    batch_key ="Batch",
    continuous_covariate_keys=["nCount_RNA", "percent_mito"]
)


# In[8]:


model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")


# In[ ]:


#notic 6M cell take 8h with batch size 4096 in CPU device
model.train(accelerator="cpu",
            max_epochs=400,batch_size=4096)


# In[ ]:


model.save(obj_path, overwrite=True, prefix=f'{celltype}_scRNA_scVI_')


# # 2. Get_latent_representation

# In[ ]:


X_scVI = model.get_latent_representation()


# In[ ]:


np.save(f"{obj_path}{celltype}_X_scVI",X_scVI)


# # 3. Get model parameter

# In[ ]:


sheet=pd.DataFrame()
for key in model.history.keys():
    temp = model.history[key]
    temp = temp.reset_index()
    sheet=pd.concat([sheet,temp],axis=1)


# In[ ]:


pd.DataFrame.to_csv(sheet,f"{obj_path}{celltype}_model_history.csv")


# In[ ]:


fig, ax = plt.subplots(1, 1)
sheet["elbo_train"].plot(ax=ax, label="train")
ax.set(title="Negative ELBO over training epochs", ylim=(0, 2000))
ax.legend()
plt.savefig(f"{obj_path}{celltype}_ELBO.png")

