#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import os
import sys
import scanpy as sc
import scvi
import mudata as md
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import torch
print("Last run with scvi-tools version:", scvi.__version__)
scvi.settings.num_threads = 48
scvi.settings.seed = 0

sc._settings.n_jobs=4
sc.settings.verbosity = 3


# In[ ]:


scvi.settings.dl_num_workers=2
sc.logging.print_header()


# # 1. Run on NVIDIA A100

# In[ ]:


#os.environ["CUDA_VISIBLE_DEVICES"] = "7"


# In[ ]:


dataset = 'CEACAM8_Neg_Neutrophil'



# In[ ]:


obj_path = f'/home/liyanguo/MyImmuCell/05_MyImmuCell_subpopulation/Level2_Refine_R3/{dataset}/'


# In[ ]:


adata = sc.read_h5ad(f"{obj_path}{dataset}_HVG.h5ad")
adt = sc.read_h5ad(f"{obj_path}{dataset}_preprocess_scADT.h5ad",backed='r')
adata.obs['Batch']=adt.obs['Batch']

print(f"Max value of protein counts that store in X: {adata.X.max()}")

# In[ ]:


scvi.model.SCVI.setup_anndata(
    adata,
    batch_key ="Batch"
)


# In[ ]:


model = scvi.model.SCVI(adata, n_layers=2, n_latent=20)


# In[ ]:


model.train(accelerator="auto",batch_size=4096,
            max_epochs=200)


# In[ ]:


model.save(obj_path, overwrite=True, prefix=f'{dataset}_scRNA_scVI_')


# # 2. Get_latent_representation

# In[ ]:


X_scVI = model.get_latent_representation()


# In[ ]:


np.save(f"{obj_path}{dataset}_X_scVI",X_scVI)


# # 3. Get model parameter

# In[ ]:


sheet=pd.DataFrame()
for key in model.history.keys():
    temp = model.history[key]
    temp = temp.reset_index()
    sheet=pd.concat([sheet,temp],axis=1)


# In[ ]:


pd.DataFrame.to_csv(sheet,f"{obj_path}{dataset}_scVI_model_history.csv")


# In[ ]:


# fig, ax = plt.subplots(1, 1)
# sheet["elbo_train"].plot(ax=ax, label="train")
# ax.set(title="Negative ELBO over training epochs", ylim=(0, 2000))
# ax.legend()
# plt.savefig(f"{obj_path}{dataset}_scVI_ELBO.png")

#nohup python 12X_CITEseq_scVI_L2L3_Train_4096.py CEACAM8_Neg_Neutrophil &