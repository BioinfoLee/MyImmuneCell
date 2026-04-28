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
import torch
print("Last run with scvi-tools version:", scvi.__version__)
scvi.settings.num_threads = 36
scvi.settings.seed = 0

sc._settings.ScanpyConfig.n_jobs=4
sc.settings.verbosity = 3


# In[2]:


scvi.settings.dl_num_workers=2
sc.logging.print_header()


# # 1. Run on NVIDIA A100

# In[3]:


os.environ["CUDA_VISIBLE_DEVICES"] = "7"


# In[4]:


dataset = sys.argv[1] #celltype


# In[5]:


#dataset='HSPC'


# In[6]:


obj_path = f'/home/liyanguo/MyImmuCell/05_MyImmuCell_subpopulation/Level2_Refine_R1/{dataset}/'


# In[7]:


adata = sc.read_h5ad(f"{obj_path}{dataset}_preprocess_scRNA.h5ad")
adt = sc.read_h5ad(f"{obj_path}{dataset}_preprocess_scADT.h5ad")


# In[8]:


mdata = md.MuData({"rna": adata, "protein": adt})
mdata.update()


# In[9]:


print(f"Max value of protein counts that store in X: {mdata.mod['protein'].X.max()}")


# In[10]:


print(f"Max value of protein counts that store in layers-counts: {mdata.mod['rna'].layers['counts'].max()}")


# In[11]:


scvi.model.TOTALVI.setup_mudata(
    mdata,
    rna_layer="counts",
    protein_layer=None,
    batch_key="Batch",
    modalities={
        "rna_layer": "rna",
        "protein_layer": "protein",
        "batch_key": "protein",
    },
)


# In[12]:


del adata


# In[13]:


model = scvi.model.TOTALVI(mdata,n_latent=30,
                           n_hidden=256, n_layers_decoder=2)


# In[14]:


model.train(accelerator="auto",lr=0.004,batch_size=4096,
            max_epochs=200)


# In[15]:


model.save(obj_path, overwrite=True, prefix=f'{dataset}_CITEseq_TOTALVI_')


# # 2. Get_latent_representation

# In[16]:


X_totalVI = model.get_latent_representation()


# In[17]:


np.save(f"{obj_path}{dataset}_X_TOTALVI",X_totalVI)


# # 3. Get model parameter

# In[18]:


sheet=pd.DataFrame()
for key in model.history.keys():
    temp = model.history[key]
    temp = temp.reset_index()
    sheet=pd.concat([sheet,temp],axis=1)


# In[19]:


pd.DataFrame.to_csv(sheet,f"{obj_path}{dataset}_TOTALVI_model_history.csv")


# In[20]:


fig, ax = plt.subplots(1, 1)
sheet["elbo_train"].plot(ax=ax, label="train")
sheet["elbo_validation"].plot(ax=ax, label="validation")
ax.set(title="Negative ELBO over training epochs", ylim=(0, 1400))
ax.legend()
plt.savefig(f"{obj_path}{dataset}_ELBO.png")


# In[21]:


fig, ax = plt.subplots(1, 1)
sheet["train_loss_epoch"].plot(ax=ax, label="train")
sheet["validation_loss"].plot(ax=ax, label="validation")
ax.set(title="loss over training epochs", ylim=(0, 1400))
ax.legend()
plt.savefig(f"{obj_path}{dataset}_loss.png")


# # 4. Denoise protein

# In[22]:


rna_denoised, protein_denoised = model.get_normalized_expression()


# In[23]:


adt.layers["denoised_protein"] = protein_denoised


# In[24]:


adt.layers["protein_foreground_prob"] = 100 * model.get_protein_foreground_probability()


# In[25]:


adt.write(f"{obj_path}{dataset}_preprocess_scADT.h5ad",compression="gzip")


# In[26]:


#shell
## nohup python 10B_CITEseq_TOTALVI_L2_Train.py Basophil &
## nohup python 10B_CITEseq_TOTALVI_L2_Train.py CEACAM8_Neg_Neutrophil &
## nohup python 10B_CITEseq_TOTALVI_L2_Train.py Platelet &

## nohup python 10B_CITEseq_TOTALVI_L2_Train.py HSPC &
## nohup python 10B_CITEseq_TOTALVI_L2_Train.py DC &
## nohup python 10B_CITEseq_TOTALVI_L2_Train.py CEACAM8_Pos_Neutrophil &
## nohup python 10B_CITEseq_TOTALVI_L2_Train.py pDC &
## nohup python 10B_CITEseq_TOTALVI_L2_Train.py Classical_Monocyte &
## nohup python 10B_CITEseq_TOTALVI_L2_Train.py Non_classical_Monocyte &
## nohup python 10B_CITEseq_TOTALVI_L2_Train_7.py CD56dimNK &
## nohup python 10B_CITEseq_TOTALVI_L2_Train_7.py CD56brightNK &

## nohup python 10B_CITEseq_TOTALVI_L2_Train_4.py B &
## nohup python 10B_CITEseq_TOTALVI_L2_Train_4.py gdT &

# nohup python 10B_CITEseq_TOTALVI_L2_Train.py MAIT &
# nohup python 10B_CITEseq_TOTALVI_L2_Train.py iNKT &
# nohup python 10B_CITEseq_TOTALVI_L2_Train_4.py CD4T &
# nohup python 10B_CITEseq_TOTALVI_L2_Train_6.py CD8T &

