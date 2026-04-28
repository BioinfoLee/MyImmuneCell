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
scvi.settings.num_threads = 36
scvi.settings.seed = 0

sc._settings.ScanpyConfig.n_jobs=4
sc.settings.verbosity = 3


# In[ ]:


scvi.settings.dl_num_workers=2
sc.logging.print_header()


# # 1. Run on NVIDIA A100

# In[ ]:


os.environ["CUDA_VISIBLE_DEVICES"] = "5"


# In[ ]:


dataset = sys.argv[1] #celltype


# In[ ]:


#dataset='HSPC'


# In[ ]:


obj_path = f'/home/liyanguo/MyImmuCell/05_MyImmuCell_subpopulation/Level2_Refine_R4/{dataset}/'


# In[ ]:


adata = sc.read_h5ad(f"{obj_path}{dataset}_preprocess_scRNA.h5ad")
adt = sc.read_h5ad(f"{obj_path}{dataset}_preprocess_scADT.h5ad")


# In[ ]:


mdata = md.MuData({"rna": adata, "protein": adt})
mdata.update()


# In[ ]:


print(f"Max value of protein counts that store in X: {mdata.mod['protein'].X.max()}")


# In[ ]:


print(f"Max value of protein counts that store in layers-counts: {mdata.mod['rna'].layers['counts'].max()}")


# In[ ]:


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


# In[ ]:


del adata


# ## 1024 batch size，lr 0.01

# In[ ]:


model = scvi.model.TOTALVI(mdata,n_latent=20,
                           n_hidden=256, n_layers_decoder=1)


# In[ ]:


#high heterogeneity, 12 celltype
#B CytotoxicCD4 HSPC NaiveCD8 Treg TemCD8 CD4_helper_memory CD8Tcm NaiveCD4  MAIT iNKT
model.train(accelerator="auto",lr=0.01,batch_size=1024,
            max_epochs=200)

model.save(obj_path, overwrite=True, prefix=f'{dataset}_CITEseq_TOTALVI_')

# # 2. Get_latent_representation

# In[ ]:


X_totalVI = model.get_latent_representation()


# In[ ]:


np.save(f"{obj_path}{dataset}_X_TOTALVI",X_totalVI)


# # 3. Get model parameter

# In[ ]:


sheet=pd.DataFrame()
for key in model.history.keys():
    temp = model.history[key]
    temp = temp.reset_index()
    sheet=pd.concat([sheet,temp],axis=1)


# In[ ]:


pd.DataFrame.to_csv(sheet,f"{obj_path}{dataset}_TOTALVI_model_history.csv")


# In[ ]:


fig, ax = plt.subplots(1, 1)
sheet["elbo_train"].plot(ax=ax, label="train")
sheet["elbo_validation"].plot(ax=ax, label="validation")
ax.set(title="Negative ELBO over training epochs", ylim=(0, 1400))
ax.legend()
plt.savefig(f"{obj_path}{dataset}_ELBO.png")


# In[ ]:


fig, ax = plt.subplots(1, 1)
sheet["train_loss_epoch"].plot(ax=ax, label="train")
sheet["validation_loss"].plot(ax=ax, label="validation")
ax.set(title="loss over training epochs", ylim=(0, 1400))
ax.legend()
plt.savefig(f"{obj_path}{dataset}_loss.png")


# # 4. Denoise protein

# In[ ]:


rna_denoised, protein_denoised = model.get_normalized_expression()


# In[ ]:


adt.layers["denoised_protein"] = protein_denoised


# In[ ]:


adt.layers["protein_foreground_prob"] = 100 * model.get_protein_foreground_probability()


# In[ ]:


adt.write(f"{obj_path}{dataset}_preprocess_scADT.h5ad",compression="gzip")
