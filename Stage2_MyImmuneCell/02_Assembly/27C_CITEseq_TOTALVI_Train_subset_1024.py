#!/usr/bin/env python
# coding: utf-8




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
scvi.settings.num_threads = 16
scvi.settings.seed = 0

sc._settings.n_jobs=16
sc.settings.verbosity = 3

scvi.settings.dl_num_workers=2
sc.logging.print_header()


# # 1. Run on NVIDIA A100
os.environ["CUDA_VISIBLE_DEVICES"] = "7"

Classification = sys.argv[1] 

dataset = sys.argv[2] #celltype

obj_path = f'/home/liyanguo/MyImmuCell/06_Finnal_raw_count/{Classification}/{dataset}/'

adata = sc.read_h5ad(f"{obj_path}{dataset}_scRNA_count_HVG.h5ad")
adt = sc.read_h5ad(f"{obj_path}{dataset}_scADT_count.h5ad")

mdata = md.MuData({"rna": adata, "protein": adt})
mdata.update()

print(f"Max value of protein counts that store in X: {mdata.mod['protein'].X.max()}")
print(f"Max value of gene counts that store in X: {mdata.mod['rna'].X.max()}")

scvi.model.TOTALVI.setup_mudata(
    mdata,
    rna_layer=None,
    protein_layer=None,
    batch_key="Batch",
    modalities={
        "rna_layer": "rna",
        "protein_layer": "protein",
        "batch_key": "protein",
    },
)

del adata
del adt
# ## 1024 batch size,lr 0.01

model = scvi.model.TOTALVI(mdata,n_latent=20,
                           n_hidden=256, n_layers_decoder=1)

model.train(accelerator="auto",lr=0.01,batch_size=1024,
            max_epochs=200)
model.save(obj_path, overwrite=True, prefix=f'{dataset}_CITEseq_TOTALVI_')

# # 2. Get_latent_representation
X_totalVI = model.get_latent_representation()

np.save(f"{obj_path}{dataset}_X_TOTALVI",X_totalVI)


# # 3. Get model parameter

sheet=pd.DataFrame()
for key in model.history.keys():
    temp = model.history[key]
    temp = temp.reset_index()
    sheet=pd.concat([sheet,temp],axis=1)

pd.DataFrame.to_csv(sheet,f"{obj_path}{dataset}_TOTALVI_model_history.csv")

fig, ax = plt.subplots(1, 1)
sheet["elbo_train"].plot(ax=ax, label="train")
sheet["elbo_validation"].plot(ax=ax, label="validation")
ax.set(title="Negative ELBO over training epochs", ylim=(0, 1400))
ax.legend()
plt.savefig(f"{obj_path}{dataset}_ELBO.png")

fig, ax = plt.subplots(1, 1)
sheet["train_loss_epoch"].plot(ax=ax, label="train")
sheet["validation_loss"].plot(ax=ax, label="validation")
ax.set(title="loss over training epochs", ylim=(0, 1400))
ax.legend()
plt.savefig(f"{obj_path}{dataset}_loss.png")