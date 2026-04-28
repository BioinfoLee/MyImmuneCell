#!/usr/bin/env python
# coding: utf-8

# In[ ]:


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
scvi.settings.num_threads = 48
scvi.settings.seed = 0
sc._settings.ScanpyConfig.n_jobs=48
sc.settings.verbosity = 4
sc.settings.set_figure_params(dpi=100, fontsize=10, dpi_save=400,
    facecolor = 'white', figsize=(8,8), format='png')
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 100)


# In[ ]:


sc.logging.print_header()


# In[ ]:


dataset = sys.argv[1] #celltype


# In[ ]:


#dataset='HSPC'#cell type


# In[ ]:


obj_path = f'/home/liyanguo/MyImmuCell/05_MyImmuCell_subpopulation/Level2_Refine_R4/{dataset}/'
sc.settings.figdir = obj_path


# In[ ]:


adata = sc.read_h5ad(f"{obj_path}{dataset}_preprocess_scRNA.h5ad")


# # 1. Load latent representation

# In[ ]:


adata.obsm["X_TOTALVI"]=np.load(f"{obj_path}{dataset}_X_TOTALVI.npy")


# In[ ]:


# sc.pp.neighbors take about 60 minutes for 3000M cells
print(f"Do sc.pp.neighbors with AnnoyTransformer,{datetime.now()}")
sc.pp.neighbors(adata, transformer=AnnoyTransformer(20), use_rep='X_TOTALVI')


# In[ ]:


print(f"Run UMAP,{datetime.now()}")
sc.tl.umap(adata, min_dist=0.5)


# In[ ]:


def run_leiden_parallel(params):
    global adata
    key_added = f"L4_leiden_TOTALVI_{params['resolution']}"
    sc.tl.leiden(adata, key_added=f'L4_leiden_TOTALVI_{params['resolution']}', 
                 **params,
                 use_weights=True,directed=False,flavor="igraph")
    return {key_added: adata.obs[key_added]}


# In[ ]:


param = [
        {"resolution": 0.1},
        {"resolution": 0.5},
        {"resolution": 1},
        {"resolution": 1.5},
        {"resolution": 0.3},
        {"resolution": 0.8},
]


# In[ ]:


with Pool(6) as p:
    leiden_columns = p.map(run_leiden_parallel, param)


# In[ ]:


for col_dict in leiden_columns:
    for key, values in col_dict.items():
        adata.obs[key] = values


# In[ ]:


del adata.obsp
del adata.layers


# In[ ]:


adata.write(f"{obj_path}{dataset}_HVG.h5ad",compression="gzip")


# In[ ]:


leiden_data = adata.obs.loc[:,['L4_leiden_TOTALVI_0.1','L4_leiden_TOTALVI_0.5', 'L4_leiden_TOTALVI_1', 'L4_leiden_TOTALVI_1.5', 'L4_leiden_TOTALVI_0.3','L4_leiden_TOTALVI_0.8']]
obsm = adata.obsm
leiden_groups=['L4_leiden_TOTALVI_0.1','L4_leiden_TOTALVI_0.5', 'L4_leiden_TOTALVI_1', 'L4_leiden_TOTALVI_1.5', 'L4_leiden_TOTALVI_0.3','L4_leiden_TOTALVI_0.8']


# In[ ]:


leiden_data.to_csv(f'{obj_path}/Level4_TOTALVI_ForClustree.csv')
os.system(f'/home/liyanguo/anaconda3/envs/R/bin/Rscript clustree_L4.R {obj_path} &')


# # 2. Get all gene count data in backed mode
#adata = sc.read_h5ad(f"{obj_path}{dataset}_HVG.h5ad",backed='r')
# In[ ]:


adata_raw = sc.read_h5ad(f"{obj_path}{dataset}_count_scRNA.h5ad",backed='r')


# # 3. Plot leiden umap

# In[ ]:


print(f"Plot umap,{datetime.now()}")
sc.pl.umap(
    adata,
    color=leiden_groups,
    legend_fontsize=6,legend_loc='on data',ncols=3,frameon=False,
    save=f'_{dataset}_TOTALVI_L4_annoy_leiden',show=False
)


# In[ ]:


print(f"Plot umap,{datetime.now()}")
sc.pl.umap(
    adata,
    color=['AIFI_L2','AIFI_L3','Immune_All_Low','Adult_Human_Blood','predicted.celltype.l2','scDblFinder.class'],
    legend_fontsize=6,legend_loc='on data',ncols=3,frameon=False,
    save=f'_{dataset}_auto_annotation',show=False
)


# # 4. Plot canonical and reference atlas marker genes

# In[ ]:


marker_dict1={
    'Immune cell': ['PTPRC'],

    'Lineage':['CD7','CD3E',
               'IL7R',
               'SPON2','KLRF1',
               'CD79A','MS4A1',
               'CD14','FCGR3A',
               'CSF3R','FCGR3B'
              ],
    
    'HSPC':['CD34','SPINK2','CYTL1','PROM1','SMIM24','EGFL7',
            'SOX4','KIT', 'DNTT','ETV6','MCL1','STAT5A','CD48'],
    'Basophil':['HDC','GATA2', 'FCER1A', 'IL3RA','ENPP3','MS4A2','IL4'],
    'Eosinophil':['ALOX15','SIGLEC10','SIGLEC8','LMO4','EPX','ITGA1','CCR3',],
    'Mast cell':['FCER1A','MS4A2','KIT'],

    'CEACAM8- Neutrophil':['FCGR3B','CSF3R','MME','G0S2','MNDA',],
    'CEACAM8+ Neutrophil':['CEACAM8','LTF','BPI','MPO','ELANE'],
    
    'Classical monocytes' :['CD14','LYZ','VCAN','FCN1',],
    'Non-classical monocytes' :['FCGR3A','CDKN1C','TCF7L2','CSF1R',],
    'DC':['ENHO','CD1C','HLA-DQA1','FCER1A','CLEC10A','CLEC9A','FLT3'],
    'pDC':['CLEC4C','IL3RA'],
    
    'T naive': ['CD3D','CD3G','LEF1','IL7R','CCR7', 'TCF7','SELL','BACH2','ITGB1'],#ITGB1-
    'CD4+ T': ['CD4','MAL','RCAN3','IL6ST','TRAT1','CAMK4'],
    'TRAV1-2- CD8+ T': ['CD8A', 'CD8B','LINC02446','CCL5','GZMH','GZMK','ZNF683','THEMIS'],
    'γδ T':['TRDV2','TRGV9','TRDC','KLRG1','TRGC1','TRGC2','CCL5','CST7','TRDV1'],
    'MAIT':['SLC4A10','KLRB1','TRAV1-2','RORA','CXCR6'],
    'Treg':['FOXP3','CTLA4','IL2RA','TIGIT','RTKN2','STAM'],

    'NK': ['KLRD1','GNLY','PRF1','GZMB','CD244','CD247',
           'IL2RB','XCL1','XCL2'],
    'CD56dim':['FCGR3A','SPON2',],#CD16+
    'CD56bright':['IL7R','GZMK','NCAM1','GATA3'],#CD16-/low
    'Adaptive NK':['KLRC2','B3GAT1',],

    'Lymphocyte relate':['RUNX1','RUNX2','KLRF1','KLRC1',#KLRC1=NKG2AA抑制   NKp80=KLRF1
                         'KLRD1','KLRG1','EOMES',
                        'NCR2','NCR3','SYNE1'],#KLRD1=CD94抑制 KLRG1抑制

    'ILC': ['IL7R','LTB','RGS1','TNFSF10','PLCG2','RUNX1','TOX','FLT3'],#lin-
    'pILC':['NFIL3',],
    'ILC1':['IL2RB','TBX21','NCR1',],#CD56-
    'ILC2':['GATA3','IL2RA','PTGDR2','KLRG1','ZBTB16','RORA','MAF'],#lin- CD4-
    'ILC3-NCR+':['AHR','KIT','RORC','NCR1','NCR2','RUNX2',],#mature ILC #lin- CD56+ CD4-
    'LTi':['CCR7','KIT','ITGA4',],#mature ILC #lin- CD56-
    'ILCreg':['SOX4',],#KLRG1-
    
    'NKT':['TRAV24','TRBV28','CLDND1'], 
    'DN T':['FXYD2','NUCB2','MYB',],
    'T Develop':['SPI1','ZBTB7B','RUNX3'],
    'T activation': ['CD69', 'CD38'],

    'B naive': ['CD79A','TCL1A',],
    'Transitional B': ['MME','CD38','CD9','FCER2','EBF1','BACH2','PAX5', 'MSI2',],
    'Atypical memory B': ['TBX21','ITGAX','FCRL5','SIGLEC6','SOX5'],
    'Plasma':['JCHAIN','IGHA1','TNFRSF17','MZB1','CD38','XBP1', 'PRDM1'],
    'CD5+ B':['IGHV3-7','LEF1','CTLA4','CD5'],

    'Platelet':['PF4','PPBP','GP9','CAVIN2'],
    'RBC':['HBB','HBA2','BLVRB'],
    'Proliferative signal':['MKI67','TOP2A','STMN1'],
    'Other':['ITGAM','ANK3','PRKCA','BCL6','TGFB1','CXCR4','GATA1','ATXN1','SELPLG',
             'ITGA6','CSF1','CSF1R','CSF2','ETS1','CD44','IRF4','STAT3','BATF','TGFBR1','TGFBR2']
}
#Protein :       CD114,  CD116,  CD123, CD124, CD125,  CD126, CD203c,  KLRB1,  CD25,   CD29,   CD73,  CD127,  CD119,  CRTH2,  PLZF, CD62L, TdT, CD57,     ThpOK,   PU.1, CRTH2
#Gene Symbol :	CSF3R,   CSF2RA, IL3RA, IL4R,  IL5RA,  IL6R, ENPP3,   CD161,  IL2RA,  ITGB1,  NT5E,  IL7R,  IFNGR1,  PTGDR2, ZBTB16, SELL, DNTT, B3GAT1, ZBTB7B,   SPI1, PTGDR2

#Protein :    CD42a, CD42b, 
#Gene Symbol :GP9, GP1BA

#TCR : Vα7.2
#Gene : TRAV1-2

#GCSF=CSF3, for neutrophils
#GM-CSF=CSF2, granulocyte-macrophage colony-stimulating factor
#M-CSF=CSF1, macrophage colony-stimulating factor


# In[ ]:


marker_dict2={
    'B': ['MS4A1','CD79A','BANK1','CD22','FCRL1','TCL1A','FCER2'],
    'Plasma|Plasmablast':['TNFRSF17','JCHAIN','IGHA1','MZB1','XBP1'],
    
    'CD4+ T':['TCF7','MAL','IL7R','RCAN3','DGKA',],
    'CD8+ T':['CD8A','CD8B','LINC02446','CCL5','THEMIS','KLRK1'],
    'MAIT':['SLC4A10','TRAV1-2','KLRB1',],
    'Platelet':['TUBB1','CAVIN2','PF4',],
    'NKT':['TRAV24','TRBV28','CLDND1','DZIP3','ORM2','GPRC5D','FGFR2'],
    'NK': ['KLRF1', 'SPON2','GNLY', 'GZMB','PRF1','KLRD1','TYROBP'],
    'non-NK ILC': ['FSTL4', 'HPGDS','PKIB', 'PTGDR2','SCART1'],
    'γδ T':['TRDC','TRDV2','TRGC1','TRGC2','TRGV9','KLRD1','CCL5','KLRC1','TRDV1'],
    
    'HPSC':['CD34','NPR3','SMIM24','SPINK2','CRHBP','EMP1','ACY3','NYNRIN'], 
    'DC':['ENHO','CLEC10A','CD1C','FCER1A','FLT3','HLA-DQA1'],
    'Classical monocyte':['FCN1','VCAN','LRP1','CD14','MS4A6A','CST3'],
    'Non-classical monocyte':['CDKN1C','CSF1R','CKB','NEURL1','LYPD2','AIF1'],
    'CEACAM8+ Neutrophil':['CEACAM8','LCN2','BPI','LTF','CRISP3','CEACAM6'],
    'CEACAM8- Neutrophil':['CSF3R','NAMPT','FCGR3B','IFITM2','SLC25A37','S100A9'],
    'Eosinophil':['EPX', 'SIGLEC8', 'IL5RA'],
    'Basophil':['HDC','GATA2','MS4A2','AKAP12','CLC','IL4','RHOXF1P1',], 
    'Mast':['STXBP6', 'PPIL6', 'KIT',],
    'pDC':['SERPINF1','PTPRS','LILRA4','TNFRSF21','PLD4','SCT'],
    
    'Proliferative signal':['MKI67','TOP2A','STMN1'],
}


# In[ ]:


all_items = []
for values in marker_dict1.values():
    all_items.extend(values)
dict1 = np.unique(np.array(all_items))


# In[ ]:


all_items = []
for values in marker_dict2.values():
    all_items.extend(values)
dict2 = np.unique(np.array(all_items))


# In[ ]:


marker_selected = np.unique(
    np.hstack((dict2,dict1))
)


# In[ ]:


len(marker_selected)


# ## Get marker_selected count data

# In[ ]:


adata_raw = adata_raw[adata.obs_names,adata_raw.var_names.isin(marker_selected)]


# In[ ]:


name = f"tmp_{dataset}_{str(uuid.uuid4())}.h5ad"
adata_raw.copy(f"/home/liyanguo/MyImmuCell/tmp/{name}")


# In[ ]:


adata_raw.file.close()


# In[ ]:


adata = sc.read(f"/home/liyanguo/MyImmuCell/tmp/{name}")


# In[ ]:


adata.obs = adata.obs.join(leiden_data, how='left')


# In[ ]:


try:
    print(f"{adata.X.max()}")
except:
    print(adata.isbacked)


# In[ ]:


os.system(f"rm /home/liyanguo/MyImmuCell/tmp/{name}")


# ## normalize and plot

# In[ ]:


sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)


# In[ ]:


for group in leiden_groups:
    sc.pl.dotplot(
        adata,
        groupby=group,
        save=f"{dataset}_{group}_marker",
        var_names=marker_dict1,
        #standard_scale="var",
        dot_min=0.1,
        show=False
    )


# In[ ]:


for group in leiden_groups:
    sc.pl.dotplot(
        adata,
        groupby=group,
        save=f"{dataset}_{group}_ref_markers",
        var_names=marker_dict2,
        #standard_scale="var",
        dot_min=0.1,
        show=False
    )


# # 5. Process and Plot adt

# In[ ]:


adt = sc.read_h5ad(f"{obj_path}{dataset}_preprocess_scADT.h5ad")


# ## add adata obsm and leiden

# In[ ]:


for key in leiden_groups:
    try:
        del adt.obs[key]
    except:
        print("A")


# In[ ]:


adt.obs = adt.obs.join(leiden_data, how='left')


# In[ ]:


adt.obsm = obsm


# In[ ]:


adt.X.max()


# In[ ]:


adt.layers["clr"] = mu.prot.pp.clr(adt,inplace= False).X.copy()


# In[ ]:


adt.write(f"{obj_path}{dataset}_preprocess_scADT.h5ad",compression="gzip")


# ## plot

# In[ ]:


sc.pl.umap(
    adt,
    color=adt.var_names,
    legend_fontsize=6,legend_loc='on data',ncols=5,frameon=False,show=False,
    layer='clr',
    save=f'_{dataset}_scADT_TOTALVI_L4_annoy_leiden'
    )


# In[ ]:


adt_dict={
    "CD45RA": ["CD45RA"],
    
    'B': ['IgD','IgM',"CD19",'CD27','CD137'],
    
    "ILC|NK|NKT|γδ T|": ["CD56","CD161",'CD16',],
    
    "T": ["CD3","CD4","CD8",'CD27','CD28','CD134','CD272'],
    "Naive/Memory T": ["CD45RA","CD127","CD62L",'CD197',],
    "Treg":['CD25','CD137','CD357',],
    "Th1|ILC":["CD183",'CD186','CD366'],
    "Th17|ILC":['CD161','CD196',],
    "Tfh|ILC|Effector":['CD185','CD278','CD279',],

    'DC':['CD11c','HLA-DR',"CD183",'CD197'],

    'Monocytes':['CD14','CD16'],
}


# In[ ]:


for group in leiden_groups:
    sc.pl.dotplot(
    adt,
    groupby=group,
    var_names=adt_dict,
    save=f"{dataset}_{group}_scADT_TOTALVI_clr",layer='clr',
    expression_cutoff=1,
    #standard_scale="group",
    dot_min=0.1,
    show=False,
    )


# In[ ]:


for group in leiden_groups:
    sc.pl.dotplot(
    adt,
    groupby=group,
    var_names=adt_dict,
    save=f"{dataset}_{group}_scADT_TOTALVI_denoised_protein",layer='denoised_protein',
    expression_cutoff=1,
    #standard_scale="group",
    dot_min=0.1,
    show=False,
    )


# In[ ]:


import session_info
session_info.show()

