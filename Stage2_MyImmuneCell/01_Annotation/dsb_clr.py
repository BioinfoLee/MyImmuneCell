# coding: utf-8
import sys
import os
import scanpy as sc
import mudata as md
import muon as mu
import anndata as ad
sc._settings.ScanpyConfig.n_jobs=4
sc.logging.print_header()


obj_path = '/home/liyanguo/MyImmuCell/02_Read_QC/'
sc.settings.figdir = obj_path

data_type = sys.argv[1]

adt = sc.read_h5ad(f"{obj_path}tmp_scADT_filtered_MyImmuCell_{data_type}.h5ad")

raw_adt = sc.read_h5ad(f"{obj_path}scADT_raw_MyImmuCell.h5ad",backed='r')
raw_adt = raw_adt[raw_adt.obs_names.difference(adt.obs_names),:].to_memory()
print('dsb')
#dsb
mu.prot.pp.dsb(adt,raw_adt,add_layer=True,denoise_counts=True)

print('save')
adt.write_h5ad(f"{obj_path}scADT_MyImmuCell_dsb_{data_type}.h5ad",compression="gzip")
print('done')