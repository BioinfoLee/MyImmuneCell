args <- commandArgs(T)

suppressPackageStartupMessages({
  library(Seurat)
  library(anndata)
  library(reticulate)
  library(openxlsx)
})

use_python("/home/liyanguo/anaconda3/envs/R/bin/python")
sc <- import("scanpy")
celltypist <- import("celltypist")

adata = read_h5ad(paste0(args[1],args[2],'_tmp.h5ad'))

scRNA_object = CreateSeuratObject(counts = t(adata$X), meta.data = adata$obs)
scRNA_object = NormalizeData(scRNA_object)
Idents(scRNA_object)=args[3]
markers_gene_seurat <- FindAllMarkers(scRNA_object,test.use = "roc",only.pos = TRUE,
                                      min.pct = 0.1,
                                      logfc.threshold =0.5,
                                      max.cells.per.ident=100000L)

write.xlsx(markers_gene_seurat, paste0(args[1],'FindAllMarkers_roc_',args[2],'_',args[3],'.xlsx'))

file.remove(paste0(args[1],args[2],'_tmp.h5ad'))