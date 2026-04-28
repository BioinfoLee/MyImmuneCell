args <- commandArgs(T)

suppressPackageStartupMessages({
  library(Seurat)
  library(anndata)
  library(reticulate)
})

use_python("/home/liyanguo/anaconda3/envs/R/bin/python")
sc <- import("scanpy")
celltypist <- import("celltypist")

adata = read_h5ad(paste0(args[1],'/',args[2],'_after_cluster_scVI.h5ad'))


print(paste0("Sample 100000 cells per group('",args[3],"')!"))
adata_sampled = celltypist$samples$downsample_adata(adata, mode = 'each', n_cells = 100000L, by = args[3], return_index = FALSE,random_state=0L)


scRNA_object = CreateSeuratObject(counts = t(adata_sampled$X), meta.data = adata_sampled$obs)
scRNA_object = NormalizeData(scRNA_object)
Idents(scRNA_object)=args[3]
markers_gene_seurat <- FindAllMarkers(scRNA_object,test.use = "roc",only.pos = TRUE,
                                      min.pct = 0.1,
                                      logfc.threshold =0.5,
                                      max.cells.per.ident=100000L)

write.csv(markers_gene_seurat,paste0(args[1],'/markers_seurat_roc',args[2],'_',args[3],'.csv'))

print(paste0("Done!"))