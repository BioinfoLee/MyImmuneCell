Doublet_detection = function(object){
  library(scDblFinder)
  library(BiocParallel)
  #识别双细胞scDblFinder,no cluster参数，没有预先分群，要按单个样本执行
  sce_object  = as.SingleCellExperiment(object)
  sce_object = scDblFinder(sce_object)
  object$scDblFinder.class = as.character(sce_object$scDblFinder.class)
  return(object)
}
