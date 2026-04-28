NormalizeData2Umap = function(scRNA_object){
  scRNA_object <- NormalizeData(scRNA_object) %>%
    FindVariableFeatures() %>% 
    ScaleData()  %>%
    RunPCA()
  pct<-scRNA_object[["pca"]]@stdev/sum(scRNA_object[["pca"]]@stdev)*100
  cumu<-cumsum(pct)
  co1<-which(cumu >80 & pct<5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  co2 <- min(co1,co2)
  scRNA_object <- RunUMAP(scRNA_object, dims = 1:co2,
                          seed.use = 1)
  return(scRNA_object)
}