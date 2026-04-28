celltypist_prediction = function(object){
  reticulate::use_python("/home/liyanguo/anaconda3/envs/R/bin/python")
  library(reticulate)
  scanpy <- import("scanpy")
  celltypist <- import("celltypist")
  #需要提前复制一下model: 
    #cp -r /localdisk/immune/softerware/celltypist/models ~/.celltypist/data
    #cp -r /home/liyanguo/OneKBloodMap/07_ref_model/celltypist_organ/* ~/.celltypist/data/models

  #> Show all available models that can be downloaded and used.使用其他软件下载，传到如下目录中
  #Show the local directory storing these models.
  # celltypist$models$models_path
  #Get an overview of the models
  # celltypist$models$models_description()

  #Select the model from the above list. If the `model` argument is not provided, will default to `Immune_All_Low.pkl`.
  # model = celltypist$models$Model$load(model = 'Immune_All_Low.pkl')
  #The model summary information.
  # model
  #Examine cell types contained in the model.
  # model$cell_types
  #Examine genes/features contained in the model.
  # model$features

  #Since the expression of each gene will be centred and scaled by matching with the mean and standard deviation of that gene in the provided model, CellTypist requires a logarithmised and normalised expression matrix stored in the AnnData (log1p normalised expression to 10,000 counts per cell). 
  adata <- AnnData(
    X = t(LayerData(object, "data"))
  )
  # scanpy$pp$filter_cells(adata, min_genes=200)
  # scanpy$pp$filter_genes(adata, min_cells=3)
  # scanpy$pp$normalize_total(adata,target_sum = 1e4)
  # scanpy$pp$log1p(adata)
  # Predict the identity of each input cell.

# celltypist models
  Immune_All_High = celltypist$annotate(filename = adata,majority_voting="True",
                                    model = 'Immune_All_High.pkl')
  Immune_All_Low = celltypist$annotate(filename = adata,majority_voting="True",
                                    model = 'Immune_All_Low.pkl')
# CellHint auto-harmonised and auto-integrated atlas
  Adult_Human_Blood = celltypist$annotate(filename = adata,majority_voting="True",
                                    model = 'Adult_Human_Blood.pkl')
  Adult_Human_Bone_marrow = celltypist$annotate(filename = adata,majority_voting="True",
                                    model = 'Adult_Human_Bone_marrow.pkl')
#
  AIFI_L1 = celltypist$annotate(filename = adata,majority_voting="True",
                                    model = 'ref_pbmc_clean_celltypist_model_AIFI_L1_2024-04-18.pkl')
  AIFI_L2 = celltypist$annotate(filename = adata,majority_voting="True",
                                    model = 'ref_pbmc_clean_celltypist_model_AIFI_L2_2024-04-19.pkl')
  AIFI_L3 = celltypist$annotate(filename = adata,majority_voting="True",
                                    model = 'ref_pbmc_clean_celltypist_model_AIFI_L3_2024-04-19.pkl')

  celltypist_predictions = data.frame(Immune_All_High = Immune_All_High$predicted_labels$majority_voting,
                                      Immune_All_Low = Immune_All_Low$predicted_labels$majority_voting,
                                      
                                      Adult_Human_Blood=Adult_Human_Blood$predicted_labels$majority_voting,
                                      Adult_Human_Bone_marrow=Adult_Human_Bone_marrow$predicted_labels$majority_voting,
                                      
                                      AIFI_L1 = AIFI_L1$predicted_labels$majority_voting,
                                      AIFI_L2 = AIFI_L2$predicted_labels$majority_voting,
                                      AIFI_L3 = AIFI_L3$predicted_labels$majority_voting,
                                      row.names = rownames(Immune_All_High$predicted_labels))
  return(celltypist_predictions)
}