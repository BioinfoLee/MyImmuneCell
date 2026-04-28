celltypist_prediction = function(object){
  reticulate::use_python("/home/liyanguo/anaconda3/envs/R/bin/python")
  library(reticulate)
  scanpy <- import("scanpy")
  celltypist <- import("celltypist")
 
  adata <- AnnData(
    X = t(LayerData(object, "data"))
  )

# Our reference atlas celltypist models
  Reference_Atlas_L1L2 = celltypist$annotate(filename = adata,majority_voting="True",
                                    model = 'model_Ref_Atlas_L1_L2_HVG.pkl')
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
#AIFI proj
  AIFI_L1 = celltypist$annotate(filename = adata,majority_voting="True",
                                    model = 'ref_pbmc_clean_celltypist_model_AIFI_L1_2024-04-18.pkl')
  AIFI_L2 = celltypist$annotate(filename = adata,majority_voting="True",
                                    model = 'ref_pbmc_clean_celltypist_model_AIFI_L2_2024-04-19.pkl')
  AIFI_L3 = celltypist$annotate(filename = adata,majority_voting="True",
                                    model = 'ref_pbmc_clean_celltypist_model_AIFI_L3_2024-04-19.pkl')

  celltypist_predictions = data.frame(Reference_Atlas_L1L2_mv = Reference_Atlas_L1L2$predicted_labels$majority_voting,
                                      Reference_Atlas_L1L2_pl = Reference_Atlas_L1L2$predicted_labels$predicted_labels,
                                      
                                      AIFI_L1 = AIFI_L1$predicted_labels$majority_voting,
                                      AIFI_L2 = AIFI_L2$predicted_labels$majority_voting,
                                      AIFI_L3 = AIFI_L3$predicted_labels$majority_voting,
                                      
                                      Immune_All_High = Immune_All_High$predicted_labels$majority_voting,
                                      Immune_All_Low = Immune_All_Low$predicted_labels$majority_voting,
                                      
                                      Adult_Human_Blood=Adult_Human_Blood$predicted_labels$majority_voting,
                                      Adult_Human_Bone_marrow=Adult_Human_Bone_marrow$predicted_labels$majority_voting,
                                      
                                      row.names = rownames(Immune_All_High$predicted_labels))
  return(celltypist_predictions)
}