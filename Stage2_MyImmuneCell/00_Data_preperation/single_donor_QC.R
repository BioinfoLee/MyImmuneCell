## The doublets  distribution, UMAP, neutrophil marker and QC results of each sample -----------------------------------------------
single_donor_QC =  function(object,type){
  library(ggplot2)
  library(cowplot)
  library(viridis)
  pal <- viridis(n = 30, option = "H")

  p3 = DimPlot_scCustom(object,group.by = "scDblFinder.class",reduction = "umap",aspect_ratio = 1,colors_use = c("#00BFFF", "#FF7256"))+
    labs(subtitle = paste0("Rate:",round(table(object$scDblFinder.class)[1]/table(object$scDblFinder.class)[2],3)),tag = "A")+
    theme(title = element_text(size = 10))

  p4=object@meta.data %>%
    ggplot(aes(x=nCount_RNA, y=nFeature_RNA, col=percent_mito)) +
    geom_point(size=0.5,alpha=0.2) +
    scale_color_gradient(low = "#36648B", high = "#CD3333") +
    scale_x_log10() +
    scale_y_log10() +
    theme_classic() +
    geom_vline(xintercept = 500) +
    geom_hline(yintercept = 500) +
    coord_fixed()+
    theme(title = element_text(size = 6))+
    labs(subtitle = paste0("sample_ID:",Project(object),"\nCell Number:",length(object$SampleID),
        "\nxintercept = 500","\nyintercept = 500"),
         tag = "B")
  p4 = ggExtra::ggMarginal(
    p = p4,
    type = 'histogram',
    margins = 'both',
    size = 3.5,
    colour = 'gray98',
    fill = 'gray74',
    bins=40
  )
  p5 = FeaturePlot_scCustom(object,"CD3D",reduction = "umap",aspect_ratio = 1,min.cutoff="q10",max.cutoff="q90",colors_use = c("#C6E2FF","#EFE2AA", "#CD4F39"))+
    labs(subtitle ="T",tag = "C")+theme(title = element_text(size = 10))
    
  p6 = FeaturePlot_scCustom(object,"CD14",reduction = "umap",aspect_ratio = 1,min.cutoff="q10",max.cutoff="q90",colors_use = c("#C6E2FF","#EFE2AA", "#CD4F39"))+
    labs(subtitle ="Monocyte",tag = "D")+theme(title = element_text(size = 10))
    
  p7 = FeaturePlot_scCustom(object,"CSF3R",reduction = "umap",aspect_ratio = 1,min.cutoff="q10",max.cutoff="q90",colors_use = c("#C6E2FF","#EFE2AA", "#CD4F39"))+
    labs(subtitle ="Neutrophil",tag = "E")+theme(title = element_text(size = 10))
  
  #有的样本没有CEACAM8
  if ("CEACAM8" %in% rownames(object)) {
    p8 = FeaturePlot_scCustom(object,"CEACAM8",reduction = "umap",aspect_ratio = 1,min.cutoff="q10",max.cutoff="q90",colors_use = c("#C6E2FF","#EFE2AA", "#CD4F39"))+
      labs(subtitle ="Neutrophil-CEACAM8",tag = "F")+theme(title = element_text(size = 10))
    ggsave2(paste0("02_Read_QC/MyImmuCell_QC/",type,"_",Project(object),"_QC.png"),
            p3+p4+p5+p6+p7+p8,width = 11,height = 9)
  }else{
    ggsave2(paste0("02_Read_QC/MyImmuCell_QC/",type,"_",Project(object),"_QC.png"),
            p3+p4+p5+p6+p7,width = 11,height = 9)
  }
}
