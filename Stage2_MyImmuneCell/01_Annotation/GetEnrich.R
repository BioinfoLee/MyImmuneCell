# args <- commandArgs(T)
# FeatureMatrix = args[1]
# R_data_path = args[2]
# groupby = args[3]

GetEnrich = function (features = FeatureMatrix, pvalue_cutoff = 0.05, 
                      R_data_path=R_data_path,groupby=groupby,
          cols = c("#F47E5D", "#CA3D74", "#7F2880", "#463873"), plot_term_number = 5,
          gmtfile=NULL)
{
  suppressPackageStartupMessages({
      library(clusterProfiler)
  library(ggplot2)
  library(dplyr)})
  c5 <- read.gmt(gmtfile)

  ego_comb <- character(0)
  cluster.names <- colnames(FeatureMatrix)
  for (i in 1:(length(cluster.names))) {
    gene_list0 <- FeatureMatrix[, cluster.names[i]]
    ego0 <- enricher(gene_list0, TERM2GENE = c5)
    ego0_result <- ego0@result
    ego0_result$GeneRatio.num <- ego0_result$Count/length(gene_list0)
    ego0_result <- ego0_result[order(ego0_result$GeneRatio.num, 
                                     decreasing = T), ]
    ego0_result <- subset(ego0_result, ego0_result$pvalue < 
                            pvalue_cutoff)
    ego0_result <- cbind(ego0_result, cluster.names[i])
    colnames(ego0_result)[ncol(ego0_result)] <- "group"
    ego_comb <- rbind(ego_comb, ego0_result)
  }
  ego_comb$Description <- factor(ego_comb$Description, levels = rev(unique(ego_comb$Description)))
  ego_comb[1:3, ]
  ego_comb_sub <- ego_comb %>% group_by(group) %>% top_n(-plot_term_number, 
                                                         pvalue)
  pal <- colorRampPalette(cols)(100)
  ego_comb_sub$logp <- -log10(ego_comb_sub$pvalue)
  ego_comb_sub$Description <- factor(ego_comb_sub$Description, 
                                     levels = rev(unique(ego_comb_sub$Description)))
  ego_comb_sub$group <- factor(ego_comb_sub$group, levels = rev(unique(ego_comb_sub$group)))
  write.csv(ego_comb_sub,paste0(R_data_path,"/GetEnrich_",groupby,".csv"))
  p=ggplot(ego_comb_sub, aes(x = group, y = Description)) + 
    geom_point(aes(size = `GeneRatio.num`, color = logp)) + theme_bw(base_size = 14) + 
    scale_color_gradientn(colours = pal) + ylab(NULL) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ggtitle("")
  ggsave(paste0(R_data_path,"/GetEnrich_",groupby,".png"),p,width=10,height=10)
  p
}
