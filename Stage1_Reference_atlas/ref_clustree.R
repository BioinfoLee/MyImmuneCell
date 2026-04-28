library(clustree)
library(cowplot)
args <- commandArgs(T)

data = read.csv(paste0(as.character(args[1]),'/Level2_scVI_ForClustree.csv'),row.names=1)
p=clustree(data,'L2_leiden_scVI_')
ggsave2(paste0(as.character(args[1]),"/L2_leiden_scVI_clustree.png"),p,width=15,height=10)
