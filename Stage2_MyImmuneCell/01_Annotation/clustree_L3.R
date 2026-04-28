library(clustree)
library(cowplot)
args <- commandArgs(T)

data = read.csv(paste0(as.character(args[1]),'/Level3_TOTALVI_ForClustree.csv'),row.names=1)
p=clustree(data,'L3_leiden_TOTALVI_')
ggsave2(paste0(as.character(args[1]),"/L3_leiden_TOTALVI_clustree.png"),p,width=15,height=10)
