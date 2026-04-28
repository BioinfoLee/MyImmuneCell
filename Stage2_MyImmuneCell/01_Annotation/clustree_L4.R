library(clustree)
library(cowplot)
args <- commandArgs(T)

data = read.csv(paste0(as.character(args[1]),'/Level4_TOTALVI_ForClustree.csv'),row.names=1)
if (nrow(data)>50000){
    data = data[sample(nrow(data), size = 50000),]
}
p=clustree(data,'L4_leiden_TOTALVI_')
ggsave2(paste0(as.character(args[1]),"/L4_leiden_TOTALVI_clustree.png"),p,width=15,height=10)