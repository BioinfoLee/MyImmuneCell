library(scRepertoire)
library(dplyr)
library(data.table)
library('Seurat')
library(ggplot2) 
library(ggpubr)
library(reshape2)
library(dplyr)
library(rlang)
library(scater)
library(tictoc)
library(stringr)
library(RColorBrewer)
library(Startrac)

---------------------------------------
##########TCR Diversity################
---------------------------------------
metadata <- fread("1kproject/metadata/metadata_final.csv")
tcr <- readRDS('Only_gex_tcr.rds')
tcr_df <- bind_rows(tcr,.id='sample')
### calculating sample_level diversity results 
clonalDiversity_all <- clonalDiversity(tcr, 
                cloneCall = "CTstrict",
                chain = "both",   
                group_by = 'sample',                              
                exportTable = TRUE,
                skip.boots = TRUE)
clonalDiversity_all_meta <- dplyr::left_join(clonalDiversity_all, metadata, by = "sample")

### prepare store_list
data_list <- list(
  "all" = clonalDiversity_all_meta 
)

### calculating cell_type level diversity results
L4_celltype <- unique(tcr_df$Classification_L4)
for (celltype in L4_celltype) {
  df <- subset(tcr_df, Classification_L4 == celltype)
  df_list <- split(df, df$sample)
  
  clonalDiversity <- clonalDiversity(df_list, 
                cloneCall = "CTstrict",
                chain = "both",                                 
                exportTable = TRUE,
                skip.boots = TRUE)
  
  clonalDiversity_meta <- dplyr::left_join(clonalDiversity, metadata, by = c("Group" = "sample"))
  data_list[[celltype]] <- clonalDiversity_meta
}

### plotting 
for (name in names(data_list)) {
  df <- data_list[[name]]
  
  p_diversity <- ggplot(df, aes(x = .data[["10_year_intervals"]], 
                               y = .data[["value"]])) +
    geom_boxplot(width = 0.6, color = "black", fill = "white", outlier.shape = NA) +
    geom_jitter(
      aes(fill = `10_year_intervals`),
      width = 0.15,
      size = 3.5,
      shape = 21,
      stroke = 0.6,
      color = "black",
      alpha = 0.9
    ) + 
    geom_smooth(
      aes(x = as.numeric(factor(`10_year_intervals`)), y = value, group = 1),
      method = "loess",
      formula = y ~ x,
      se = TRUE,
      color = "black",
      linewidth = 0.5,
      fill = "grey70",
      alpha = 0.2,
      inherit.aes = FALSE
    ) +
    labs(
      x = "Age Group", 
      y = "Shannon index",
      title = ""
    )  +
    theme_pubr(base_size = 16) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.title.y = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 14, color = "black"),
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3)))
  
  # 保存
  ggsave(paste0('1kproject/plot/diversity/', name, '_Shannon_index.pdf'), 
         p_diversity, height = 5, width = 5)
  fwrite(df, paste0("1kproject/csv/diversity/", name, "_Shannon_index.csv"))
}


---------------------------------------
##########BCR Diversity################
---------------------------------------
bcr <- readRDS('Only_gex_bcr.rds')
bcr_df <- bind_rows(bcr,.id='sample')
### calculating sample_level diversity results 
clonalDiversity_all <- clonalDiversity(bcr, 
                cloneCall = "CTstrict",
                chain = "both",   
                group_by = 'sample',                              
                exportTable = TRUE,
                skip.boots = TRUE)
clonalDiversity_all_meta <- dplyr::left_join(clonalDiversity_all, metadata, by = "sample")

### prepare store_list
data_list <- list(
  "all" = clonalDiversity_all_meta 
)

### calculating cell_type level diversity results
L4_celltype <- unique(bcr_df$Classification_L4)
for (celltype in L4_celltype) {
  df <- subset(bcr_df, Classification_L4 == celltype)
  df_list <- split(df, df$sample)
  
  clonalDiversity <- clonalDiversity(df_list, 
                cloneCall = "CTstrict",
                chain = "both",                                 
                exportTable = TRUE,
                skip.boots = TRUE)
  
  clonalDiversity_meta <- dplyr::left_join(clonalDiversity, metadata, by = c("Group" = "sample"))
  data_list[[celltype]] <- clonalDiversity_meta
}

### plotting 
for (name in names(data_list)) {
  df <- data_list[[name]]
  
  p_diversity <- ggplot(df, aes(x = .data[["10_year_intervals"]], 
                               y = .data[["value"]])) +
    geom_boxplot(width = 0.6, color = "black", fill = "white", outlier.shape = NA) +
    geom_jitter(
      aes(fill = `10_year_intervals`),
      width = 0.15,
      size = 3.5,
      shape = 21,
      stroke = 0.6,
      color = "black",
      alpha = 0.9
    ) + 
    geom_smooth(
      aes(x = as.numeric(factor(`10_year_intervals`)), y = value, group = 1),
      method = "loess",
      formula = y ~ x,
      se = TRUE,
      color = "black",
      linewidth = 0.5,
      fill = "grey70",
      alpha = 0.2,
      inherit.aes = FALSE
    ) +
    labs(
      x = "Age Group", 
      y = "Shannon index",
      title = ""
    )  +
    theme_pubr(base_size = 16) +
    theme(
      legend.position = "none",
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
      axis.title.y = element_text(size = 15, face = "bold"),
      axis.text = element_text(size = 14, color = "black"),
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.3)))
  
  # 保存
  ggsave(paste0('1kproject/plot/diversity/', name, '_Shannon_index.pdf'), 
         p_diversity, height = 5, width = 5)
  fwrite(df, paste0("1kproject/csv/diversity/", name, "_Shannon_index.csv"))
}
---------------------------------------
#######TCR cloneSize########
---------------------------------------
tcr <- readRDS('Only_gex_tcr.rds')
tcr_df <- bind_rows(tcr,.id='sample')
results_list <- vector("list", length(tcr))
cloneCall <- "CTstrict"
cloneSize <- c("Single"= 1,"2-5" = 5,"6-20" = 20,"21-30" = 30,"31-50" = 50,"51-100"=100,"100-200"=200,">200"=Inf)
for (i in seq_along(tcr)) {
        data <- data.frame(tcr[[i]], stringsAsFactors = FALSE)
        data2 <- unique(data[,c("barcode", cloneCall)])
        # 计算克隆比例和频率
        data2 <- data2 %>% 
            group_by(data2[,cloneCall]) %>%
            summarise(clonalProportion = dplyr::n()/nrow(data2), 
                      clonalFrequency = dplyr::n())
            colnames(data2)[1] <- cloneCall
            data <- merge(data, data2, by = cloneCall, all = TRUE)
    
if(max(na.omit(data[,"clonalFrequency"])) > cloneSize[length(cloneSize)]) {
  cloneSize[length(cloneSize)] <- max(na.omit(data[,"clonalFrequency"]))
}
data$cloneSize <- NA
break_points <- c(0, cloneSize)  
labels <- names(cloneSize)
data$cloneSize <- cut(data$clonalFrequency,
                       breaks = break_points,
                       labels = labels,
                       include.lowest = FALSE,
                       right = TRUE)
results_list[[i]] <- data
}
Con.df <- bind_rows(results_list)
df_tcr <- Con.df
saveRDS(df_tcr,"scTCR_Clonesize.rds",compress="bzip2")

colorblind_vector_tcr <- c(
  `>200` =  "#5a2f25",
  `100-200` =  "#884738",
  `51-100` =  "#b55e4a",
  `31-50` = "#e2765d",
  `21-30` = "#e8917d",
  `6-20` = "#eead9e",
  `2-5` = "#f3c8be",
  `Single` = "#f9e4df",
  "NA" = "#B7BBC4"
)
clone_size_order <- c(
  "NA",  
  "Single",
  "2-5",
  "6-20",
  "21-30",
  "31-50",
  "51-100",
  "100-200",
  ">200"
)


#######barplot#######
df_age_tcr <- df_tcr
freq_table_age_tcr <- table(df_age_tcr$Classification_L4, df_age_tcr$cloneSize, useNA = "no")
freq_percent_age_tcr <- as.data.frame.matrix(prop.table(freq_table_age_tcr, margin = 1) )
freq_percent_age_tcr$Classification_L4 <- rownames(freq_percent_age_tcr)
freq_percent_age_tcr <- freq_percent_age_tcr[ , !names(freq_percent_age_tcr) %in% c("NA")]
freq_long_age_tcr <- freq_percent_age_tcr %>%
  melt(id.vars = "Classification_L4", 
       variable.name = "cloneSize", 
       value.name = "abundance") 
freq_long_age_tcr <- freq_long_age_tcr %>%
  mutate(cloneSize = factor(cloneSize, levels = clone_size_order))
pdf("./plot/Clonesize_barplot/Age_clonesize_barplot_TCR.pdf", width = 7.6, height = 5.5)
ggplot(freq_long_age_tcr, aes(x = Classification_L4, y = abundance, fill = cloneSize)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.25) + 
  scale_fill_manual(name = "Clone Size", values = colorblind_vector_tcr) + 
  labs(
    title = "",
    x = "",
    y = "Relative Percentage of Cells" 
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
 # scale_y_break(c(0.2,0.8), scales=0.2, ticklabels=c(0.8,0.1), space=0.2)
dev.off()


## Umap cloneSize
df_tcr <- df_tcr %>%
mutate(
    cloneSize_plot = factor(
      case_when(
        is.na(cloneSize) ~ "NA",
        TRUE ~ as.character(cloneSize)
      ),
      levels = clone_size_order  
    )
  )
p <- ggplot() +
  geom_point(
    data = subset(df_tcr, is.na(cloneSize_plot)),
    aes(x = UMAP_1, y = UMAP_2, color =cloneSize_plot),
    size = 2,
    alpha = 1  ) +
  geom_point(
    data = subset(df_tcr, cloneSize_plot != "NA"),
    aes(x = UMAP_1, y = UMAP_2, color = cloneSize_plot),
    size = 2,
    alpha = 1
  ) +
scale_color_manual(
    name = "Clone Size",
    values = colorblind_vector_tcr,
    breaks = clone_size_order,  
    drop = FALSE  
  ) +
  labs(x = "UMAP 1", y = "UMAP 2") +
 theme_classic() +
  theme(
    axis.line = element_line(color = "black", size = 0.5),
    axis.line.x.top = element_blank(),    
    axis.line.y.right = element_blank(), 
    axis.ticks = element_blank(),         
    axis.text = element_blank(),  
    panel.background = element_rect(fill = 'white'),
    plot.background = element_rect(fill = "white"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = "white"))
      ggsave(
    filename = paste0("./plot/Clonesize_umap/clonesize_Umap_TCR.png"),
    plot = p,
    width = 5.5, 
    height = 5.5,
    #device = "pdf"
  ) 
---------------------------------------
######BCR cloneSize########
---------------------------------------
bcr <- readRDS('Only_gex_bcr.rds')
bcr_df <- bind_rows(bcr,.id='sample')
results_list <- vector("list", length(bcr))
cloneCall <- "CTstrict"
cloneSize <- c("Single"= 1,"2-5" = 5,"6-20" = 20,"21-30" = 30,"31-50" = 50,"51-100"=100,"100-200"=200,">200"=Inf)
for (i in seq_along(bcr)) {
        data <- data.frame(bcr[[i]], stringsAsFactors = FALSE)
        data2 <- unique(data[,c("barcode", cloneCall)])
        # 计算克隆比例和频率
        data2 <- data2 %>% 
            group_by(data2[,cloneCall]) %>%
            summarise(clonalProportion = dplyr::n()/nrow(data2), 
                      clonalFrequency = dplyr::n())
            colnames(data2)[1] <- cloneCall
            data <- merge(data, data2, by = cloneCall, all = TRUE)
    
if(max(na.omit(data[,"clonalFrequency"])) > cloneSize[length(cloneSize)]) {
  cloneSize[length(cloneSize)] <- max(na.omit(data[,"clonalFrequency"]))
}
data$cloneSize <- NA
break_points <- c(0, cloneSize)  
labels <- names(cloneSize)
data$cloneSize <- cut(data$clonalFrequency,
                       breaks = break_points,
                       labels = labels,
                       include.lowest = FALSE,
                       right = TRUE)
results_list[[i]] <- data
}
bcr_merged <- bind_rows(results_list)
df_bcr <- bcr_merged
saveRDS(df_bcr,"scBCR_Clonesize.rds",compress="bzip2")


colorblind_vector_bcr <- c(
  `>200` =  "#8b00fc",
  `100-200` =  "#0008fc",
  `51-100` =  "#0061fc",
  `31-50` = "#008ffc",
  `21-30` = "#00bdfc",
  `6-20` = "#00e7fc",
  `2-5` = "#00fcf0",
  `Single` = "#007056",
  "NA" = "#B7BBC4"
)
clone_size_order <- c(
  "NA",  
  "Single",
  "2-5",
  "6-20",
  "21-30",
  "31-50",
  "51-100",
  "100-200",
  ">200"
)


#######barplot
df_age_bcr <- df_bcr
freq_table_age_bcr <- table(df_age_bcr$Classification_L4, df_age_bcr$cloneSize, useNA = "no")
freq_percent_age_bcr <- as.data.frame.matrix(prop.table(freq_table_age_bcr, margin = 1) )
freq_percent_age_bcr$Classification_L4 <- rownames(freq_percent_age_bcr)
freq_percent_age_bcr <- freq_percent_age_bcr[ , !names(freq_percent_age_bcr) %in% c("NA")]
freq_long_age_bcr <- freq_percent_age_bcr %>%
  melt(id.vars = "Classification_L4", 
       variable.name = "cloneSize", 
       value.name = "abundance") 
freq_long_age_bcr <- freq_long_age_bcr %>%
  mutate(cloneSize = factor(cloneSize, levels = clone_size_order))
pdf("1kproject/plot/Clonesize_barplot/Age_clonesize_barplot_BCR.pdf", width = 7.6, height = 5.5)
ggplot(freq_long_age_bcr, aes(x = Classification_L4, y = abundance, fill = cloneSize)) +
  geom_bar(stat = "identity", position = "fill", color = "black", linewidth = 0.25) + 
  scale_fill_manual(name = "Clone Size", values = colorblind_vector_bcr) + 
  labs(
    title = "",
    x = "",
    y = "Relative Percentage of Cells" 
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
 # scale_y_break(c(0.2,0.8), scales=0.2, ticklabels=c(0.8,0.1), space=0.2)
dev.off()


## Umap cloneSize
df_bcr <- df_bcr %>%
mutate(
    cloneSize_plot = factor(
      case_when(
        is.na(cloneSize) ~ "NA",
        TRUE ~ as.character(cloneSize)
      ),
      levels = clone_size_order  
    )
  )
p <- ggplot() +
  geom_point(
    data = subset(df_bcr, is.na(cloneSize_plot)),
    aes(x = UMAP_1, y = UMAP_2, color =cloneSize_plot),
    size = 2,
    alpha = 1  ) +
  geom_point(
    data = subset(df_bcr, cloneSize_plot != "NA"),
    aes(x = UMAP_1, y = UMAP_2, color = cloneSize_plot),
    size = 2,
    alpha = 1
  ) +
scale_color_manual(
    name = "Clone Size",
    values = colorblind_vector_bcr,
    breaks = clone_size_order,  
    drop = FALSE  
  ) +
  labs(x = "UMAP 1", y = "UMAP 2") +
 theme_classic() +
  theme(
    axis.line = element_line(color = "black", size = 0.5),
    axis.line.x.top = element_blank(),    
    axis.line.y.right = element_blank(), 
    axis.ticks = element_blank(),         
    axis.text = element_blank(),  
    panel.background = element_rect(fill = 'white'),
    plot.background = element_rect(fill = "white"),
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 14),
    legend.text = element_text(size = 12),
    legend.key = element_rect(fill = "white"))
      ggsave(
    filename = paste0("./plot/Clonesize_umap/clonesize_Umap_BCR.png"),
    plot = p,
    width = 5.5, 
    height = 5.5,
    #device = "pdf"
  ) 


##########STARTRAC anlaysis function############
calculate_entropy <- function(x) {
  p <- x / sum(x)
  p <- p[p > 0]  
  -sum(p * log2(p))
}

calculate_sample_expa <- function(df) {
  patients <- unique(df$patient)
  results <- list()
  
  for (p in patients) {
    sample_data <- df %>% filter(patient == p)
    
    sample_matrix <- sample_data %>%
      group_by(clone.id) %>%
      summarise(count = n(), .groups = 'drop') %>%
      as.data.frame()
    
    rownames(sample_matrix) <- sample_matrix$clone.id
    sample_matrix$clone.id <- NULL
    sample_mat <- as.matrix(sample_matrix)
    
    valid_sample <- colSums(sample_mat) > 0
    sample_mat_valid <- sample_mat[, valid_sample, drop = FALSE]
    
    if (ncol(sample_mat_valid) == 0) {
      next
    }
    
    entropy_values <- apply(sample_mat_valid, 2, calculate_entropy)
    max_entropy_values <- log2(colSums(sample_mat_valid > 0))
    expa_values <- 1 - entropy_values / max_entropy_values
    
    # 创建结果dataframe
    result_df <- data.frame(
      patient = p,
      expa = expa_values,
      NCells = colSums(sample_mat_valid),
      stringsAsFactors = FALSE
    )
    
    results[[p]] <- result_df
  }
  
  # 合并所有样本的结果
  final_result <- do.call(rbind, results)
  return(final_result)
}
--------------------------------------------
####BCR  STARTRAC analysis###################
--------------------------------------------
bcr_dataframe <- readRDS("scBCR_Clonesize.rds")
subB.meta <- as.data.table(bcr_dataframe)

## convert data
setnames(subB.meta,
         old = c("barcode", "CTstrict", "sample", "Classification_L4", 
                "Classification_L3", "Classification_L2", "Classification_L1", 
                "10_year_intervals"),
         new = c("Cell_Name", "clone.id", "patient", "majorCluster", 
                 "majorCluster_L3", "majorCluster_L2", "majorCluster_L1", 
                 "age_group"))
subB.meta[, clone.status := ifelse(cloneSize == "single", "NoClonal", "Clonal")]
subB.meta[, loc := "pbmc"]

final_columns <- c("Cell_Name", "clone.id", "patient", "clone.status",
                   "majorCluster", "majorCluster_L3", "majorCluster_L2", 
                   "majorCluster_L1", "loc","age_group")
subB.meta <- subB.meta[, ..final_columns] 


expa_results <- calculate_sample_expa(subB.meta)
expa_results_meta <- dplyr::left_join(expa_results,metadata,by=c('patient'="sample"))
fwrite(expa_results_meta,'./csv/BCR_sample_level_STARTRAC_expansion_score.csv')
## cell_type level STARTRAC index
tic("Startrac.run")
suppressWarnings(out <- Startrac.run(subB.meta, proj="1k",verbose=F))
cluster_dynamics <- out@cluster.data
data_Startrac <- cluster_dynamics %>% filter(aid != "1k")
metadata  <- fread("1kproject/metadata/metadata_final.csv")
data_Startrac <- dplyr::left_join(data_Startrac,metadata,by=c("aid"="sample"))
fwrite(data_Startrac,'1kproject/csv/BCR_sample_celltype_STARTRAC_score.csv')

b_cell_development_order <- c(
  'Naïve B cells',
  'Atypical naïve B cells',
  'Early memory B cells',
  'Non-switched memory B cells',
  'Switched Memory B cells',
  'CD95 memory B cells',
  'CD27+ IgD+ atypical memory B cells',
  'CD27+ IgD- atypical memory B cells', 
  'CD27- IgD- atypical memory B cells',
  'Plamsablasts',
  'IGKChi Plasma cells',
  'IGLL5hi Plasma cells'
)
data_Startrac <- data_Startrac %>%
          mutate(majorCluster= factor(majorCluster, levels = b_cell_development_order))

colors <- c("#305DA9","#0A8BCD","#4BB3D8","#A2D8E9","#DDEEF3","#FEE3B3","#FFB37D","#F27A5A","#CE3C4D","#B51A44","#800000","#F0C040")
L4_celltype <- unique(data_Startrac$majorCluster)
for (celltype in L4_celltype) {
  df <- subset(data_Startrac, majorCluster == celltype)
  p <- ggplot(df, aes(x = `10_year_intervals`, y = expa, fill = `10_year_intervals`)) +
    geom_boxplot(
      color = "black",
      width = 0.6,
      outlier.shape = NA
    ) +
    geom_jitter(
      width = 0.15,
      size = 1,
      color = "black"
    ) +
    geom_smooth( 
      aes(x = as.numeric(factor(`10_year_intervals`)), y = expa, group = 1), 
      method = "loess", 
      formula = y ~ x, 
      se = TRUE, 
      color = "black", 
      linewidth = 0.5,
      fill = "grey70", 
      alpha = 0.2, 
      inherit.aes = FALSE 
    ) +
    scale_fill_manual(values = colors) +
    labs(
      x = "",
      y = "STARTRAC expansion score",
      title = celltype
    ) +
    theme_pubr() +
    theme(
      plot.title=element_text(hjust=0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      panel.background = element_rect(fill = "white"),
      legend.position = "none"
    ) +
    coord_cartesian(ylim = c(0, 0.8))
  ggsave(paste0("1kproject/plot/STARTRAC-expa/", celltype, "_age_STARTRAC-expa.pdf"), plot = p, width = 7.5, height = 5.5)
}

p <- ggplot(expa_results_meta, aes(x = `10_year_intervals`, y = expa,fill = `10_year_intervals`)) +
  geom_boxplot(      
    color = "black",            
    width = 0.6,               
    outlier.shape = NA        
  ) +
  geom_jitter(
    width = 0.15,
    size = 1,      
    color = "black"    
  )  +
    geom_smooth(
      aes(x = as.numeric(factor(`10_year_intervals`)), y = expa, group = 1), 
      method = "loess", 
      formula = y ~ x,
      se = TRUE, 
      color = "black",
      linewidth = 0.5,
      fill = "grey70", 
      alpha = 0.2, 
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = colors) +
  labs(        
    x = "",                     
    y = "STARTRAC expansion score",
    title = "Total B-cell subsets"
  ) +
  theme_pubr() +               
  theme(
    plot.title=element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),    
    axis.text.y = element_text(size = 10),                             
    panel.background = element_rect(fill = "white"),
    legend.position = "none"
  ) #+ coord_cartesian(ylim = c(0, 0.8)) 

ggsave("./plot/STARTRAC-expa/sample_level_STARTRAC-expa.pdf",
           plot=p,
           width=7.5,
           height=5.5)

--------------------------------------------
#####TCR STARTRAC analysis###################
--------------------------------------------
TCR_dataframe <- readRDS("scTCR_Clonesize.rds")
subT.meta <- as.data.table(TCR_dataframe)
setnames(subT.meta,
         old = c("barcode", "CTstrict", "sample.x", "Classification_L4", 
                "Classification_L3", "Classification_L2", "Classification_L1", 
                "10_year_intervals"),
         new = c("Cell_Name", "clone.id", "patient", "majorCluster", 
                 "majorCluster_L3", "majorCluster_L2", "majorCluster_L1", 
                 "age_group"))
subT.meta[, clone.status := ifelse(cloneSize == "single", "NoClonal", "Clonal")]
subT.meta[, loc := "pbmc"]
final_columns <- c("Cell_Name", "clone.id", "patient", "clone.status",
                   "majorCluster", "majorCluster_L3", "majorCluster_L2", 
                   "majorCluster_L1", "loc","age_group")
subT.meta <- subT.meta[, ..final_columns] 
tic("calculate_sample_expa")
expa_results <- calculate_sample_expa(subT.meta)
metadata  <- fread("1kproject/metadata/metadata_final.csv")
expa_results_meta <- dplyr::left_join(expa_results,metadata,by=c('patient'="sample"))
fwrite(expa_results,"1kproject/csv/TCR_sample_level_STARTRAC-expa.csv")
tic("Startrac.run")
suppressWarnings(out <- Startrac.run(subT.meta, proj="1k",verbose=F))
cluster_dynamics <- out@cluster.data
data_Startrac_tcr <- cluster_dynamics %>% filter(aid != "1k")
data_Startrac_tcr <- dplyr::left_join(data_Startrac_tcr,metadata,by=c("aid"="sample"))
fwrite(data_Startrac_tcr,'1kproject/csv/TCR_celltype_STARTRAC_score.csv')

colors <- c("#305DA9","#0A8BCD","#4BB3D8","#A2D8E9","#DDEEF3","#FEE3B3","#FFB37D","#F27A5A","#CE3C4D","#B51A44")
L4_celltype <- unique(data_Startrac_tcr$majorCluster)
for (celltype in L4_celltype ){
  df <- subset(data_Startrac_tcr,majorCluster== celltype)
p <- ggplot(df, aes(x = `10_year_intervals`, y = expa,fill = `10_year_intervals`)) +
  geom_boxplot(      
    color = "black",            
    width = 0.6,               
    outlier.shape = NA        
  ) +
  geom_jitter(
    width = 0.15,
    size = 1,      
    color = "black"    
  )  +
     geom_smooth(
      aes(x = as.numeric(factor(`10_year_intervals`)), y = expa, group = 1), 
      method = "loess", 
      formula = y ~ x,
      se = TRUE, 
      color = "black",
      linewidth = 0.5,
      fill = "grey70", 
      alpha = 0.2, 
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = colors) +
  labs(        
    x = "",                     
    y = "STARTRAC expansion score",
    title = celltype
  ) +
  theme_pubr() +               
  theme(
    plot.title=element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),    
    axis.text.y = element_text(size = 10),                             
    panel.background = element_rect(fill = "white"),
    legend.position = "none"
  ) + coord_cartesian(ylim = c(0, 0.8)) 
    ggsave(paste0("./plot/STARTRAC expansion score/celltype/",celltype,"_age_STARTRAC-expa.pdf"),
           plot=p,
           width=7.5,
           height=5.5)
}

colors <- c("#305DA9","#0A8BCD","#4BB3D8","#A2D8E9","#DDEEF3","#FEE3B3","#FFB37D","#F27A5A","#CE3C4D","#B51A44")
p <- ggplot(expa_results_meta, aes(x = `10_year_intervals`, y = expa,fill = `10_year_intervals`)) +
  geom_boxplot(      
    color = "black",            
    width = 0.6,               
    outlier.shape = NA        
  ) +
  geom_jitter(
    width = 0.15,
    size = 1,      
    color = "black"    
  )  +
   geom_smooth(
      aes(x = as.numeric(factor(`10_year_intervals`)), y = expa, group = 1), 
      method = "loess", 
      formula = y ~ x,
      se = TRUE, 
      color = "black",
      linewidth = 0.5,
      fill = "grey70", 
      alpha = 0.2, 
      inherit.aes = FALSE
    ) +
    scale_fill_manual(values = colors) +
  labs(        
    x = "",                     
    y = "STARTRAC expansion score",
    title = "Total αβ T cell subsets"
  ) +
  theme_pubr() +               
  theme(
    plot.title=element_text(hjust=0.5),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),    
    axis.text.y = element_text(size = 10),                             
    panel.background = element_rect(fill = "white"),
    legend.position = "none"
  ) 
    
ggsave(paste0("./plot/STARTRAC expansion score/sample_level_STARTRAC-expa.pdf"),
           plot=p,
           width=7.5,
           height=5.5)