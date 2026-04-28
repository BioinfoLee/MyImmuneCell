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
library(stringr)

extract_CTstrict <- function(text) {
  ighv <- str_extract(text$IGH, "IGHV[^.]*")
  ighj <- str_extract(text$IGH, "IGHJ[^.]*")
  igkv <- str_extract(text$IGK, "IGKV[^.]*")
  igkj <- str_extract(text$IGK, "IGKJ[^.]*")
  iglv <- str_extract(text$IGL, "IGLV[^.]*")
  iglj <- str_extract(text$IGL, "IGLJ[^.]*")
  
  light_chain_v <- ifelse(!is.na(igkv) & igkv != "", igkv, iglv)
  light_chain_j <- ifelse(!is.na(igkj) & igkj != "", igkj, iglj)
   
  text$CTstrict<-paste0(
    ighv, ".", ighj, ";",
    text$cdr3_nt1, "_",
    light_chain_v, ".", light_chain_j, ";",
    text$cdr3_nt2
  )
    return(text)
}

---------------------------------
###### TCR preprocessing ########
---------------------------------
##读取数据 
csv_files <- list.files(
  path = "/media/AnalysisDisk2/Lifupeng/Lifupeng/1k_project/TCR_data/data/", 
  pattern = "\\.csv$",      
  full.names = TRUE)
data_list <- lapply(csv_files, fread, header = TRUE, stringsAsFactors = FALSE,fill=TRUE)

names(data_list) <- sapply(csv_files, function(path) { 
  sub("_TCR\\.csv$", "", basename(path))  
}) 
## 去除单链与保留多链中TCR表达量最高的一条配对链
combined.TCR <- combineTCR(data_list, 
                           samples =  names(data_list),
                           removeNA = TRUE, 
                           removeMulti = FALSE, 
                           filterMulti = TRUE)

tcr_df <- bind_rows(combined.TCR,.id="sample")
tcr_df <- tcr_df[!grepl("OR", tcr_df$CTgene), ]
T_rna_celltype_information <- fread('1kproject/cell_type_annotation/T_cell_annotation.csv')
metadata <- fread("1kproject/1k_matadata.csv")
GEX_TCR_all  <- dplyr::left_join(T_rna_celltype_information,tcr_df,by='barcode')
GEX_TCR_all  <- dplyr::left_join(GEX_TCR_all,metadata,by='sample')
gex_tcr_all <- split(GEX_TCR_all,GEX_TCR_all$sample)
saveRDS(gex_tcr_all,'All_gex_tcr.rds')
gex_tcr_only <- GEX_TCR_all %>% filter(!is.na(CTgene)) %>% split(.,.$sample)
saveRDS(gex_tcr_only,'Only_gex_tcr.rds')

---------------------------------
###### BCR preprocessing ########
---------------------------------
csv_files <- list.files(
  path = "/media/AnalysisDisk2/Lifupeng/Lifupeng/1k_project/BCR_data/data/", 
  pattern = "\\.csv$",      
  full.names = TRUE)
data_list <- lapply(csv_files, fread, header = TRUE, stringsAsFactors = FALSE,fill=TRUE)
names(data_list) <- sapply(csv_files, function(path) { 
  sub("_BCR\\.csv$", "", basename(path))  
})

## 去除单链与保留多链中TCR表达量最高的一条配对链
combined.BCR <- combineBCR(data_list, 
                           samples = names(data_list),
                           call.related.clones = FALSE,
                           sequence = "nt",
                           removeNA = TRUE, 
                           removeMulti = FALSE,
                           filterMulti = TRUE)

combined.BCR <- lapply(combined.BCR,function(x){
    extract_CTstrict(x)
})
bcr_df <- bind_rows(combined.BCR,.id='sample')
##去除假基因
bcr_df <- bcr_df[!grepl("OR", bcr_df$CTgene), ]
metadata <- fread("1kproject/1k_matadata.csv")
B_rna_celltype_information <- fread('1kproject/cell_type_annotation/B_cell_annotation.csv')
GEX_BCR_all  <- dplyr::left_join(B_rna_celltype_information,bcr_df,by='barcode')
GEX_BCR_all  <- dplyr::left_join(GEX_BCR_all,metadata,by='sample')
gex_bcr_all <- split(GEX_BCR_all,GEX_BCR_all$sample)
saveRDS(gex_bcr_all,'All_gex_bcr.rds')
gex_bcr_only <- GEX_BCR_all %>% filter(!is.na(CTgene)) %>% split(.,.$sample)
saveRDS(gex_bcr_only,'Only_gex_bcr.rds')
