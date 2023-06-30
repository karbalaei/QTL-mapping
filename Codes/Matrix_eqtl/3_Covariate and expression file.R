library(SingleCellExperiment)
library(dreamlet)
library(tidyverse)
library(readr)
library(readxl)
library(scran) # for PCA
library(data.table)


##### expression and covariate files #####

setwd("D:/Fernando/Single cell analysis/")

set.seed(100) # See below.

pb_mdd_k3000_10_broad <- readRDS("pb_mdd_k3000_10_broad.rds")

pb_mdd_k3000_10_broad

(cell_type_names <-  assayNames(pb_mdd_k3000_10_broad) %>% str_replace_all("/" , "_") )


rowData_pb <-  rowData(pb_mdd_k3000_10_broad) %>% as.data.frame() %>% dplyr::select(1,2,4)

Single_Cell_sample_id_data <- read_excel("Single_Cell_sample_updated_pheno_wSNPpcs.xlsx", 
                                         sheet = "data") %>%
  dplyr::select(c(1, 2))%>% dplyr::rename(Id = SampleID , SampleId = genoSample.x)



All_cells <-  readRDS("res_processed_dreamlet_k3000_10_broad.rds")

names(All_cells) <- names(All_cells) %>% str_replace_all("/" , "_")  

names(All_cells)


pd <-  colData(pb_mdd_k3000_10_broad) %>% as.data.frame() %>% dplyr::select(c(1 , 8 , 23:27))%>% 
  dplyr::rename(PrimaryDx = group_id) %>% 
  mutate(Id =colnames(pb_mdd_k3000_10_broad) )%>% 
  left_join(Single_Cell_sample_id_data , by = "Id")%>%
  relocate( SampleId) %>% dplyr::select(-9)


#cell <-  "excit_mixed"

for ( i in 1: length(cell_type_names)) {
  
  dir.create(paste0("Tensorqtl/Broad_cell_type/Matrix_eQTL/",cell_type_names[i] , "_data" ))
  
}

setwd("D:/Fernando/Single cell analysis")

i = 1

for ( i in 1: length(cell_type_names)) {
  
  dir <-  paste0("D:/Fernando/Single cell analysis/Tensorqtl/Broad_cell_type/Matrix_eQTL/" , 
                 cell_type_names[i] , "_data" ) 
  setwd(dir)
  
  logcount_data <- All_cells[[cell_type_names[i]]]  %>% as.data.frame()
  
  colnames(logcount_data) <-   colnames(logcount_data) %>% as.data.frame()%>%
    set_names("Id") %>% left_join(Single_Cell_sample_id_data , by = "Id") %>%
    pull(2)
  
  logcount_data%>%
    mutate(geneid = rownames(logcount_data)) %>%
    relocate(geneid) %>%
    fwrite(paste0(cell_type_names[i], "_GE.txt"),  sep = "\t" , quote = F , row.names = F) # write a compressed file read.table(gzfile("test.dat.gz"),row.names=1)# read it back in
  
  ## PCA ##
  
  pd_ready <- pd %>% dplyr::filter(SampleId %in% colnames(logcount_data))
  
  mod <- model.matrix(~ PrimaryDx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5, data = pd_ready)
  
  pca <- prcomp(t(logcount_data))
  
  message(paste("start num sv for" , cell_type_names[i]))
  
  vfilter <-  nrow(logcount_data)
  
  k <- sva::num.sv(logcount_data, mod, vfilter = vfilter)
  
  PCs <- pca$x[, 1:k] %>% as.data.frame()
  
  PCs$SampleId  <-  rownames(PCs)
  
  covariates_Final <-  pd_ready %>% left_join(PCs , by = "SampleId")%>% 
    mutate(Sex = factor(if_else( Sex =="F" , 1 , 2))) %>%
    mutate(PrimaryDx = factor(if_else( PrimaryDx =="Control" , 0 , 1)) )
  
  #rownames(covariates_Final) <-  covariates_Final %>% pull(1)
  
  
  # covariates_Final = setNames(data.frame(t(covariates_Final[,-1])), 
  #                             covariates_Final[,1])
  t(covariates_Final) %>%
    write.table(paste0(cell_type_names[i] ,"_covariates.txt") , row.names = T , col.names = F, sep = "\t" , quote = F)
  
}
