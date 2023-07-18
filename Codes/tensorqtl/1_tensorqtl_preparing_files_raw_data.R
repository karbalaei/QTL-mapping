#### Library list ####

library(SingleCellExperiment)
library(dreamlet)
library(tidyverse)
library(readr)
library(readxl)
library(scran) # for PCA
library(data.table)
### Load data ####


setwd("/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/Tensorqtl")

options(scipen = 999) #You can turn off scientific notation for numbers using the option


#### load the TSS file ####

# it already made in python using pyqtl (io.gtf_to_tss_bed) : https://github.com/broadinstitute/pyqtl


TSS_df <- read_csv("/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/Tensorqtl/TSS_df.csv") %>%
  dplyr::rename(ID = gene_id)

#### Single cell data ######

set.seed(100) # See below.

pb_mdd_k3000_10_subtype <- readRDS("pb_mdd_k3000_10_subtype.rds")

pb_mdd_k3000_10_subtype

(cell_type_names <-  assayNames(pb_mdd_k3000_10_subtype) %>% str_replace_all("/" , "_") )


rowData_pb <-  rowData(pb_mdd_k3000_10_subtype) %>% as.data.frame() %>% dplyr::select(1,2,4)

Single_Cell_sample_id_data <- read_xlsx("Single_Cell_sample_updated_pheno_wSNPpcs.xlsx", 
                                         sheet = "data") %>%
  dplyr::select(c(1, 2))%>% dplyr::rename(Id = SampleID , SampleId = genoSample.x)



All_cells <-  readRDS("res_processed_dreamlet_k3000_10_subtype.rds")

names(All_cells) <- names(All_cells) %>% str_replace_all("/" , "_")  

names(All_cells)


pd <-  colData(pb_mdd_k3000_10_subtype) %>% as.data.frame() %>% dplyr::select(c(1 , 8 , 23:27))%>% 
  dplyr::rename(PrimaryDx = group_id) %>% 
  mutate(Id =colnames(pb_mdd_k3000_10_subtype) )%>% 
  left_join(Single_Cell_sample_id_data , by = "Id")%>%
  relocate( SampleId) %>% dplyr::select(-9)


# for ( i in 1: length(cell_type_names)) {
#   
#   dir.create(paste0("Tensorqtl/",cell_type_names[i] , "_results" ))
#   
# }

pca_df <-  data.frame()

#setwd("/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/Tensorqtl/")

i = 1

for ( i in 1: length(cell_type_names)) {
  
  setwd("/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/Tensorqtl/")
  
  pb_mdd_k3000_10_subtype <- readRDS("pb_mdd_k3000_10_subtype.rds")
  
  colnames(pb_mdd_k3000_10_subtype) <-   colnames(pb_mdd_k3000_10_subtype) %>% as.data.frame()%>%
    set_names("Id") %>% left_join(Single_Cell_sample_id_data , by = "Id") %>%
    pull(2)
  
  
  dir <-  paste0("/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/Tensorqtl/" , 
                 cell_type_names[i] , "_results" ) 
  setwd(dir)
  
  Normalized_data <- All_cells[[cell_type_names[i]]]  %>% as.data.frame()
  
  colnames(Normalized_data) <-   colnames(Normalized_data) %>% as.data.frame()%>%
    set_names("Id") %>% left_join(Single_Cell_sample_id_data , by = "Id") %>%
    pull(2)
  
  Sample_names = colnames(Normalized_data)
  
  Gene_names = rownames(Normalized_data)
  
  
  
  counts(pb_mdd_k3000_10_subtype) <- assay(pb_mdd_k3000_10_subtype, cell_type_names[i])
  
  
  pb_mdd_k3000_10_subtype <-  pb_mdd_k3000_10_subtype[ names(pb_mdd_k3000_10_subtype) %in% Gene_names,
                                                                       colnames(pb_mdd_k3000_10_subtype) %in% Sample_names]
  
  
  logcounts(pb_mdd_k3000_10_subtype) <-
    edgeR::cpm(edgeR::calcNormFactors(pb_mdd_k3000_10_subtype),
               log = TRUE,
               prior.count = 1)
  
  
  logcount_data = assay(pb_mdd_k3000_10_subtype, "logcounts") %>% as.data.frame()
  
  
  logcount_data%>%
    merge(rowData_pb, 
          by = 'row.names')%>% 
    left_join(TSS_df , by = "ID" )%>%
    dplyr::select(!(c("CHR", "Symbol" ,  "Row.names")))%>%
    dplyr::rename(gene_id = ID ) %>%
    relocate(c(  chr, start , end , gene_id) )%>%
    dplyr::filter(chr %in% c(1:22 , "X" )) %>% #MT should add to list
    mutate(chr = factor(chr ,  levels = c(1:22 ,  "X") )) %>%
    arrange(chr , start , end) %>%
    dplyr::rename("#chr" = chr) %>% 
    #format(round(8), nsmall = 8) %>%
    write.table(gzfile(paste0(cell_type_names[i], "_edgeR.bed.gz")),  sep = "\t" , quote = F , row.names = F) # write a compressed file read.table(gzfile("test.dat.gz"),row.names=1)# read it back in
  
  ###### PCA #####
  
  pd_ready <- pd %>% dplyr::filter(SampleId %in% colnames(logcount_data))
  
  mod <- model.matrix(~ PrimaryDx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5, data = pd_ready)
  
  pca <- prcomp(t(logcount_data))
  
  message(paste("start num sv for" , cell_type_names[i]))
  
  vfilter <-  nrow(logcount_data)
  
  k <- sva::num.sv(logcount_data, mod, vfilter = vfilter)
  
  pca_df[i , 1 ] <-  cell_type_names[i]
  pca_df[i , 2 ] <-  k
    
  PCs <- pca$x[, 1:k] %>% as.data.frame()
  
  PCs$SampleId  <-  rownames(PCs)
  
  covariates_Final <-  pd_ready %>% left_join(PCs , by = "SampleId")%>% 
    mutate(Sex = factor(if_else( Sex =="F" , 1 , 2))) %>%
    mutate(PrimaryDx = factor(if_else( PrimaryDx =="Control" , 0 , 1)) )
  
  rownames(covariates_Final) <-  covariates_Final %>% pull(1)
  

  covariates_Final = setNames(data.frame(t(covariates_Final[,-1])), 
                              covariates_Final[,1])

  
   write.table(covariates_Final ,  paste0(cell_type_names[i] ,"_edgeR.covariates.txt") , row.names = T , sep = "\t" , quote = F)
 
}





