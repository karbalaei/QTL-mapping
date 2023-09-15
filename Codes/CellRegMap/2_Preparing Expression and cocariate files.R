
#### Library list ####

library(SingleCellExperiment)

library(pak) # to install package
#pkg_install("zenith")
#pkg_install("bedtorch")

# Install package and dependencies

#devtools::install_github("GabrielHoffman/dreamlet")
library(dreamlet)

library(tidyverse)
library(readr)
library(readxl)
library(scran) # for PCA
library(data.table)

### CellRegMap ####



##### load GTF file ####
setwd("D:/Fernando/Single cell analysis/")


genes <- read_delim("genes.gtf", delim = "\t", 
                    escape_double = FALSE, col_names = FALSE, 
                    comment = "#", trim_ws = TRUE) %>%
  setNames(c("CHR" , "Version" , "Type" , "Start" , "End" , "NA" , "Strand" , 
             "NA2" , "Info")) %>%dplyr::filter(Type =="gene") %>%
  separate("Info" , sep = ";" , into = c("ID" , "gene_version" , "Symbol" , "gene_source" , "gene_biotype" )) %>%
  mutate_at("ID", str_remove, "gene_id " ) %>%
  mutate_at("gene_version", str_remove, "gene_version " ) %>%
  mutate_at("Symbol", str_remove, "gene_name " ) %>%
  mutate_at("gene_source", str_remove, "gene_source " ) %>%
  mutate_at("gene_biotype", str_remove, "gene_biotype " )%>%
  dplyr::select((c(1 , 4, 5, 9 , 11)))

genes <-  as.data.frame(sapply(genes, function(x) gsub("\"", "", x))) # remove quotes

genes$Start  <-   genes$Start %>% as.numeric() #start and end supposed to be a numeric
genes$End  <-   genes$End %>% as.numeric()

options(scipen = 999) #You can turn off scientific notation for numbers using the option

write.table(genes , "annotation_table.txt" , sep = "\t" , row.names = F , quote = F)


##### Count data ####

setwd("D:/Fernando/Single cell analysis/")

Single_Cell_sample_id_data <- read_excel("Single_Cell_sample_updated_pheno_wSNPpcs.xlsx", 
                                         sheet = "data") %>%
  dplyr::select(c(1, 2))%>% dplyr::rename(Id = SampleID , SampleId = genoSample.x)

pb_mdd_k3000_10_broad <- readRDS("pb_mdd_k3000_10_broad.rds")

pd <-  colData(pb_mdd_k3000_10_broad) %>% as.data.frame() %>% dplyr::select(c(1 , 8 , 23:27))%>% 
  dplyr::rename(PrimaryDx = group_id) %>% 
  mutate(Id =colnames(pb_mdd_k3000_10_broad) )%>% 
  left_join(Single_Cell_sample_id_data , by = "Id")%>%
  relocate( SampleId) %>% dplyr::select(-9)

colnames(pb_mdd_k3000_10_broad) <-   colnames(pb_mdd_k3000_10_broad) %>% as.data.frame()%>%
  set_names("Id") %>% left_join(Single_Cell_sample_id_data , by = "Id") %>%
  pull(2)


All_cells <-  readRDS("res_processed_dreamlet_k3000_10_broad.rds")

cell_type_names <-  names(All_cells)

i = 1

for ( i in 1 : length(cell_type_names)) {
  
  setwd("D:/Fernando/Single cell analysis/")
  
  
  Normalized_data <- All_cells[[cell_type_names[i]]]  %>% as.data.frame()
  
  colnames(Normalized_data) <-   colnames(Normalized_data) %>% as.data.frame()%>%
    set_names("Id") %>% left_join(Single_Cell_sample_id_data , by = "Id") %>%
    pull(2)
  
  Sample_names = colnames(Normalized_data)
  
  Gene_names = rownames(Normalized_data)
  
  pb_mdd_k3000_10_broad <- readRDS("pb_mdd_k3000_10_broad.rds")
  
  counts(pb_mdd_k3000_10_broad) <- assay(pb_mdd_k3000_10_broad, cell_type_names[i] )
  
  colnames(pb_mdd_k3000_10_broad) <-   colnames(pb_mdd_k3000_10_broad) %>% as.data.frame()%>%
    set_names("Id") %>% left_join(Single_Cell_sample_id_data , by = "Id") %>%
    pull(2)
  
  pb_mdd_k3000_10_broad <-  pb_mdd_k3000_10_broad[ names(pb_mdd_k3000_10_broad) %in% Gene_names,
                                                   colnames(pb_mdd_k3000_10_broad) %in% Sample_names]
  
  count_data = assay(pb_mdd_k3000_10_broad, "counts") %>% as.matrix()  %>% as.data.frame()
  
  GE_data <-  rowData(pb_mdd_k3000_10_broad) %>% as.data.frame() %>% dplyr::select(c(1, 2)) %>%
    merge(count_data, 
          by = 'row.names') %>% select(!(c(1, 3)))
  
  
  logcounts(pb_mdd_k3000_10_broad) <-
    edgeR::cpm(edgeR::calcNormFactors(pb_mdd_k3000_10_broad),
               log = TRUE,
               prior.count = 1)
  
  logcount_data = assay(pb_mdd_k3000_10_broad, "logcounts") %>% as.data.frame()
  
  
  
  annotaton_data <-  genes %>% dplyr::filter(ID %in% GE_data$ID) 
  
  setwd("D:/Fernando/Single cell analysis/Tensorqtl/CellRegMap")
  
  colnames(annotaton_data) <-  c("chromosome" , "start" , "end" , "feature_id" , "Symbol")
  
  write.table(annotaton_data , paste0(cell_type_names[i] ,"_annotation.txt") , sep = "\t" , row.names = F , quote = F)
  
  
  GE_data[GE_data == 0] = 1
  
  write.table(GE_data , paste0(cell_type_names[i] ,"_GE.txt") , sep = "\t" , row.names = F , quote = F)
  
  
  #### sample info file ####
  
  sample_mapping_file <- read_csv("sample_mapping_file.csv") %>% 
    dplyr::filter(genotype_individual_id %in% Sample_names)
  
  
  write.table(sample_mapping_file , paste0(cell_type_names[i] ,"_sample_mapping_file.csv") , sep = "," ,
              row.names = F , quote = F)
  
  ##### PCA ######
  
  pd_ready <- pd %>% dplyr::filter(SampleId %in% colnames(logcount_data))
  
  mod <- model.matrix(~ PrimaryDx + Sex + snpPC1 + snpPC2 + snpPC3 + snpPC4 + snpPC5, data = pd_ready)
  
  pca <- prcomp(t(logcount_data))
  
  message(paste("start num sv for" , cell_type_names[i]))
  
  vfilter <-  nrow(logcount_data)
  
  k <- sva::num.sv(logcount_data, mod, vfilter = vfilter)
  

  
  PCs <- pca$x[, 1:k] %>% as.data.frame()
  
  covariates_Final <-  pd_ready %>%
    mutate(Sex = factor(if_else( Sex =="F" , 1 , 2))) %>%
    mutate(PrimaryDx = factor(if_else( PrimaryDx =="Control" , 1 , 2)) )
  
  rownames(covariates_Final) <-  covariates_Final %>% pull(1)
  
  
  write.table(covariates_Final[,  c(2,3)] ,  paste0(cell_type_names[i] ,"_CellRegMap_covariate.txt") , row.names = T , sep = "\t" , quote = F)
  
  write.table(PCs ,  paste0(cell_type_names[i] ,"_CellRegMap_PCAs.txt") , row.names = T , sep = "\t" , quote = F)
  
  Test = matrix(1, nrow= nrow(PCs) , ncol = 3)
  
  colnames(Test) = paste0("PCA" , c(1:3)) 
  
  row.names(Test) = rownames(PCs)
  
  write.table(Test ,  paste0(cell_type_names[i] ,"_CellRegMap_PCAs_3.txt") , row.names = T , sep = "\t" , quote = F)
  
}


