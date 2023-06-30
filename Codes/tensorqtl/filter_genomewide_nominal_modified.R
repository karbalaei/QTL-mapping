library(tidyverse)
library(jaffelab)
library(miniparquet)
library(sessioninfo)
library(here)
library(data.table)


cell_type_names <- c("astro","endo","excit_L2_3_A", "excit_L4_5_A",  "excit_L4_5_B", "excit_L5_6_A", "excit_L5_6_B" ,
                     "excit_L5_6_C","excit_L6_A" , "excit_L6_B" ,"excit_mixed","inhib_Lamp5_A", "inhib_Lamp5_B",
                     "inhib_PVALB_A","inhib_PVALB_B","inhib_SST", "inhib_VIP" ,"micro_A" , "micro_B" ,"oligo" ,"opc"   )

dir <-  ("/dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/Tensorqtl/")


for (  i in 1 : length(cell_type_names)) {
  
  setwd(paste0(dir , cell_type_names[i] , "_results"))
  
  list_modified <- grep(list.files(), pattern= paste0(cell_type_names[i], "_modified.cis_qtl_pairs.*"), invert=F, value=TRUE)  
  
list_modified_data <- do.call("rbind", map(list_modified, parquet_read))

  list_modified_FDR_af <- list_modified_data %>% dplyr::filter(af > 0.05) %>%
    mutate(FDR = p.adjust(pval_nominal, "fdr")) 
  
  
  fwrite(list_modified_FDR_af , paste0(cell_type_names[i] , "_results_excluding_af_all_modified.txt") , sep = "\t" , quote = F , row.names = F)

   list_modified_FDR_af %>% dplyr::filter(FDR < 0.1) %>%
    fwrite(paste0(cell_type_names[i] , "_results_excluding_af_0.1_modified.txt") , sep = "\t" , quote = F , row.names = F)


  list_modified_FDR <- list_modified_data %>% 
    mutate(FDR = p.adjust(pval_nominal, "fdr")) 
  
   fwrite(list_modified_FDR , paste0(cell_type_names[i] , "_results_all_modified.txt") , sep = "\t" , quote = F , row.names = F)

  
    list_modified_FDR %>% dplyr::filter(FDR < 0.1  & af > 0.05) %>% 
	fwrite(paste0(cell_type_names[i] , "_results_FDR_0.1_modified.txt") , sep = "\t" , quote = F , row.names = F)

    list_modified_FDR %>% dplyr::filter(FDR < 0.05  & af > 0.05) %>%
	 fwrite(paste0(cell_type_names[i] , "_results_FDR_0.05_modified.txt") , sep = "\t" , quote = F , row.names = F)

       list_modified_FDR %>% dplyr::filter(FDR < 0.1  & af > 0.01) %>% 
	fwrite(paste0(cell_type_names[i] , "_results_FDR_0.1_af_0.01_modified.txt") , sep = "\t" , quote = F , row.names = F)

    list_modified_FDR %>% dplyr::filter(FDR < 0.05  & af > 0.01) %>%
	 fwrite(paste0(cell_type_names[i] , "_results_FDR_0.05_af_0.01_modified.txt") , sep = "\t" , quote = F , row.names = F)

       
    png("hist_list_FDR_modified.png" , height=1200 , width=1800)
    
    hist(list_modified_FDR$FDR , breaks= 300)
    
    dev.off()


    png("hist_list_FDR_af_modified.png" , height=1200 , width=1800)
    
    hist(list_modified_FDR_af$FDR , breaks= 300)
    
    dev.off()

     
   png("Scatter_list_FDR_modified.png" , height=1200 , width=1800)
    
    plot(-log10(list_modified_FDR$pval_nominal) ,  list_modified_FDR$af)
    
    dev.off()


     png("Scatter_list_FDR_modified.png" , height=1200 , width=1800)
    
    plot(-log10(list_modified_FDR_af$pval_nominal) ,  list_modified_FDR_af$af)
    
    dev.off()




}





## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()