library(dplyr)
library(purrr)
library(data.table)


cell_type_names <- c("astro","endo","excit", "inhib" ,"micro" ,"oligo" ,"opc")


dir <-  ("/data/fgoes1/RK/CellRegMap/")

i = 1
for (  i in 1 : length(cell_type_names)) {
  
  setwd(paste0(dir , cell_type_names[i]))
  
  print(cell_type_names[i])
  list_csv <- grep(list.files(), pattern= "*_Association.csv", invert=F, value=TRUE)  
  

 # %>% mutate(FDR = p.adjust(Pvalues, "fdr"))
#list_csv_files <- list.files(.)
  #df2 <- readr::read_csv(list_csv_data, id = "file_name")
#df2

#df <-
 # list.files(path = "/data/fgoes1/RK/CellRegMap/endo", pattern = "*_Association.csv") %>% 
  #map_df(~fread(.))
   setwd(dir)
   
 fwrite(list_csv_data , paste0(cell_type_names[i] , "_results_all.txt") , sep = "\t" , quote = F , row.names = F)

   list_csv_data %>% dplyr::filter(FDR < 0.1) %>%
   fwrite(paste0(cell_type_names[i] , "_results_FDR_0.1.txt") , sep = "\t" , quote = F , row.names = F)


  list_csv_data %>% dplyr::filter(FDR < 0.05) %>%
     fwrite(paste0(cell_type_names[i] , "_results_FDR_0.05.txt") , sep = "\t" , quote = F , row.names = F)
	 
	png(paste0(cell_type_names[i] ,"_hist_FDR.png" ), height=1200 , width=1800)
    
    hist(list_csv_data$FDR , breaks= 300)
    
    dev.off()
	
	png(paste0(cell_type_names[i] ,"_hist_Pvalues.png" ), height=1200 , width=1800)
    
    hist(list_csv_data$Pvalues , breaks= 300)
    
    dev.off()
 
	 
}


   fwrite(list_base_FDR , paste0(cell_type_names[i] , "_base_results_all.txt") , sep = "\t" , quote = F , row.names = F)

  
    list_base_FDR %>% dplyr::filter(FDR < 0.1 & af > 0.05 ) %>% fwrite(paste0(cell_type_names[i] , "_base_results_FDR_0.1.txt") , sep = "\t" , quote = F , row.names = F)

    list_base_FDR %>% dplyr::filter(FDR < 0.05  & af > 0.05 ) %>% fwrite(paste0(cell_type_names[i] , "_base_results_FDR_0.05.txt") , sep = "\t" , quote = F , row.names = F)

     list_base_FDR %>% dplyr::filter(FDR < 0.1 & af > 0.01 ) %>% fwrite(paste0(cell_type_names[i] , "_base_results_FDR_0.1_af_0.01.txt") , sep = "\t" , quote = F , row.names = F)

    list_base_FDR %>% dplyr::filter(FDR < 0.05  & af > 0.01 ) %>% fwrite(paste0(cell_type_names[i] , "_base_results_FDR_0.05_af_0.01.txt") , sep = "\t" , quote = F , row.names = F)


       
    png("hist_list_base_FDR.png" , height=1200 , width=1800)
    
    hist(list_base_FDR$FDR , breaks= 300)
    
    dev.off()


    png("hist_list_base_FDR_af.png" , height=1200 , width=1800)
    
    hist(list_base_FDR_af$FDR , breaks= 300)
    
    dev.off()

     
   png("Scatter_list_base_FDR.png" , height=1200 , width=1800)
    
    plot(-log10(list_base_FDR$pval_nominal) ,  list_base_FDR$af)
    
    dev.off()


     png("Scatter_list_base_FDR_af.png" , height=1200 , width=1800)
    
    plot(-log10(list_base_FDR_af$pval_nominal) ,  list_base_FDR_af$af)
    
    dev.off()




}

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()