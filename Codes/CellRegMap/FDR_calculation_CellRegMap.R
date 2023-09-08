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
  
  list_csv_data <- do.call("rbind", map(list_csv, readr::read_csv)) %>% mutate(FDR = p.adjust(Pvalues, "fdr"))
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

