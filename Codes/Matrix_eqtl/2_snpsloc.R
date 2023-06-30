library(tidyverse)
library(readr)

##### SNPlocation file #####

MDD_control_sc_maf <-  read.table("MDD_control_sc_maf1.txt", row.names = 1, sep = "\t", header = T)

nams <-  rownames(MDD_control_sc_maf) %>% as.data.frame() %>% 
  setNames("snpid") %>%separate(col = snpid, into =c("chr" ,"pos",  "Allele","Alt_Allele"), 
                                sep = ":" , remove = F , extra = "merge") %>%
  mutate(pos = as.numeric(pos)) 

nams %>% dplyr::select(!(c(4, 5))) %>%
  write.table( "snpsloc.txt" , row.names = F , sep = "\t" , quote = F)
