library(tidyverse)
library(readr)

##### genepos file #####
#### using GTF file 


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

genes[14671,  3 ] <-  23000000

options(scipen = 999) #You can turn off scientific notation for numbers using the option

setwd("D:/Fernando/Single cell analysis/Tensorqtl/MDD_eQTL/")

genes %>% dplyr::select(-4) %>%
  dplyr::rename(geneid = Symbol ,left = Start , right = End , chr = CHR) %>% 
  relocate(geneid) %>% 
  write.table("genepos.txt" , sep = "\t" , row.names = F , quote = F)
