#create Matrix of genotyping fro vcf file
# ml r4.2.0
#R
library(vcfR)
#library(tidyverse)
library(data.table)
#library(readr)
library(dplyr)
library(stringr)

setwd("/data/fgoes1/RK/tensorqtl")
vcf <-  read.vcfR( "MDD_control_sc_maf1.vcf", verbose = FALSE )

gt <-  extract.gt(vcf, element = 'GT', as.numeric = F)


gt_df <-  gt %>% as.data.frame() %>%
  mutate_at(vars(), funs(((function(x) {
  if (x=="0/0")
    return(0)
  else if (x=="1/1")
    return(2)
  else if (x=="0/1")
      return(1)
  else
    return(1)
})(.))))


gt %>% as.data.frame() %>%
  mutate_all(funs(str_replace(., "0/0", "0"))) %>%
  mutate_all(funs(str_replace(., "0/1", "1"))) %>%
  mutate_all(funs(str_replace(., "1/0", "1"))) %>%
  mutate_all(funs(str_replace(., "1/1", "2"))) %>%
  mutate_all(gt_df2, function(x) as.numeric(as.character(x))) %>%
  setNames(nams) %>%
  fwrite("MDD_control_sc_maf1_df.txt", sep = "\t" , quote = F , row.names = T)

  


df2 <- mutate_all(gt_df2, function(x) as.numeric(as.character(x)))




nams <-  colnames(gt)

pattens <- c("R01C02_","R01C01_","A_","R07C01_","R02C01_","R05C01_","R03C01_","R08C01_","R04C01_","R06C01_")

for( i in 1:length(pattens)){
  nams <-  str_remove_all(nams , pattens[i])
}


gt_df %>%setNames(nams) %>%
  fwrite("MDD_control_sc_maf1.txt", sep = "\t" , quote = F , row.names = T)
