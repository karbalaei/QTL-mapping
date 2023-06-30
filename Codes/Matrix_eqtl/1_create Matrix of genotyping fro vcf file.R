#create Matrix of genotyping fro vcf file
# ml r4.2.0
#R
library(vcfR)
library(tidyverse)
library(data.table)
library(readr)


setwd(datafgoes1RKtensorqtl)
vcf - read.vcfR( MDD_control_sc_maf1.vcf, verbose = FALSE )

gt - extract.gt(vcf, element = 'GT', as.numeric = TRUE)


gt_df = gt %>% as.data.frame()

nams <-  colnames(gt_df)

pattens <- c("R01C02_","R01C01_","A_","R07C01_","R02C01_","R05C01_","R03C01_","R08C01_","R04C01_","R06C01_")

for( i in 1:length(pattens)){
  nams <-  str_remove_all(nams , pattens[i])
}


gt_df %>%setNames(nams) %>%
  fwrite("MDD_control_sc_maf1.txt", sep = "\t" , quote = F , row.names = T)
