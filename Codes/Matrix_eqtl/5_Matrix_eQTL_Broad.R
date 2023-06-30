library(MatrixEQTL)
library(tidyverse)
library(data.table)

dir = ("/data/fgoes1/RK/tensorqtl/MDD_eQTL/")

setwd(dir)

MDD_control_sc_maf_original <-  read.table("MDD_control_sc_maf1.txt", row.names = 1, sep = "\t", header = T)


cell_type_names <- c("astro","endo","excit", "inhib" ,"micro" ,"oligo" ,"opc")

useModel = modelLINEAR



for (  i in 1 : length(cell_type_names)) {

  setwd(paste0(dir , "Broad/", cell_type_names[i] , "_data"))
  
   name_GE = paste0(cell_type_names[i] ,"_GE.txt")
  GE = read.table(name_GE, sep = "\t", header = T)
  
  Sample_Ids = colnames(GE)[-1]

    MDD_control_sc_maf_original %>% dplyr::select(all_of(Sample_Ids)) %>%
    fwrite(paste0(cell_type_names[i] ,"_MDD_control_sc_maf.txt"), sep = "\t" , quote = F , row.names = T)
    
    setwd(dir)

SNP_file_name =paste0("Broad/",cell_type_names[i] ,"_data/", cell_type_names[i] ,"_MDD_control_sc_maf.txt")

expression_file_name=paste0("Broad/",cell_type_names[i] ,"_data/", cell_type_names[i] ,"_GE.txt")

covariates_file_name=paste0("Broad/",cell_type_names[i] ,"_data/", cell_type_names[i] ,"_covariates.txt")



snps_location_file_name = ("snpsloc.txt")
gene_location_file_name = ("genepos.txt")

# Output file name

output_file_name_cis = paste0("Broad/",cell_type_names[i] ,"_data/", cell_type_names[i] ,"_results_cis")

output_file_name_tra = paste0("Broad/",cell_type_names[i] ,"_data/", cell_type_names[i] ,"_results_trans")

output_file_name=paste0("Broad/",cell_type_names[i] ,"_data/", cell_type_names[i] ,"_results")


pvOutputThreshold_cis = 1
pvOutputThreshold_tra = 0
errorCovariance = numeric()
cisDist = 1e6;

snps = SlicedData$new()
snps$fileDelimiter = "\t"
snps$fileOmitCharacters = "NA"
snps$fileSkipRows = 1
snps$fileSkipColumns = 1
snps$fileSliceSize = 200000
snps$LoadFile( SNP_file_name )

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);


## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
cvrt$LoadFile(covariates_file_name)
## Run the analysis
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
snps = snps,
gene = gene,
cvrt = cvrt,
output_file_name     = output_file_name_tra,
pvOutputThreshold     = pvOutputThreshold_tra,
useModel = useModel,
errorCovariance = errorCovariance,
verbose = TRUE,
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos,
genepos = genepos,
cisDist = cisDist,
pvalue.hist = FALSE ,
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE)

# 0.5 Mbp distance

# Output file name

output_file_name_cis = paste0("Broad/",cell_type_names[i] ,"_data/", cell_type_names[i] ,"_results_cis_2")

output_file_name_tra = paste0("Broad/",cell_type_names[i] ,"_data/", cell_type_names[i] ,"_results_trans_2")

output_file_name=paste0("Broad/",cell_type_names[i] ,"_data/", cell_type_names[i] ,"_results_2")


pvOutputThreshold_cis = 1
pvOutputThreshold_tra = 0
errorCovariance = numeric()
cisDist = 5e5;

me_2 = Matrix_eQTL_main(
snps = snps,
gene = gene,
cvrt = cvrt,
output_file_name     = output_file_name_tra,
pvOutputThreshold     = pvOutputThreshold_tra,
useModel = useModel,
errorCovariance = errorCovariance,
verbose = TRUE,
output_file_name.cis = output_file_name_cis,
pvOutputThreshold.cis = pvOutputThreshold_cis,
snpspos = snpspos,
genepos = genepos,
cisDist = cisDist,
pvalue.hist = FALSE ,
min.pv.by.genesnp = FALSE,
noFDRsaveMemory = FALSE)


}