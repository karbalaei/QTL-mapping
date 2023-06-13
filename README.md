# Tensorqtl results
Here you may find the Flowchart of running [tensorqtl](https://github.com/broadinstitute/tensorqtl/tree/master) on snRNA-seq data :

![Flowchart](https://github.com/karbalaei/tensorqtl/blob/main/Graph/Flowchart.jpg)


## Preparing TSS file:

In python using recommended [function](https://github.com/karbalaei/tensorqtl/blob/main/Codes/gtf_tss_code.py) from tensorqtl tuturial ([pyqtl](https://github.com/broadinstitute/pyqtl)), the TSS file was created which containing chr, start, end and gene_id. 

![TSS file](https://github.com/karbalaei/tensorqtl/blob/main/Graph/TSS_file.jpg)

## Preparing BED and Covariate files in R

Using TSS file, Raw and normalized data of snRNA-Seq, the BED files and Covariate files were created in [R](https://github.com/karbalaei/tensorqtl/blob/main/Codes/tensorqtl_preparing_files.R). Breafly speaking, I used metadata from raw file of snRNA-Seq data and combin it with TSS file. Then, after loading normalized data, I created the sorted BED file. For covariate file, I used *num.sv* from *SVA* package to calculae how many PCA should i include in my  analysis, which for 21 cell types it ranges from 15 to 20.  

Here is an example of *covariate file* from *astro* cell type : 
![Covariate_file](https://github.com/karbalaei/tensorqtl/blob/main/Graph/Covariate_file.jpg)

and Here is an example of *BED* file from *astro* cell type : 
![BED_file](https://github.com/karbalaei/tensorqtl/blob/main/Graph/BED_file.jpg)

## Run Tensorqtl in BASH

using a for loop, I run tensorqtl with three different parameters :

1. Default parameters :

> 	for f in ./*/*bed.gz;
	do
	echo "location of analyses of $f" 
	python3 -m tensorqtl "MDD_control_sc_maf1" \
	$f \
	${f/_df.bed.gz/} \
    --covariates ${f/df.bed.gz/covariates.txt} \
    --mode cis_nominal

done

2. Add a window of 500,000 instead of  defalut value(1Mbp).

3. Add maf of 0.05 instead of  defalut value(0).