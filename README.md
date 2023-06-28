# 1. Tensorqtl 

Tensorqtl algorithm [tensorqtl](https://github.com/broadinstitute/tensorqtl/tree/master) applied on both 21 different cell types and seven different broad cell types resulted from snRNA-seq data. 

![Flowchart](https://github.com/karbalaei/tensorqtl/blob/main/Graph/Flowchart.jpg)


## 1.1. Preparing TSS file:

In Python, using recommended [function](https://github.com/karbalaei/tensorqtl/blob/main/Codes/gtf_tss_code.py) from Tensorqtl tutorial ([pyqtl](https://github.com/broadinstitute/pyqtl)), the TSS file was created which containing chr, start, end and gene_id. 

![TSS file](https://github.com/karbalaei/tensorqtl/blob/main/Graph/TSS_file.jpg)

## 1.2. Preparing BED and Covariate files in R

Using the TSS file, Raw and normalized data of snRNA-Seq, the BED files and the Covariate files were created in [R](https://github.com/karbalaei/tensorqtl/blob/main/Codes/tensorqtl_preparing_files.R). Briefly speaking, I used metadata from the raw file of snRNA-Seq data and merged it with the TSS file. Then, after loading normalized data, I created the sorted BED file. For the covariate file, I used *num.sv* from *SVA* package to calculate how many PCA should I include in my analysis, which for 21 cell types it ranges from 15 to 20.  

Here is an example of *covariate file* from *astro* cell type : 
![Covariate_file](https://github.com/karbalaei/tensorqtl/blob/main/Graph/Covariate_file.jpg)

and Here is an example of *BED* file from *astro* cell type : 
![BED_file](https://github.com/karbalaei/tensorqtl/blob/main/Graph/BED_file.jpg)

## 1.3. Runnig Tensorqtl in BASH

Using a for loop, I run Tensorqtl with three different parameters :

1. Default parameters :

> 	for f in ./*/*bed.gz \
		do \
		echo "Start analyses of $f" \
		python3 -m tensorqtl "MDD_control_sc_maf1" \
		$f \
		${f/_df.bed.gz/} \
		--covariates ${f/df.bed.gz/covariates.txt} \
		--mode cis_nominal \
	done

2. Add a window of 500,000 instead of a defalut value(1Mbp).

> for f in ./*/*bed.gz\
		do \
		echo "Start analyses of $f" \
		python3 -m tensorqtl "MDD_control_sc_maf1" \
		$f \
		${f/_df.bed.gz/_modified} \
		--covariates ${f/df.bed.gz/covariates.txt} \
		--mode cis_nominal --window=500000 \
	done 

3. Add a maf=0.05 instead of a defalut value(0).

> for f in ./*/*bed.gz\
	do \
		echo "Start analyses of $f" \
		python3 -m tensorqtl "MDD_control_sc_maf1" \
		$f \
	${f/_df.bed.gz/_modified_maf} \
    --covariates ${f/df.bed.gz/covariates.txt} \
    --mode cis_nominal --window=500000 --maf_threshold=0.05 \
	done 
	
##1.4.  Loading Parquet file and applying p-value correction :

In the final step, Parquet files from different analyses were loaded in R and after calculation FDR, SNPs with *FDR < 0.05* were selected for further investigation. 

#2. Matrix eQTL

![Flowchart](https://github.com/karbalaei/tensorqtl/blob/main/Graph/Matrix_eqtl.jpg)
