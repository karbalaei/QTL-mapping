#!/bin/bash
#$ -cwd
#$ -l mem_free=100G,h_vmem=100G,h_fsize=100G
#$ -N Tensorqtl_maf_nominal
#$ -o o_tensorqtl_genomewide_nominal_maf.txt
#$ -e e_tensorqtl_genomewide_nominal_maf.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"



## List current modules for reproducibility
module load tensorqtl
module list

cd /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/Tensorqtl

for f in ./*/*norm_df.bed.gz;
	do
	echo "location of analysys of $f" 
	ls
	python3 -m tensorqtl "MDD_control_sc_maf1" \
	$f \
	${f/_norm_df.bed.gz/_maf} \
    --covariates ${f/df.bed.gz/covariates.txt} \
    --mode cis_nominal --maf_threshold=0.05

done

echo "**** Job ends ****"
date