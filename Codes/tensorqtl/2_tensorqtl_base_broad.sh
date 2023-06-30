#!/bin/bash
#$ -cwd
#$ -l mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -N Tensorqtl_Base_nominal
#$ -o o_tensorqtl_genomewide_nominal_base_broad.txt
#$ -e e_tensorqtl_genomewide_nominal_base_broad.txt
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

for f in ./Broad_cell_type/*/*.bed.gz;
	do
	echo "Start analyses of $f" 
	python3 -m tensorqtl "MDD_control_sc_maf1" \
	$f \
	${f/_broad.bed.gz/_broad} \
    --covariates ${f/bed.gz/covariates.txt} \
    --mode cis_nominal
    
done

echo "**** Job ends ****"
date