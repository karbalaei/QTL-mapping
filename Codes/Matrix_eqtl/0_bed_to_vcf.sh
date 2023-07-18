#!/bin/bash
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G,h_fsize=10G
#$ -N plink
#$ -o o_plink_bed_to_vcf.txt
#$ -e e_plink_bed_to_vcf.txt
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
module load plink/1.90b6.6
module list

cd /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/Tensorqtl/MDD_eQTL

plink --bfile MDD_control_sc_maf1 --recode vcf --out MDD_control_sc_maf1

echo "**** Job ends ****"
date