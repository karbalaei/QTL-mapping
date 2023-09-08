#!/bin/bash
#$ -cwd
#$ -l mem_free=50G,h_vmem=50G,h_fsize=50G
#$ -N plink_Geno
#$ -o o_plink_geno.txt
#$ -e e_plink_geno.txt
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
module load plink/2
module list

cd /dcs04/lieber/lcolladotor/hydeGoes_LIBD3010/hydeGoes_scSeq_mdd/Tensorqtl

plink2 --bfile MDD_control_sc_maf1 --maf 0.05 --hwe 1e-6 --geno 0.02 --make-bed --out MDD_control_sc_maf_filtered_05_geno
plink2 --bfile MDD_control_sc_maf_filtered_05_geno --indep-pairwise 250 50 0.2 --bad-ld --out MDD_control_sc_maf_info_05_geno
plink2 --bfile MDD_control_sc_maf_filtered_05_geno --extract MDD_control_sc_maf_info_05_geno.prune.in --make-king square --out MDD_control_ibd_out__05_geno

echo "**** Job ends ****"
date