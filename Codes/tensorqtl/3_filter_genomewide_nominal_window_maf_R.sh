#!/bin/bash
#$ -cwd
#$ -l mem_free=200G,h_vmem=200G,h_fsize=100G
#$ -N filter_R_Win-maf
#$ -o o_filter_genomewide_nominal_window_maf.txt
#$ -e e_filter_genomewide_nominal_window_maf.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${JOB_ID}"
echo "Job name: ${JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SGE_TASK_ID}"

## Load the R module (absent since the JHPCE upgrade to CentOS v7)
module load conda_R

## List current modules for reproducibility
module list

## Edit with your job command
Rscript filter_genomewide_nominal_window_maf.R

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/