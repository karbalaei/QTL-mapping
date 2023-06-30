#!/bin/bash -l
#SBATCH --job-name=Matrix_eQTL
#SBATCH --output=Matrix_eQTL
#SBATCH --partition=defq
#SBATCH -t 02-00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --export=ALL


ml r/4.0.2
ml r-MatrixEQTL

Rscript Matrix_eQTL.R

