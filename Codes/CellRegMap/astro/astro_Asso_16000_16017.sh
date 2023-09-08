#!/bin/bash -l
#SBATCH --job-name=Asso_astro_16000_16017
#SBATCH --output=CellRegMap_astro_Association_16000_16017
#SBATCH --partition=defq
#SBATCH -t 72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=32
#SBATCH --export=ALL
cd /data/fgoes1/RK/CellRegMap/

ml python/3.8.6
ml py-pip
source CellRegMap/bin/activate


python astro_Asso_16000_16017.py
