#!/bin/bash -l
#SBATCH --job-name=Asso_micro_10000_10999
#SBATCH --output=CellRegMap_micro_Association_10000_10999
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


python micro_Asso_10000_10999.py
