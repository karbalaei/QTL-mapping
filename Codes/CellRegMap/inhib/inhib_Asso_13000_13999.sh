#!/bin/bash -l
#SBATCH --job-name=Asso_inhib_13000_13999
#SBATCH --output=CellRegMap_inhib_Association_13000_13999
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


python inhib_Asso_13000_13999.py
