#!/bin/bash

#SBATCH -p cicese
#SBATCH --job-name=conodictor
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 06-00:00:00

# 
NPROCS=$SLURM_NPROCS


module load conda-2024
source activate conodictor

FASTA=$1

conodictor --out ${FASTA%.*}_dir --cpus $NPROCS $FASTA

exit

# conodictor --out outfolder --cpus $NPROCS --mlen 51 file.fa