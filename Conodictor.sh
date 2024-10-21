#!/bin/bash

#SBATCH -p cicese
#SBATCH --job-name=conodictor
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 06-00:00:00

# 
NPROCS=$SLURM_NPROCS


EXPORT=/LUSTRE/apps/bioinformatica/conodictor/conodictor/
export PATH=$PATH:$EXPORT
 
# Load correct python module conda-2024_py3.11

FASTA=$1

conodictor --out ${FASTA%.*}_dir --cpus $NPROCS $FASTA

# conodictor --out outfolder --cpus $NPROCS --mlen 51 file.fa