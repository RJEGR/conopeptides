#!/bin/bash

#SBATCH -p cicese
#SBATCH --job-name=Lace
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 06-00:00:00

# 

module load conda-2024-2
source activate lace
export PATH=/LUSTRE/apps/bioinformatica/Lace-1.14.1/Lace:$PATH

FASTA=$1

OUTDIR=LACE_${FASTA%.*}_DIR

mkdir -p $OUTDIR

Lace --cores $SLURM_NPROCS $FASTA clusters.txt -o $OUTDIR

exit
