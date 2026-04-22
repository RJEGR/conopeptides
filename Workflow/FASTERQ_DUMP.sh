#!/bin/bash
#SBATCH --job-name=fasterq-dump
#SBATCH -N 1
#SBATCH --mem=50GB
#SBATCH --ntasks-per-node=20
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=rgomez@uabc.edu.mx

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/sratoolkit.3.1.1-ubuntu64/bin
export PATH=$PATH:$EXPORT

THREADS=$SLURM_NPROCS

for i in $(cat Sra.list); do fasterq-dump --split-files $i --skip-technical -p -e $THREADS;done

exit