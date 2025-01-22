#!/bin/bash
#SBATCH --job-name=Trinity
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20
#SBATCH -t 6-00:00:00

module load trinityrnaseq-v2.15.1

CPU=$SLURM_NPROCS
MEM=$SLURM_MEM_PER_NODE
DIR=$SLURM_SUBMIT_DIR

Trinity --genome_guided_bam rnaseq.coordSorted.bam \
         --genome_guided_max_intron 10000 \
         --max_memory ${MEM}G --CPU $CPU 

exit
