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

# Bam file most be sorted!

BAM_FILE=$1

mkdir -p S2_TRINITY_DIR

Trinity --genome_guided_bam $BAM_FILE \
         --genome_guided_max_intron 10000 \
         --max_memory 100G --CPU $CPU --output S2_TRINITY_DIR

exit
