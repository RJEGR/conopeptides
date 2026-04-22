#!/bin/bash
#SBATCH --job-name=spades
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20
#SBATCH --error=slurm-%j.err

EXPORT=/LUSTRE/apps/bioinformatica/SPAdes-3.15.5-Linux/bin/
export PATH=$PATH:$EXPORT

CPU=$SLURM_NPROCS
MEM=$SLURM_MEM_PER_NODE
DIR=$SLURM_SUBMIT_DIR

# zcat *merged.fq.gz > MERGED.fq

REFERENCE=$1

REF_PREFIX=${REFERENCE%.*}

WD=$DIR/Rnaspades_${REF_PREFIX}_out

mkdir -p $WD

# --merged  MERGED.fq 

rnaspades.py --s1  MERGED.fq --trusted-contigs $REFERENCE -t $CPU -m $MEM -o $WD

exit