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

rnaspades.py $left_reads $right_read \
  -k 21,25,49,73 \
  --pe1-fr -t $CPU -m $MEM \
  --checkpoints last -o $DIR/Rnaspades_out

exit