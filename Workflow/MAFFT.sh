#!/bin/bash
#SBATCH --job-name=mafft
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24

echo "Fecha inicio: `date`"
echo "Ejecutandose con $SLURM_JOB_CPUS_PER_NODE"
echo "Numero de nodos: $SLURM_NNODES y CPU por nodo: $SLURM_CPUS_ON_NODE"
echo "CPUs totales= $SLURM_NPROCS"
echo "Los nodos utilizados son: $SLURM_NODELIST"

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/mafft-linux64
export PATH=$PATH:$EXPORT


input=$1

output1=${input%.f*}_default.aln
output2=${input%.f*}_localpair.aln
output3=${input%.f*}_globalpair.aln

mafft --thread $SLURM_NPROCS $input > $output
mafft --thread $SLURM_NPROCS --maxiterate 1000 --localpair  $input > $output2
mafft --thread $SLURM_NPROCS --maxiterate 1000 --globalpair  $input > $output3

exit
