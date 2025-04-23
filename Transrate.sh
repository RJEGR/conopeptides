#!/bin/bash
#SBATCH --job-name=trate
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 6-00:00:00


TRANSRATE=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin
export PATH=$PATH:$TRANSRATE

export PATH=/LUSTRE/apps/bioinformatica/.local/bin:$PATH
export PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin:$PATH
export LD_LIBRARY_PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/lib:$LD_LIBRARY_PATH

fasta=$1
left=$2
right=$3

output=${fasta%.*}_transrate_dir

transrate --assembly $fasta --left $left --right $right --threads $SLURM_NPROCS --output $output

exit