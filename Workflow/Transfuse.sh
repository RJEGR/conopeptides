#!/bin/sh
## Directivas
#SBATCH --job-name=MMseqs
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20


EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/transfuse-master/bin/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin
export PATH=$PATH:$EXPORT

export PATH=/LUSTRE/apps/bioinformatica/.local/bin:$PATH
export PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/bin:$PATH

export LD_LIBRARY_PATH=/LUSTRE/apps/bioinformatica/ruby2/ruby-2.2.0/lib:$LD_LIBRARY_PATH


transfuse --assemblies Trinity.fasta,Rnaspades.fasta --left tmp.R1.fq.gz --right tmp.R2.fq.gz --output merged_assemblies.fa --threads $SLURM_NPROCS

transrate --assembly merged_assemblies.fa --left tmp.R1.fq.gz --right tmp.R2.fq.gz --threads $SLURM_NPROCS --output Merged_dir

exit