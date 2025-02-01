#!/bin/sh
## Directivas
#SBATCH --job-name=MMseqs
#SBATCH -N 2
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20


EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/mmseqs/bin
export PATH=$PATH:$EXPORT

thread_count=$SLURM_NPROCS

INPUT=$1 # Merged.tpm
ID=0.9
OUTFILE=${INPUT%.*}.${ID}.mmseq.fasta

# s1)
mmseqs createdb $INPUT sequenceDB 

mmseqs clusthash sequenceDB clusthashDB --threads $thread_count --min-seq-id $ID
mmseqs clust sequenceDB clusthashDB clusterDB --threads $thread_count


mmseqs createsubdb clusterDB sequenceDB createsubDB
mmseqs convert2fasta createsubDB $OUTFILE

exit
