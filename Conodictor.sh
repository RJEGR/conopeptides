#!/bin/bash

#SBATCH -p cicese
#SBATCH --job-name=conopred
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=24
#SBATCH -t 06-00:00:00

# The program requires only the assembled transcriptome or raw reads file, 
# in either DNA or amino acid format. ConoDictor 2 automatically recognizes the alphabet used.
NPROCS=$SLURM_NPROCS


module load conda-2024
source activate conodictor

FASTA=$1

conodictor --out ${FASTA%.*}_dir --cpus $NPROCS $FASTA

# Run conosorter

export PATH=/LUSTRE/apps/bioinformatica/hmmer-3.3.2/bin:$PATH
export PATH=/LUSTRE/apps/bioinformatica/ConoSorter_v1.1:$PATH

ConoSorter -d $FASTA

exit

# conodictor --out outfolder --cpus $NPROCS --mlen 51 file.fa
# Fasta file must contain only ATCG characters
EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/
export PATH=$PATH:$EXPORT

grep -v "^>" $FASTA  | grep -c "NNN" # 

cat $FASTA | seqkit grep -v -s -i -p NNN > ${FASTA%.*}_no_redundant_NNN.fasta


cat $FASTA | seqkit grep -s -i -p NNN | seqkit fx2tab -l 