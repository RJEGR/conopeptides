#!/bin/sh
## Directivas
#SBATCH --job-name=Decoder
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20

# Trinotate

FILE=$1 #transcriptome.fasta
bs=`basename ${FILE%.f*}`

EXPORT=/LUSTRE/apps/bioinformatica/TransDecoder-v5.7.0/
export PATH=$PATH:$EXPORT

output_dir=${bs}.transdecoder_dir

#TransDecoder.LongOrfs  --complete_orfs_only                   yields only complete ORFs (peps start with Met (M), end with stop (*))
# TransDecoder.LongOrfs -m minimum protein length (default: 100)

TransDecoder.LongOrfs -t $FILE --output_dir $output_dir --complete_orfs_only -m 50

ln -s ${output_dir}/longest_orfs.pep ${bs}_longest_orfs.pep

# count number of unique contigs with orfs
grep "^>" ${bs}_longest_orfs.pep | awk '{gsub(/\.p[0-9]+$/,"",$1); print $1}' | uniq 

exit

