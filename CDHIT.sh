#!/bin/bash
#SBATCH --job-name=CDHIT
#SBATCH -N 1
#SBATCH --mem=100GB
#SBATCH --ntasks-per-node=20

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/apps/cdhit/
export PATH=$PATH:$EXPORT

identity=$3  # ranges between 0 to 1, 0.95 means 95% identity threshold
word_length=11  # CD-HIT word length parameter
threads=20  # Number of threads for parallel processing

cat $1 $2 > Merged.tmp

REFERENCE=Merged.tmp


OUTFILE=${REFERENCE%.*}.${identity}.cdhit.fasta

# Run CD-HIT-EST
cd-hit-est -i "${REFERENCE}" -o "${OUTFILE}" -c "${identity}" -n "${word_length}" -T "${threads}" -M 20000

echo "CD-HIT-EST completed. Output saved to ${OUTFILE}"

rm *.tmp

exit