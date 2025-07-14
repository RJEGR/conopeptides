#!/bin/sh
## Directivas
#SBATCH --job-name=kallisto
#SBATCH -N 2
#SBATCH --mem=70GB
#SBATCH --ntasks-per-node=20

# sbatch kallisto.sh Merged_polyA_hisat_SuperDuper.fasta.transdecoder.cds samples.txt

EXPORT=/LUSTRE/bioinformatica_data/genomica_funcional/rgomez/Software/kallisto
export PATH=$PATH:$EXPORT

# kallisto index builds an index from a FASTA formatted file of target sequences. The arguments for the index command are:

REFERENCE=$1 
REF_PREFIX=${REFERENCE%.*} 

Manifest=$2

mkdir -p INDEX

if [ ! -f "INDEX/${REF_PREFIX}_kallisto" ]; then
    kallisto index --index=INDEX/${REF_PREFIX}_kallisto $REFERENCE
else
    echo "index '$REF_PREFIX' already exists for Kallisto."
fi



WD=KALLISTO_${REF_PREFIX}_DIR

mkdir -p $WD

while read -r line; do
    left_file=$(echo "$line" | awk '{print $3}')
    right_file=$(echo "$line" | awk '{print $4}')
    
    bs=${left_file##*/}
    bs="${bs%_R1.fq.gz}"

    output=$WD/${bs}_transcripts_quant

    checkpoint_file="$WD/${bs}.kallisto.chkpt"

    if [ ! -f "$WD/${bs}.kallisto.chkpt" ]; then

    call="kallisto quant --index=INDEX/${REF_PREFIX}_kallisto --threads=$SLURM_NPROCS --output-dir=$output $left_file $right_file"

    eval $call

    touch "$checkpoint_file"

else
    echo "File $bs already exists in $WD directory"

fi


done < "$Manifest"

exit