#!/bin/sh
## Directivas
#SBATCH --job-name=Salmon
#SBATCH -N 1
#SBATCH --mem=70GB
#SBATCH --ntasks-per-node=20

# sbatch SALMON.sh Merged_polyA_hisat_SuperDuper.fasta.transdecoder.cds samples.txt


EXPORT=/LUSTRE/apps/bioinformatica/salmon/bin/
export PATH=$PATH:$EXPORT


REFERENCE=$1
REF_PREFIX=${REFERENCE%.*}

Manifest=$2

thread_count=$SLURM_NPROCS

mkdir -p INDEX

if [ ! -f "INDEX/SALMON_${REF_PREFIX}/ctable.bin" ]; then
   salmon index --index INDEX/SALMON_${REF_PREFIX} --transcripts $REFERENCE --threads $thread_count
   
   ln -s $(realpath $REFERENCE) INDEX/SALMON_${REF_PREFIX}/$(basename $REFERENCE)

else
    echo "index '$REF_PREFIX' already exists for Salmon."
fi

WD=SALMON_${REF_PREFIX}_DIR

mkdir -p $WD

while read -r line; do
    left_file=$(echo "$line" | awk '{print $3}')
    right_file=$(echo "$line" | awk '{print $4}')
    bs=${left_file##*/}
    bs="${bs%_R1.fq.gz}"

    output=$WD/${bs}_transcripts_quant

    checkpoint_file="$WD/${bs}.salmon.chkpt"

    if [ ! -f "$WD/${bs}.salmon.chkpt" ]; then

    call="salmon quant --index INDEX/SALMON_${REF_PREFIX} --libType A -1 $left_file -2 $right_file --threads $thread_count  --output $output"
    #echo $call

    eval $call

    touch "$checkpoint_file"

else
    echo "File $bs already exists in $WD directory"

fi


done < "$Manifest"

exit