#!/bin/sh
## Directivas
#SBATCH --job-name=MMap
#SBATCH -N 1
#SBATCH --mem=70GB
#SBATCH --ntasks-per-node=20

EXPORT=/LUSTRE/apps/bioinformatica/bowtie2/bin/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/samtools-1.3/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/salmon/bin/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/hisat2/
export PATH=$PATH:$EXPORT


# BOWTIE2

# 1) INDEX THE REFERENCE

REFERENCE=$1 
REF_PREFIX=${REFERENCE%.*} 

thread_count=$SLURM_NPROCS

Manifest=$2

mkdir -p INDEX

WD=S1_BOWTIE2_BAM_FILES_${REF_PREFIX}_DIR

mkdir -p $WD
mkdir -p $WD/STATS


# Build Bowtie2 index if not already present
if [ ! -f "INDEX/${REF_PREFIX}.1.bt2" ]; then
    bowtie2-build --threads $thread_count $REFERENCE INDEX/$REF_PREFIX
else
    echo "index '$REF_PREFIX' already exists."
fi

echo "Continue with step two"

WD=S1_BOWTIE2_BAM_FILES_${REF_PREFIX}_DIR

mkdir -p $WD

while read -r line; do
    left_file=$(echo "$line" | awk '{print $3}')
    right_file=$(echo "$line" | awk '{print $4}')
    
    bs=${left_file##*/}
    bs="${bs%_R1.fq.gz}"

    bam_file=${bs}.sorted.bam 

    checkpoint_file="$WD/${bs}.sorted.bam.chkpt"

    if [ ! -f "$WD/${bs}.sorted.bam.chkpt" ]; then

    # Corset and Lace needs report multimap alignments between reads and transcriptome. 
    # This step can be handle using -k INT and --all flag in either, bowtie and hisat. 
    # Nevetheless, Corset authors also encourage to limit the number of reported alignments when large dataset are used.  
    # A finite number of -k 40 instead of --all alignments was used for bowtie2 or hisat2
    
    aligner_params="--no-mixed --no-discordant --gbar 1000 --end-to-end -k 40"
    read_type="-q"
    max_ins_size=800


    call="bowtie2 $aligner_params $read_type -X $max_ins_size \
    -x INDEX/$REF_PREFIX -1 $left_file -2 $right_file -p $thread_count 2> $WD/${bs}.stderr | \
    samtools view -@ $thread_count -F 4 -S -b | samtools sort -@ $thread_count -n -o $WD/$bam_file"

    eval $call

    call="samtools flagstat $WD/$bam_file > $WD/STATS/${bs}.flagstats.txt"

    eval $call

    touch "$checkpoint_file"

else
    echo "File $bs already exists in $WD directory"

fi


done < "$Manifest"

exit


# Salmon

# HISAT


if [ ! -f "INDEX/${REF_PREFIX}.1.ht2" ]; then
hisat2-build -p $thread_count $REFERENCE INDEX/$REF_PREFIX

else
    echo "index for '$REFERENCE' already exists."
fi

WD=S1_HISAT2_SAM_BAM_FILES_${REF_PREFIX}_DIR

mkdir -p $WD

for i in $(ls *_R1.fq.gz)
do
bs="${i%*_R1.fq.gz}"

# Test if the alignment was previously done!

if [ ! -f "$WD/${bs}.sorted.bam.chkpt" ]; then

left_file=${bs}_R1.fq.gz
right_file=${bs}_R2.fq.gz

    hisat2  --phred33 --dta -p $thread_count \
        -x INDEX/$REF_PREFIX -1 $left_file -2 $right_file \
        --rg-id=${bs} --rg SM:${bs} -S $WD/${bs}.sam \
        --summary-file $WD/${bs}.summary.txt --met-file $WD/${bs}.met.txt

    samtools sort -@ $thread_count -o $WD/${bs}.sorted.bam $WD/${bs}.sam

    rm $WD/${bs}.sam

    touch $WD/${bs}.sorted.bam.chkpt

else 
    echo "'$WD/${bs}.sorted.bam' already exists."
    echo "Continue with Corset."
fi

unlink $left_file
unlink $right_file

done
####
# Run Corset
####



exit
