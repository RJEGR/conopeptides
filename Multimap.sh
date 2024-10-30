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

mkdir -p INDEX

# Build Bowtie2 index if not already present
if [ ! -f "INDEX/${REF_PREFIX}.1.bt2" ]; then
    bowtie2-build --threads $thread_count $REFERENCE INDEX/$REF_PREFIX
else
    echo "index '$REF_PREFIX' already exists."
fi

echo "Continue with step two"

WD=S1_BOWTIE2_BAM_FILES_${REF_PREFIX}_DIR

mkdir -p $WD

for i in $(ls *_R1.fq.gz)
do
bs="${i%*_R1.fq.gz}"

aligner_params="--no-mixed --no-discordant --gbar 1000 --end-to-end -k 200"
read_type="-q"
max_ins_size=800


if [ ! -f "$WD/${bs}.sorted.bam.chkpt" ]; then
    
    echo "Aligning reads back to reference"


    left_file=${bs}_R1.fq.gz
    right_file=${bs}_R2.fq.gz  

    bam_file=${bs}.sorted.bam 
    
    bowtie2 --met-file $met_file $aligner_params $read_type -X $max_ins_size \
    -x INDEX/$REF_PREFIX -1 $left_file -2 $right_file -p $thread_count 2> $WD/${bs}.stderr | \
    samtools view -@ $thread_count -F 4 -S -b | samtools sort -@ $thread_count -n -o $WD/$bam_file

    touch $WD/${bs}.sorted.bam.chkpt

else
    echo "File $bam_file already exists in $WD directory"

fi

done

# Salmon

if [ ! -f "INDEX/SALMON_${REF_PREFIX}/ctable.bin" ]; then
   salmon index --index INDEX/SALMON_${REF_PREFIX} --transcripts $REFERENCE --threads $thread_count 

else
    echo "index '$REF_PREFIX' already exists for Salmon."
fi

WD=S1_SALMON_BAM_FILES_${REF_PREFIX}_DIR

mkdir -p $WD

for i in $(ls *_R1.fq.gz)
do
bs="${i%*_R1.fq.gz}"

if [ ! -f "$WD/${bs}.salmon.chkpt" ]; then
    
    echo "Aligning reads back to reference"


    left_file=${bs}_R1.fq.gz
    right_file=${bs}_R2.fq.gz  

    salmon_file=${bs}.out 
    
    call="salmon quant --index INDEX/SALMON_${REF_PREFIX} --libType A --dumpEq --hardFilter --skipQuant -1 $left_file -2 $right_file --threads $thread_count --output $WD"

    echo $call

    eval $call


    touch $WD/${bs}.salmon.chkpt

else
    echo "File $salmon_file already exists in $WD directory"

fi

done

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
