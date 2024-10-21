#!/bin/sh
## Directivas
#SBATCH --job-name=corset
#SBATCH -N 1
#SBATCH --mem=70GB
#SBATCH --ntasks-per-node=20

EXPORT=/LUSTRE/apps/bioinformatica/bowtie2/bin/
export PATH=$PATH:$EXPORT

EXPORT=/LUSTRE/apps/bioinformatica/RSEM/bin/samtools-1.3/
export PATH=$PATH:$EXPORT


EXPORT=/LUSTRE/apps/bioinformatica/salmon/bin/
export PATH=$PATH:$EXPORT

# BOWTIE2

# 1) INDEX THE REFERENCE

REFERENCE=$1 
REF_PREFIX=${REFERENCE%.*} 

thread_count=$SLURM_NPROCS

mkdir -p INDEX

# Build Bowtie2 index if not already present
if [ ! -f "INDEX/$REF_PREFIX.1.bt2" ]; then
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
    
    call="bowtie2 --met-file $met_file $aligner_params $read_type -X $max_ins_size \
    -x INDEX/$REF_PREFIX -1 $left_file -2 $right_file -p $thread_count 2> ${bs}.stderr | \
    samtools view -@ $thread_count -F 4 -S -b | samtools sort -@ $thread_count -n -o S1_BOWTIE2_BAM_FILES/$bam_file"

    echo $call

    eval $call


    touch $WD/${bs}.sorted.bam.chkpt

else
    echo "File $bam_file already exists in S1_BOWTIE2_BAM_FILES directory"

fi

#unlink $left_file
#unlink $right_file

done

# Salmon

if [ ! -f "INDEX/$REF_PREFIX.1.bt2" ]; then
   salmon index --index INDEX/SALMON_${REF_PREFIX} --transcripts $REFERENCE

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
    
    call="salmon quant --index INDEX/SALMON_${REF_PREFIX} --libType A --dumpEq --hardFilter --skipQuant -1 $left_file -2 $right_file --output $WD/$salmon_file"

    echo $call

    eval $call


    touch $WD/${bs}.salmon.chkpt

else
    echo "File $salmon_file already exists in $WD directory"

fi

unlink $left_file
unlink $right_file

done

exit
