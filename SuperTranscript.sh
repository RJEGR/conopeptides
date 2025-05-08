#!/bin/sh
## Directivas

#SBATCH -p d30
#SBATCH --job-name=SuperDuper
#SBATCH -N 2
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


# 1) INDEX THE REFERENCE

REFERENCE=$1
REF_PREFIX=${REFERENCE%.*}

thread_count=$SLURM_NPROCS

mkdir -p INDEX

# # HISAT


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

WDD=S2_CORSET_FILES_${REF_PREFIX}_DIR

mkdir -p $WDD

cd $WD

EXPORT=/LUSTRE/apps/bioinformatica/corset-1.09-linux64/
export PATH=$PATH:$EXPORT


for f in $(ls *.sorted.bam)
do
basename=${f##*/}
bs="${basename%.sorted.bam}"

if [ ! -f "${bs}.sorted.bam.corset-reads" ]; then

    call="corset -f true -r true ${bs}.sorted.bam"

    echo "Running corset in '${bs}' sample."

    echo $call

    eval $call

else
    echo "'${bs}.sorted.bam.corset-reads' already exists."
    echo "Continue with next sample."
fi

    echo "Continue with next step for corset."

done

corset  -f true -i corset *.corset-reads

mv *corset-reads ../$WDD
mv *txt ../$WDD

# Run Lace


module load conda-2024-2
source activate lace
export PATH=/LUSTRE/apps/bioinformatica/Lace-1.14.1/Lace:$PATH

FASTA=$REFERENCE

OUTDIR=S3_LACE_${FASTA%.*}_DIR

mkdir -p $OUTDIR

Lace ../$FASTA ../$WDD/clusters.txt -o $OUTDIR

mv  $OUTDIR ../

exit